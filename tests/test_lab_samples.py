"""Sample identity is sufficient for genotyping; visit metadata stays optional."""

import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
RT = ROOT / "runtime"
sys.path.insert(0, str(RT / "src"))

from simplseq.metadata import enrich_rows_with_metadata, inspect_metadata
from simplseq.samplesheet import write_samples_csv


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


analysis = load_module("lab_analysis", RT / "workflow/scripts/prepare_analysis_table.py")
cdhit = load_module("lab_cdhit", RT / "workflow/scripts/cdhit_clustering.py")
gui = load_module("lab_gui", RT / "gui/flask_app.py")


class SampleMetadataTests(unittest.TestCase):
    def test_illumina_numeric_sample_ids_match_metadata_without_false_participants(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for read in (1, 2):
                (root / f"7551_S1_L001_R{read}_001.fastq.gz").touch()
            metadata = root / "metadata.csv"
            metadata.write_text("sample_id,strain\n7551,3D7\n")
            output = root / "samples.csv"
            self.assertEqual(write_samples_csv(root, output, metadata_path=metadata), (1, []))
            row = pd.read_csv(output, keep_default_na=False).iloc[0]
            self.assertEqual((str(row.sample_id), str(row.biological_sample_id)), ("7551", "7551"))
            self.assertEqual(row.metadata_match_status, "matched")
            self.assertEqual((row.participant_id, row.collection_date), ("", ""))

    def test_sample_only_metadata_preserves_lab_attributes_without_inventing_visits(self):
        with tempfile.TemporaryDirectory() as tmp:
            metadata = Path(tmp) / "metadata.csv"
            metadata.write_text("sample_id,strain,concentration\nLab-A,3D7,10\n")
            catalog = inspect_metadata(metadata)
            self.assertFalse([issue for issue in catalog.issues if issue.severity == "error"])
            self.assertEqual(catalog.summary()["sample_records"], 1)
            self.assertEqual(catalog.summary()["longitudinal_records"], 0)
            row = enrich_rows_with_metadata([{"sample_id": "Lab-A"}], metadata)[0]
            self.assertEqual(row["metadata_match_status"], "matched")
            self.assertEqual(json.loads(row["metadata_attributes"]), {"strain": "3D7", "concentration": "10"})
            self.assertFalse(row.get("participant_id"))
            self.assertFalse(row.get("collection_date"))

    def test_full_fastq_filename_links_to_its_parsed_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fastq = "mpg_L28896-1_LIB-1-LabA_Rep1_S1_R1_001.fastq.gz"
            (root / fastq).touch()
            (root / fastq.replace("_R1_", "_R2_")).touch()
            metadata = root / "metadata.csv"
            metadata.write_text(f"filename,sample_type,strain\n{fastq},sample,3D7\n")
            output = root / "samples.csv"
            count, duplicates = write_samples_csv(root, output, metadata_path=metadata)
            self.assertEqual((count, duplicates), (1, []))
            row = pd.read_csv(output, keep_default_na=False).iloc[0]
            self.assertEqual(row.metadata_match_status, "matched")
            self.assertEqual(row.collection_date, "")
            self.assertEqual(json.loads(row.metadata_attributes)["strain"], "3D7")

    def test_manual_column_mapping_and_sample_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            metadata = Path(tmp) / "metadata.csv"
            metadata.write_text("tube,kind,specimen,repeat\nLabA,positive_control,C1,2\n")
            row = enrich_rows_with_metadata([{"sample_id": "LabA"}], metadata, column_overrides={
                "sample_id": "tube", "sample_type": "kind", "biological_sample_id": "specimen", "replicate": "repeat",
            })[0]
            self.assertEqual((row["sample_type"], row["biological_sample_id"], row["replicate"]), ("positive", "C1", "2"))

    def test_sample_ids_keep_punctuation_and_do_not_borrow_other_visits(self):
        with tempfile.TemporaryDirectory() as tmp:
            metadata = Path(tmp) / "metadata.csv"
            metadata.write_text("sample_id,participant_id,collection_date\nA-1,P1,2022-03-01\n")
            row = enrich_rows_with_metadata([{"sample_id": "A1", "participant_id": "P1"}], metadata)[0]
            self.assertEqual(row["metadata_match_status"], "missing")
            self.assertFalse(row.get("collection_date"))

    def test_conflicting_sample_rows_are_not_guessed(self):
        with tempfile.TemporaryDirectory() as tmp:
            metadata = Path(tmp) / "metadata.csv"
            metadata.write_text("sample_id,strain\nLabA,3D7\nLabA,Dd2\n")
            row = enrich_rows_with_metadata([{"sample_id": "LabA"}], metadata)[0]
            self.assertEqual(row["metadata_match_status"], "ambiguous")
            self.assertEqual(row["metadata_attributes"], "")

    def test_sample_match_supplies_real_participant_and_exact_date(self):
        with tempfile.TemporaryDirectory() as tmp:
            metadata = Path(tmp) / "metadata.csv"
            metadata.write_text("sample_id,participant_id,collection_date,pcr\nLabA,P1,2022-03-12,Pos\n")
            row = enrich_rows_with_metadata([{"sample_id": "LabA"}], metadata)[0]
            self.assertEqual((row["participant_id"], row["collection_date"]), ("P1", "2022-03-12"))

    def test_sample_metadata_month_conflicts_are_not_guessed(self):
        with tempfile.TemporaryDirectory() as tmp:
            metadata = Path(tmp) / "metadata.csv"
            metadata.write_text("sample_id,participant_id,visit_month\nLabA,P1,May\nLabA,P1,June\n")
            row = enrich_rows_with_metadata([{"sample_id": "LabA"}], metadata)[0]
            self.assertEqual(row["metadata_match_status"], "ambiguous")


class LabAnalysisTests(unittest.TestCase):
    def make_run(self, root):
        run = root / "run"
        (run / "run_dada2").mkdir(parents=True)
        samples = pd.DataFrame([
            ["LabA_Rep1", "LabA", "sample", "1"],
            ["LabA_Rep2", "LabA", "sample", "2"],
            ["LabB", "LabB", "sample", ""],
            ["NTC", "NTC", "negative", ""],
        ], columns=["sample_id", "biological_sample_id", "sample_type", "replicate"])
        samples.to_csv(run / "samples.csv", index=False)
        pd.DataFrame([
            ["LabA_Rep1", 1000, 1000, 50], ["LabA_Rep2", 1000, 1000, 0],
            ["LabB", 1000, 1000, 0], ["NTC", 1000, 1000, 1000],
        ], columns=["sample", "AMA,.", "KELT,.", "AMA,31A"]).to_csv(run / "run_dada2/seqtab_cigar.tsv", sep="\t", index=False)
        return run

    def test_controls_excluded_but_undated_samples_remain_eligible(self):
        samples = pd.DataFrame({"sample_id": ["LabA", "LabB", "NTC"], "sample_type": ["sample", "sample", "negative"]})
        kept, audit = analysis.biological_metadata(samples)
        self.assertEqual(kept.sample_id.tolist(), ["LabA", "LabB"])
        self.assertEqual(audit["controls_excluded"], 1)
        self.assertEqual(audit["unassigned_or_undated_excluded"], 0)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samples.csv"
            samples.to_csv(path, index=False)
            ids, _ = cdhit.biological_sample_ids(path)
            self.assertEqual(ids, {"LabA", "LabB"})

    def test_gui_shared_input_and_downstream_readiness_without_dates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = self.make_run(root)
            with patch.dict(os.environ, {"SIMPLSEQ_CONFIG_DIR": str(root / "config")}):
                client = gui.create_app(RT, root).test_client()
                response = client.post("/api/analysis-table/build", json={"outdir": str(run), "mode": "primary"})
                self.assertEqual(response.status_code, 200, response.get_json())
                summary = response.get_json()["summary"]
                self.assertEqual(summary["output"]["biological_samples"], 2)
                self.assertEqual(summary["output"]["participant_visits"], 0)
                self.assertEqual(summary["output"]["participants"], 0)
                self.assertEqual(summary["longitudinal"]["calendar_visits"], 0)
                calls = pd.read_csv(run / "analysis_input/primary_alleles.tsv", sep="\t")
                self.assertEqual(set(calls.sample_id), {"LabA", "LabB"})
                self.assertNotIn("31A", set(calls.allele))
                readiness = gui.analysis_readiness(run)
                self.assertTrue(readiness["dcifer_ready"])
                self.assertFalse(readiness["dinemites_ready"])
                response = client.post("/api/dinemites/run", json={"outdir": str(run)})
                self.assertEqual(response.status_code, 400)
                self.assertIn("participant ID and collection dates", response.get_json()["error"])

    def test_sample_only_scan_is_ready_and_metadata_template_is_generic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "LabA_R1.fastq.gz").touch()
            (root / "LabA_R2.fastq.gz").touch()
            metadata = root / "metadata.csv"
            metadata.write_text("sample_id,strain\nLabA,3D7\n")
            client = gui.create_app(RT, root).test_client()
            inspected = client.post("/api/metadata/inspect", json={"path": str(metadata)}).get_json()
            self.assertTrue(inspected["ok"])
            response = client.post("/api/scan", json={"fastq_dir": str(root), "metadata_path": str(metadata), "write_samples": False})
            self.assertEqual(response.status_code, 200, response.get_json())
            row = response.get_json()["scan_rows"][0]
            self.assertEqual((row["review_status"], row["metadata_match_status"]), ("ready", "matched"))
            self.assertEqual(row["collection_date"], "")
            self.assertTrue(client.get("/api/metadata/template").data.startswith(b"sample_id,"))

    def test_partial_visit_identity_cannot_split_one_biological_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = self.make_run(root)
            samples = pd.read_csv(run / "samples.csv", keep_default_na=False)
            samples["participant_id"] = ["P1", "", "", ""]
            samples["collection_date"] = ["2022-03-12", "", "", ""]
            samples.to_csv(run / "samples.csv", index=False)
            client = gui.create_app(RT, root).test_client()
            response = client.post("/api/analysis-table/build", json={"outdir": str(run)})
            self.assertEqual(response.status_code, 400)
            self.assertIn("conflicting participants or collection dates", response.get_json()["error"])

    def test_longitudinal_alleles_are_unique_across_loci(self):
        from simplseq.longitudinal import build_dinemites_inputs, build_visit_calendar
        samples = pd.DataFrame({"sample_id": ["a", "b"], "participant_id": ["P1", "P1"], "collection_date": ["2022-03-01", "2022-04-01"]})
        calls = pd.DataFrame([["P1", date, locus, "."] for date in samples.collection_date for locus in ("AMA", "KELT")], columns=["participant_id", "collection_date", "locus", "allele"])
        calendar, _, _ = build_visit_calendar(calls, samples)
        observed, _ = build_dinemites_inputs(calls, calendar)
        self.assertEqual(set(observed.allele), {"AMA:.", "KELT:."})


@unittest.skipUnless(shutil.which("Rscript"), "Workflow R is required")
class LabRIntegrationTests(unittest.TestCase):
    def test_numeric_ids_keep_leading_zeros_through_dcifer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ids = ["007001", "007002", "007003"]
            pd.DataFrame({"sample_id": ids}).to_csv(root / "samples.csv", index=False)
            pd.DataFrame({"sample": ids, "AMA,.": [1000, 0, 600], "AMA,31A": [0, 1000, 400],
                          "KELT,.": [500, 1000, 0], "KELT,42A": [500, 0, 1000]}).to_csv(root / "counts.tsv", sep="\t", index=False)
            bridge = subprocess.run(["Rscript", str(RT / "workflow/scripts/simplseq_to_dcifer.R"), "--cigar", str(root / "counts.tsv"),
                                     "--samples", str(root / "samples.csv"), "--out", str(root / "dcifer.tsv")],
                                    cwd=root, text=True, capture_output=True, timeout=90)
            self.assertEqual(bridge.returncode, 0, bridge.stdout + bridge.stderr)
            rows = pd.read_csv(root / "dcifer.tsv", sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual(set(rows.sample_id), set(ids))
            estimated = subprocess.run(["Rscript", str(RT / "workflow/scripts/run_dcifer.R"), "--input", str(root / "dcifer.tsv"),
                                        "--outdir", str(root / "dcifer"), "--ibd_grid_nr", "100"],
                                       cwd=root, text=True, capture_output=True, timeout=90)
            self.assertEqual(estimated.returncode, 0, estimated.stdout + estimated.stderr)

    def test_reference_mapping_accepts_undated_samples_and_excludes_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with (RT / "reference/amplicons_noprimers.fasta").open() as handle:
                ref = next(record for record in SeqIO.parse(handle, "fasta") if record.id == "AMA")
            pd.DataFrame({"sample": ["LabA", "LabB", "NTC"], str(ref.seq): [1000, 1000, 5000]}).to_csv(root / "counts.tsv", sep="\t", index=False)
            pd.DataFrame({"sample_id": ["LabA", "LabB", "NTC"], "sample_type": ["sample", "sample", "negative"]}).to_csv(root / "samples.csv", index=False)
            result = subprocess.run(["Rscript", str(RT / "workflow/scripts/postProc_dada2.R"), "-s", str(root / "counts.tsv"), "--samples", str(root / "samples.csv"), "-ref", str(RT / "reference/amplicons_noprimers.fasta"), "--snv_filter", str(RT / "reference/snv_filters.txt"), "-o", str(root / "mapped.tsv"), "--fasta"], cwd=root, text=True, capture_output=True, timeout=90)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            row = pd.read_csv(root / "mapped.tsv", sep="\t").iloc[0]
            self.assertEqual((row.total_reads, row.total_samples, row.qc_total_reads_all_rows), (2000, 2, 7000))

    def test_dcifer_bridge_keeps_distinct_undated_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame({"sample_id": ["LabA", "LabB", "NTC"], "sample_type": ["sample", "sample", "negative"]}).to_csv(root / "samples.csv", index=False)
            pd.DataFrame({"sample": ["LabA", "LabB", "NTC"], "AMA,.": [1000, 1000, 5000], "KELT,.": [1000, 1000, 5000]}).to_csv(root / "counts.tsv", sep="\t", index=False)
            result = subprocess.run(["Rscript", str(RT / "workflow/scripts/simplseq_to_dcifer.R"), "--cigar", str(root / "counts.tsv"), "--samples", str(root / "samples.csv"), "--out", str(root / "dcifer.tsv")], cwd=root, text=True, capture_output=True, timeout=90)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            rows = pd.read_csv(root / "dcifer.tsv", sep="\t", keep_default_na=False)
            self.assertEqual(set(rows.sample_id), {"LabA", "LabB"})
            self.assertTrue(rows.participant_id.eq("").all())
            self.assertTrue(rows.collection_date.eq("").all())


if __name__ == "__main__":
    unittest.main()
