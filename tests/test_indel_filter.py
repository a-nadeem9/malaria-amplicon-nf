"""Exercise the production R postprocessor and its downstream Python consumer.

Run with the workflow environment on PATH: python -m unittest discover -s tests
"""

import csv
import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / "runtime"
sys.path.insert(0, str(RUNTIME / "src"))
from simplseq.runner import run_nextflow


class IndelParameterTests(unittest.TestCase):
    def test_app_records_numeric_ratio_cutoff(self):
        for profile in ("local", "reproducible"):
            with self.subTest(profile=profile), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                run_nextflow(root / "samples.csv", root / "out", root=RUNTIME,
                             profile=profile, dry_run=True)
                data = json.loads((root / "out/parameters.json").read_text())
                self.assertEqual(data["analysis_parameters"]["indel_filter"], 0.90)
                self.assertIn("indel_filter", data["pinned_parameter_keys"])
                self.assertEqual(data["analysis_parameters"]["trim_right"], "0,0")
                self.assertEqual(data["analysis_parameters"]["dada2_pool"], "1")


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is required; use the workflow environment")
class IndelPostprocessorTests(unittest.TestCase):
    def run_r(self, *args):
        result = subprocess.run(["Rscript", *map(str, args)], capture_output=True,
                                text=True, timeout=120,
                                env={**os.environ, "SIMPLSEQ_THREADS": "1"})
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_filter_decisions_and_invalid_inputs(self):
        self.run_r(ROOT / "tests/test_asv_filters.R", RUNTIME)

    def test_postprocessor_defaults_override_and_cigar_consumer(self):
        rng = random.Random(7841)
        reference = "".join(rng.choices("ACGT", k=200))
        sequences = [reference[:179], reference[:180], reference,
                     reference + "A" * 22, reference + "A" * 23]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "reference.fa").write_text(">Toy\n" + reference + "\n")
            # A secondary reference must not overwrite the primary target length.
            (root / "secondary.fa").write_text(">Toy\n" + reference + "C" * 40 + "\n")
            (root / "filters.tsv").write_text("id\tSNP\nToy\t100\n")
            (root / "samples.csv").write_text(
                "sample_id,biological_sample_id,sample_type,participant_id,collection_date\n"
                "sample1,visit1,sample,person1,2026-01-01\n")
            with (root / "seqtab.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample", *sequences])
                writer.writerow(["sample1", *([100] * len(sequences))])

            def process(name, *extra):
                output = root / name
                self.run_r(RUNTIME / "workflow/scripts/postProc_dada2.R",
                           "-s", root / "seqtab.tsv", "--samples", root / "samples.csv",
                           "-ref", root / "reference.fa", "-ref2", root / "secondary.fa",
                           "--strain2", "secondary", "--snv_filter", root / "filters.tsv",
                           "-o", output, "--fasta", *extra)
                with output.open(newline="") as handle:
                    return output, {row["hapid"]: row for row in csv.DictReader(handle, delimiter="\t")}

            output, rows = process("default.tsv")
            self.assertEqual([rows[f"ASV{i}"]["indel_filter"] for i in range(1, 6)],
                             ["FAIL", "PASS", "PASS", "PASS", "FAIL"])
            self.assertTrue(all(float(row["reference_length"]) == 200 for row in rows.values()))
            self.assertTrue(all(row["indel_filter_method"] == "broad_length_ratio" for row in rows.values()))
            self.assertEqual(float(rows["ASV2"]["indel_length_ratio"]), 0.90)
            self.assertEqual(sum(line.startswith(">") for line in (root / "ASVSeqs.fasta").read_text().splitlines()), 5)

            spec = importlib.util.spec_from_file_location("cigar_indel_test", RUNTIME / "workflow/scripts/ASV_to_CIGAR.py")
            cigar = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(cigar)
            self.assertEqual(set(cigar.parse_asv_table(output)["Toy"]), {"ASV2", "ASV3", "ASV4"})

            # The explicit cutoff reaches the actual R CLI and changes decisions.
            _, strict = process("strict.tsv", "--indel_filter", "1")
            self.assertEqual([key for key, row in strict.items() if row["indel_filter"] == "PASS"], ["ASV3"])


if __name__ == "__main__":
    unittest.main()
