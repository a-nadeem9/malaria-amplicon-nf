"""Exercise INDEL filtering through the actual R postprocessor on artificial DNA.

Run in the workflow environment: python -m unittest discover -s tests
"""

import csv
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

RUNTIME = Path(__file__).resolve().parents[1] / "runtime"
RSCRIPT = shutil.which("Rscript") or str(Path(sys.executable).with_name("Rscript"))
sys.path.insert(0, str(RUNTIME / "src"))
from simplseq.runner import run_nextflow


def write_fixture(root):
    rng = random.Random(417)
    dna = lambda n: "".join(rng.choices("ACGT", k=n))
    long_ref, short_ref, inserted = dna(200), dna(100), dna(21)
    cases = [("reference", long_ref, 0, 200, "PASS"),
             ("four_base_deletion", long_ref[:80] + long_ref[84:], 4, 200, "PASS")]
    for distance in (19, 20, 21):
        outcome = "FAIL" if distance > 20 else "PASS"
        cases.append((f"insertion_{distance}", long_ref[:80] + inserted[:distance] + long_ref[80:],
                      distance, 200, outcome))
        cases.append((f"deletion_{distance}", long_ref[:80] + long_ref[80 + distance:],
                      distance, 200, outcome))
    # Eleven inserted and eleven deleted bases leave the total length unchanged.
    cases.append(("opposing_events", long_ref[:40] + inserted[:11] + long_ref[40:140] + long_ref[151:],
                  22, 200, "FAIL"))
    cases.extend((f"short_deletion_{n}", short_ref[:40] + short_ref[40+n:], n, 100,
                  "PASS" if n == 10 else "FAIL") for n in (10, 11))
    (root / "reference.fasta").write_text(f">Test200\n{long_ref}\n>Test100\n{short_ref}\n")
    (root / "filters.tsv").write_text("id\tSNP\nTest200\t1000\nTest100\t1000\n")
    with (root / "seqtab.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + [case[1] for case in cases])
        writer.writerow(["TestSample"] + [1000] * len(cases))
    (root / "samples.csv").write_text(
        "sample_id,biological_sample_id,sample_type,participant_id,collection_date\n"
        "TestSample,TestVisit,sample,TestParticipant,2026-01-01\n")
    return cases


class IndelParameterTests(unittest.TestCase):
    def test_app_records_ten_percent_in_both_profiles(self):
        for profile in ("local", "reproducible"):
            with self.subTest(profile=profile), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                run_nextflow(root / "samples.csv", root / "out", root=RUNTIME,
                             profile=profile, dry_run=True)
                recorded = json.loads((root / "out/parameters.json").read_text())
                self.assertEqual(recorded["analysis_parameters"]["indel_filter"], "0.10")
                self.assertIn("indel_filter", recorded["pinned_parameter_keys"])


@unittest.skipUnless(Path(RSCRIPT).is_file(), "Rscript from the workflow environment is required")
class IndelPostprocessorTests(unittest.TestCase):
    def test_distance_fraction_boundaries_and_opposing_events(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cases = write_fixture(root)
            command = [RSCRIPT, str(RUNTIME / "workflow/scripts/postProc_dada2.R"),
                       "-s", str(root / "seqtab.tsv"), "--samples", str(root / "samples.csv"),
                       "-ref", str(root / "reference.fasta"), "-snv", str(root / "filters.tsv"),
                       "-o", str(root / "mapped.tsv")]
            result = subprocess.run(command, capture_output=True, text=True, timeout=60,
                                    env={**os.environ, "SIMPLSEQ_THREADS": "1"})
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / "mapped.tsv").open() as handle:
                mapped = {row["hapid"]: row for row in csv.DictReader(handle, delimiter="\t")}
            for index, (name, _, distance, reference_length, outcome) in enumerate(cases, 1):
                with self.subTest(case=name):
                    row = mapped[f"ASV{index}"]
                    self.assertEqual(float(row["indel_dist_from_3D7"]), distance)
                    self.assertEqual(float(row["reference_length"]), reference_length)
                    self.assertAlmostEqual(float(row["indel_fraction"]), distance / reference_length)
                    self.assertEqual(float(row["indel_fraction_threshold"]), 0.10)
                    self.assertEqual(row["indel_filter_rule"], "inserted_deleted_bases_over_reference_length")
                    self.assertEqual(row["indel_filter"], outcome)
                    self.assertEqual(row["mapping_status"], "MAPPED" if outcome == "PASS" else "UNMAPPED")

    def test_rejects_absolute_count_as_fraction(self):
        result = subprocess.run(
            [RSCRIPT, str(RUNTIME / "workflow/scripts/postProc_dada2.R"), "--indel_filter", "3"],
            capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("fraction between 0 and 1", result.stderr)


if __name__ == "__main__":
    unittest.main()
