"""Check pooling configuration and command forwarding without sequencing data.

Run with the workflow environment: python -m unittest discover -s tests
"""

import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

RUNTIME = Path(__file__).resolve().parents[1] / "runtime"
sys.path.insert(0, str(RUNTIME / "src"))

from simplseq.runner import analysis_parameters, run_nextflow


def load_script(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


helper = load_script("nf_helpers", RUNTIME / "workflow/bin/simplseq_nf_helpers.py")
pipeline = load_script("amplicon_pipeline", RUNTIME / "workflow/scripts/AmpliconPipeline.py")


class PoolingTests(unittest.TestCase):
    def write_config(self, root, mode=None):
        args = ["helper", "write-pipeline-json", "--meta", str(root / "meta.tsv"),
                "--out", str(root / "inputs.json")]
        for key, value in (("max-ee", "5,5"), ("trim-right", "0,0"),
                           ("min-len", "50"), ("trunc-q", "5,5"),
                           ("max-consist", "10"), ("omega-a", "1e-120"),
                           ("just-concatenate", "0"), ("primers-fwd", str(root / "primers.fa")),
                           ("primers-rev", str(root / "primers.fa")),
                           ("overlap-primers-fwd", str(root / "primers.fa")),
                           ("overlap-primers-rev", str(root / "primers.fa"))):
            args.extend(["--" + key, value])
        if mode is not None:
            args.extend(["--dada2-pool", mode])
        with patch.object(sys, "argv", args):
            helper.main()
        return json.loads((root / "inputs.json").read_text())

    def test_app_records_pooling_and_no_end_clipping(self):
        for profile in ("local", "reproducible"):
            with self.subTest(profile=profile), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                run_nextflow(root / "samples.csv", root / "out", root=RUNTIME,
                             profile=profile, dry_run=True)
                recorded = json.loads((root / "out/parameters.json").read_text())
                self.assertEqual(recorded["dada2_pool"], "1")
                self.assertIn("dada2_pool", recorded["pinned_parameter_keys"])
                self.assertEqual(recorded["analysis_parameters"], analysis_parameters(profile))
                self.assertEqual(recorded["analysis_parameters"]["trim_right"], "0,0")

    def test_forwarding_to_every_dada2_route(self):
        for mode in (None, "0", "1"):
            for route in ("general", "iseq", "demux_by_amp"):
                with self.subTest(mode=mode, route=route), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    (root / "meta.tsv").write_text("CtrlTest\tunused_R1\tunused_R2\n")
                    (root / "primers.fa").write_text(">TestTarget\nACGT\n")
                    config = self.write_config(root, mode)
                    self.assertEqual(config["dada2_pool"], mode or "1")
                    commands = []

                    def capture(command, label):
                        commands.append(command)
                        if Path(command[1]).name == "runDADA2.R":
                            dest = Path(command[command.index("-d") + 1])
                            name = command[command.index("-o") + 1]
                            (dest / name).write_text("sample\tACGT\nCtrlTest\t10\n")

                    def create_meta(_source, destination, **_patterns):
                        Path(destination).write_text("CtrlTest\tunused_R1\tunused_R2\n")

                    def capture_tee(*streams):
                        for log in streams[1:]:
                            log.close()
                        return io.StringIO()

                    args = ["pipeline", "--json", str(root / "inputs.json"),
                            "--skip_preprocess", "--workers", "1"]
                    if route == "demux_by_amp":
                        args.append("--demux_by_amp")
                    else:
                        args.append("--keep_primers")
                        if route == "iseq":
                            args.append("--iseq")
                    with patch.object(sys, "argv", args), \
                         patch.object(sys, "stdout", io.StringIO()), \
                         patch.object(sys, "stderr", io.StringIO()), \
                         patch.object(pipeline, "Tee", side_effect=capture_tee), \
                         patch.object(pipeline, "run_subprocess", side_effect=capture), \
                         patch.object(pipeline, "run_parallel_jobs"), \
                         patch.object(pipeline, "create_meta", side_effect=create_meta):
                        pipeline.main()
                    self.assertEqual(len(commands), 2 if route == "iseq" else 1)
                    for command in commands:
                        self.assertEqual(command[command.index("--pool") + 1], mode or "1")
                        self.assertEqual(command[command.index("-tR") + 1], "0,0")


if __name__ == "__main__":
    unittest.main()
