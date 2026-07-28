#!/usr/bin/env python3

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class BlastWrapperTests(unittest.TestCase):
    def test_config_mode_runs_with_python_yaml_without_rscript(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_root = root / "output"
            ref_root = root / "ref"
            fake_bin = root / "bin"
            fake_pythonpath = root / "pythonpath"

            comparison_dir = output_root / "BlastAnalysis" / "TestComparison" / "Taxon_A"
            comparison_dir.mkdir(parents=True)
            (comparison_dir / "Taxon_A_ASV.fasta").write_text(">asv1\nACGT\n")

            (ref_root / "testdb" / "blast_format_db").mkdir(parents=True)
            (ref_root / "testdb.fasta").write_text(">ref1\nACGT\n")
            (ref_root / "testdb" / "seqid2taxid.map").write_text("ref1\t1\n")
            (ref_root / "testdb" / "blast_format_db" / "testdb_blast.nsq").write_text("")

            fake_bin.mkdir()
            blastn = fake_bin / "blastn"
            blastn.write_text("#!/usr/bin/env bash\nprintf 'asv1\\tref1\\t100\\n'\n")
            blastn.chmod(0o755)

            rscript = fake_bin / "Rscript"
            rscript.write_text("#!/usr/bin/env bash\necho 'Rscript should not be called' >&2\nexit 99\n")
            rscript.chmod(0o755)

            fake_pythonpath.mkdir()
            (fake_pythonpath / "yaml.py").write_text(
                "def safe_load(handle):\n"
                "    result = {}\n"
                "    section = None\n"
                "    pending_list = None\n"
                "    for raw in handle.read().splitlines():\n"
                "        if not raw.strip() or raw.lstrip().startswith('#'):\n"
                "            continue\n"
                "        indent = len(raw) - len(raw.lstrip(' '))\n"
                "        line = raw.split('#', 1)[0].strip()\n"
                "        if not line:\n"
                "            continue\n"
                "        if indent == 0 and ':' in line:\n"
                "            key, value = line.split(':', 1)\n"
                "            key = key.strip()\n"
                "            if value.strip():\n"
                "                result[key] = _parse_scalar(value.strip())\n"
                "                section = None\n"
                "            else:\n"
                "                result[key] = {}\n"
                "                section = key\n"
                "            pending_list = None\n"
                "            continue\n"
                "        if section and indent == 2 and ':' in line:\n"
                "            key, value = line.split(':', 1)\n"
                "            key = key.strip()\n"
                "            if value.strip():\n"
                "                result[section][key] = _parse_scalar(value.strip())\n"
                "                pending_list = None\n"
                "            else:\n"
                "                result[section][key] = []\n"
                "                pending_list = key\n"
                "            continue\n"
                "        if section and pending_list and indent > 2 and line.startswith('-'):\n"
                "            result[section][pending_list].append(_parse_scalar(line[1:].strip()))\n"
                "    return result\n"
                "\n"
                "def _parse_scalar(value):\n"
                "    value = value.strip().strip('\\\"').strip(\"'\")\n"
                "    if value.lower() in {'null', 'none', '~'}:\n"
                "        return None\n"
                "    return value\n"
            )

            config = root / "analysis.yaml"
            config.write_text(
                "project:\n"
                f"  output_dir: {output_root}\n"
                "blast_confirmation:\n"
                "  candidate_comparisons:\n"
                "    - TestComparison\n"
                "  reference_db: testdb\n"
                f"  ref_base_dir: {ref_root}\n"
                "  num_threads: 2\n"
                "  max_target_seqs: 7\n"
            )

            python_dir = Path(shutil.which("python3") or "").parent
            env = {
                "PATH": os.pathsep.join([str(fake_bin), str(python_dir), "/usr/bin", "/bin"]),
                "PYTHONPATH": str(fake_pythonpath),
            }
            completed = subprocess.run(
                ["/bin/bash", "scripts/BLASTWrapper.sh", "--analysis-config", str(config)],
                cwd=REPO_ROOT,
                env=env,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stdout)
            hits = comparison_dir / "Taxon_A_blast_hits.tsv"
            self.assertEqual(hits.read_text(), "asv1\tref1\t100\n")


if __name__ == "__main__":
    unittest.main()
