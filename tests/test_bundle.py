from hashlib import sha256
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile

from syutyu import export_vmd
from syutyu.report.bundle import iter_dump_blocks
from support import tiny_dump


class BundleTests(unittest.TestCase):
    def test_complete_bundle_allowlist_and_original_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dump = tiny_dump(tmp)
            (tmp / 'credentials.toml').write_text('DO_NOT_EXPORT')
            original = dump.read_bytes()
            output = tmp / 'package with spaces.zip'
            record = export_vmd(dump, output)
            self.assertEqual(record['source_frames'], 3)
            with ZipFile(output) as archive:
                self.assertIsNone(archive.testzip())
                self.assertEqual(set(archive.namelist()), {'trajectory.lammpstrj', 'manifest.json', 'load-vmd.tcl', 'README.txt'})
                self.assertEqual(archive.read('trajectory.lammpstrj'), original)
                self.assertEqual(sha256(original).hexdigest(), record['trajectory_sha256'])
                self.assertIn('file dirname [file normalize [info script]]', archive.read('load-vmd.tcl').decode())
            with self.assertRaises(FileExistsError):
                export_vmd(dump, output)
            self.assertEqual(dump.read_bytes(), original)

    def test_stride_counts_whole_frames(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dump = tiny_dump(tmp)
            expected = list(iter_dump_blocks(dump))[::2]
            record = export_vmd(dump, tmp / 'stride.zip', stride=2)
            self.assertEqual(record['exported_frames'], 2)
            with ZipFile(tmp / 'stride.zip') as archive:
                self.assertEqual(archive.read('trajectory.lammpstrj'), b''.join(expected))
            with self.assertRaises(ValueError):
                export_vmd(dump, tmp / 'bad.zip', stride=0)

    def test_malformed_dump_leaves_no_published_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dump = tmp / 'incomplete.dump'
            dump.write_text('ITEM: TIMESTEP\n0\n')
            with self.assertRaises(ValueError):
                export_vmd(dump, tmp / 'bad.zip')
            self.assertFalse((tmp / 'bad.zip').exists())


if __name__ == '__main__':
    unittest.main()
