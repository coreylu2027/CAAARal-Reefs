import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import drive_download_utils as utils


class FakeGdown(object):
    def __init__(self):
        self.folder_calls = []
        self.file_calls = []

    def download_folder(self, **kwargs):
        self.folder_calls.append(kwargs)

    def download(self, **kwargs):
        self.file_calls.append(kwargs)


class ExtractDriveIdTests(unittest.TestCase):
    def test_extracts_folder_id(self):
        url = "https://drive.google.com/drive/folders/abc123_DEF?usp=sharing"
        self.assertEqual(utils.extract_drive_id(url), "abc123_DEF")

    def test_extracts_file_id(self):
        url = "https://drive.google.com/file/d/abc123_DEF/view?usp=sharing"
        self.assertEqual(utils.extract_drive_id(url), "abc123_DEF")

    def test_returns_none_for_non_drive_url(self):
        self.assertIsNone(utils.extract_drive_id("https://example.com/data.csv"))


class BuildDownloadJobsTests(unittest.TestCase):
    def test_uses_default_target_for_plain_urls(self):
        scripts_dir = Path("/tmp/project/scripts")
        jobs = utils.build_download_jobs(
            ["https://drive.google.com/drive/folders/abc123"],
            scripts_dir,
            "../data",
        )

        self.assertEqual(
            jobs,
            [
                (
                    "https://drive.google.com/drive/folders/abc123",
                    Path("/tmp/project/data").resolve(),
                )
            ],
        )

    def test_allows_per_job_target_override(self):
        scripts_dir = Path("/tmp/project/scripts")
        jobs = utils.build_download_jobs(
            [
                {
                    "url": "https://drive.google.com/file/d/abc123/view?usp=sharing",
                    "target_dir": "../custom",
                }
            ],
            scripts_dir,
            "../data",
        )

        self.assertEqual(
            jobs,
            [
                (
                    "https://drive.google.com/file/d/abc123/view?usp=sharing",
                    Path("/tmp/project/custom").resolve(),
                )
            ],
        )

    def test_requires_url_key(self):
        with self.assertRaises(ValueError):
            utils.build_download_jobs(
                [{"target_dir": "../data"}],
                Path("/tmp/project/scripts"),
                "../data",
            )


class DownloadUrlTests(unittest.TestCase):
    def test_dispatches_folder_downloads(self):
        fake_gdown = FakeGdown()
        with tempfile.TemporaryDirectory() as tmp_dir:
            with contextlib.redirect_stdout(io.StringIO()):
                utils.download_url(
                    "https://drive.google.com/drive/folders/abc123?usp=sharing",
                    Path(tmp_dir),
                    gdown_module=fake_gdown,
                )

        self.assertEqual(len(fake_gdown.folder_calls), 1)
        self.assertEqual(fake_gdown.folder_calls[0]["id"], "abc123")
        self.assertEqual(fake_gdown.file_calls, [])

    def test_dispatches_file_downloads(self):
        fake_gdown = FakeGdown()
        with tempfile.TemporaryDirectory() as tmp_dir:
            with contextlib.redirect_stdout(io.StringIO()):
                utils.download_url(
                    "https://drive.google.com/file/d/abc123/view?usp=sharing",
                    Path(tmp_dir),
                    gdown_module=fake_gdown,
                )

        self.assertEqual(len(fake_gdown.file_calls), 1)
        self.assertEqual(fake_gdown.file_calls[0]["id"], "abc123")
        self.assertEqual(fake_gdown.folder_calls, [])

    def test_rejects_unsupported_urls(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            with self.assertRaises(ValueError):
                utils.download_url(
                    "https://example.com/file.csv",
                    Path(tmp_dir),
                    gdown_module=FakeGdown(),
                )


if __name__ == "__main__":
    unittest.main()
