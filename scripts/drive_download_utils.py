#!/usr/bin/env python3
"""Shared helpers for the Google Drive downloader scripts."""

import re
from pathlib import Path
from typing import Iterable, List, Mapping, Optional, Tuple, Union

DriveUrlItem = Union[str, Mapping[str, str]]
DownloadJob = Tuple[str, Path]

_DRIVE_ID_PATTERNS = (
    r"/folders/([a-zA-Z0-9_-]+)",
    r"/file/d/([a-zA-Z0-9_-]+)",
    r"[?&]id=([a-zA-Z0-9_-]+)",
)


def import_gdown():
    """Import gdown lazily so config-only code can still be tested."""
    try:
        import gdown
    except ImportError as exc:
        raise RuntimeError(
            "Missing dependency 'gdown'. Install it with: python3 -m pip install gdown"
        ) from exc
    return gdown


def extract_drive_id(url: str) -> Optional[str]:
    for pattern in _DRIVE_ID_PATTERNS:
        match = re.search(pattern, url)
        if match:
            return match.group(1)
    return None


def is_folder_url(url: str) -> bool:
    return "/drive/folders/" in url or "/folders/" in url


def is_file_url(url: str) -> bool:
    return "/file/d/" in url or "id=" in url


def target_has_files(path: Path) -> bool:
    return path.exists() and any(path.iterdir())


def download_drive_folder(url: str, output_dir: Path, gdown_module=None) -> None:
    folder_id = extract_drive_id(url)
    if not folder_id:
        raise ValueError("Could not extract folder ID from Google Drive URL.")

    gdown_module = gdown_module or import_gdown()
    output_dir.mkdir(parents=True, exist_ok=True)
    print("Downloading folder to: {0}".format(output_dir))
    gdown_module.download_folder(
        id=folder_id,
        output=str(output_dir),
        quiet=False,
        use_cookies=False,
        remaining_ok=True,
    )


def download_drive_file(url: str, output_dir: Path, gdown_module=None) -> None:
    file_id = extract_drive_id(url)
    if not file_id:
        raise ValueError("Could not extract file ID from Google Drive URL.")

    gdown_module = gdown_module or import_gdown()
    output_dir.mkdir(parents=True, exist_ok=True)
    print("Downloading file to: {0}".format(output_dir))
    gdown_module.download(
        id=file_id,
        output=str(output_dir),
        quiet=False,
        fuzzy=True,
        use_cookies=False,
    )


def download_url(url: str, output_dir: Path, gdown_module=None) -> None:
    if is_folder_url(url):
        download_drive_folder(url, output_dir, gdown_module=gdown_module)
        return

    if is_file_url(url):
        download_drive_file(url, output_dir, gdown_module=gdown_module)
        return

    raise ValueError("Unsupported Google Drive URL. Use shared file/folder links.")


def build_download_jobs(
    drive_urls: Iterable[DriveUrlItem],
    scripts_dir: Path,
    default_target_dir: str,
) -> List[DownloadJob]:
    jobs = []
    default_target = (scripts_dir / default_target_dir).resolve()

    for item in drive_urls:
        if isinstance(item, str):
            jobs.append((item, default_target))
            continue

        if "url" not in item:
            raise ValueError("Each dict item in DRIVE_URLS must include a 'url' key.")

        target_dir = item.get("target_dir", default_target_dir)
        jobs.append((item["url"], (scripts_dir / target_dir).resolve()))

    return jobs
