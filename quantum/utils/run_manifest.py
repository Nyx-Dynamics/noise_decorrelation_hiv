"""
Run manifest utilities: unique run IDs, environment/git capture, and JSON writing.
This supports per-run provenance and non-overwriting output directories.
"""
from __future__ import annotations

import json
import os
import platform
import subprocess
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class Manifest:
    run_id: str
    timestamp_utc: str
    model_variant: str
    ratio: str
    era: str
    input_files: Dict[str, str]
    cli_args: Dict[str, str]
    mechanism: Dict[str, float]
    data_summary: Dict[str, int]
    environment: Dict[str, str]
    outputs: Dict[str, str]
    validation_notes: Dict[str, str]


def make_run_id() -> str:
    now = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    rnd = os.urandom(4).hex()
    return f"{now}_{rnd}"


def _cmd(args: List[str]) -> Optional[str]:
    try:
        return subprocess.check_output(args, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


def get_git_info() -> Dict[str, str]:
    return {
        'commit': _cmd(['git', 'rev-parse', 'HEAD']) or 'unknown',
        'branch': _cmd(['git', 'rev-parse', '--abbrev-ref', 'HEAD']) or 'unknown',
        'dirty': '1' if ((_cmd(['git', 'status', '--porcelain']) or '') != '') else '0',
    }


def file_sha256(path: Path) -> str:
    h = sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(manifest_path: Path, manifest: Manifest) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, 'w') as f:
        json.dump(asdict(manifest), f, indent=2)


def base_environment() -> Dict[str, str]:
    py = platform.python_version()
    try:
        import numpy, scipy  # type: ignore
        npv, sciv = numpy.__version__, scipy.__version__
    except Exception:
        npv, sciv = 'unknown', 'unknown'
    git = get_git_info()
    return {
        'python': py,
        'numpy': npv,
        'scipy': sciv,
        **{f'git_{k}': v for k, v in git.items()}
    }
