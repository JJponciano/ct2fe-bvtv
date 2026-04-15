#!/usr/bin/env python3
"""Download the public CT demo volume used by the repository examples."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path


DEFAULT_URL = "https://neurolabusc.github.io/niivue-images/CT_Abdo.nii.gz"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Download the public CT_Abdo NIfTI demo image.")
    parser.add_argument("--url", default=DEFAULT_URL, help="Source URL for the demo CT volume.")
    parser.add_argument("--out", default="data/CT_Abdo.nii.gz", help="Output path.")
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing file.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_path = Path(args.out)
    if out_path.exists() and not args.overwrite:
        print(f"Exists: {out_path}")
        return 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {args.url}")
    urllib.request.urlretrieve(args.url, out_path)
    print(f"Wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
