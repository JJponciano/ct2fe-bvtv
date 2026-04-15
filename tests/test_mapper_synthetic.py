#!/usr/bin/env python3
"""Small regression test for BV/TV mapping."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from map_bvtv_to_fe import Mesh, Volume, map_bvtv  # noqa: E402


def main() -> int:
    data = np.zeros((7, 7, 7), dtype=np.float32)
    data[2:5, 2:5, 2:5] = 100.0
    volume = Volume(data=data, affine=np.eye(4), spacing=np.ones(3))
    mesh = Mesh(
        nodes={
            1: np.array([2.5, 2.5, 2.5]),
            2: np.array([3.5, 2.5, 2.5]),
            3: np.array([3.5, 3.5, 2.5]),
            4: np.array([2.5, 3.5, 2.5]),
            5: np.array([2.5, 2.5, 3.5]),
            6: np.array([3.5, 2.5, 3.5]),
            7: np.array([3.5, 3.5, 3.5]),
            8: np.array([2.5, 3.5, 3.5]),
            9: np.array([0.0, 0.0, 0.0]),
            10: np.array([1.0, 0.0, 0.0]),
            11: np.array([1.0, 1.0, 0.0]),
            12: np.array([0.0, 1.0, 0.0]),
            13: np.array([0.0, 0.0, 1.0]),
            14: np.array([1.0, 0.0, 1.0]),
            15: np.array([1.0, 1.0, 1.0]),
            16: np.array([0.0, 1.0, 1.0]),
        },
        elements={
            1: [1, 2, 3, 4, 5, 6, 7, 8],
            2: [9, 10, 11, 12, 13, 14, 15, 16],
        },
    )

    results, threshold = map_bvtv(volume, mesh, threshold=50.0, sphere_diameter_mm=3.0)
    by_id = {row.element_id: row for row in results}
    assert threshold == 50.0
    assert by_id[1].bvtv == 1.0, by_id[1]
    assert by_id[1].youngs_modulus == 8534.64, by_id[1]
    assert by_id[2].bvtv == 0.0, by_id[2]
    assert by_id[2].youngs_modulus == 0.0, by_id[2]
    print("Synthetic BV/TV mapping test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
