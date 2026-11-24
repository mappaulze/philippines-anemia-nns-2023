"""Unit tests for the WHO 2024 rules. Run: pytest -q"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import nns_anemia_who2024 as m


def test_elevation_adjustment_matches_who_table4():
    assert m.elevation_adjustment(499) == 0                      # no adjustment below 500 m
    assert abs(m.elevation_adjustment(750) - 4.4) < 0.1          # WHO table: 4 g/L for 500-999 m
    assert abs(m.elevation_adjustment(1750) - 11) < 0.5          # 11 g/L for 1500-1999 m
    assert abs(m.elevation_adjustment(4750) - 33.6) < 0.5        # 33 g/L for 4500-4999 m


def test_cutoffs():
    f = lambda a, fem, psc: float(m.who2024_cutoff([a], [fem], [psc])[0])
    assert f(1.0, False, 99) == 105 and f(3.0, True, 0) == 110 and f(8, False, 99) == 115
    assert f(13, True, 0) == 120 and f(30, True, 0) == 120 and f(30, False, 99) == 130
    assert [f(30, True, p) for p in (2, 5, 8)] == [110, 105, 110]


def test_severity_bounds():
    mild, mod = m.who2024_severity_bounds([1.0, 30, 30], [99, 0, 5])
    assert list(mild) == [95, 110, 95] and list(mod) == [70, 80, 70]


def test_classify_adjusts_for_elevation_and_smoking():
    import pandas as pd
    d = pd.DataFrame({"hemoglobin": [12.2, 12.2, 12.2], "age": [30.0] * 3, "female": [True] * 3, "psc": [0] * 3, "smoker": [0, 0, 1]})
    c = m.classify(d, np.array([0, 1500, 0]))
    assert list(c.anemia24) == [0, 1, 1]       # 122 g/L ok at sea level; <120 after 9.1 g/L (elevation) or 3 g/L (smoking)
