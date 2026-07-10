import numpy as np
import pandas as pd
import pytest

from genetics.genomics.variant.popdep import _percentile_upper_scatter_limits
from infra.utils.graph import adaptive_axis_max, axis_limits_with_adaptive_upper


def test_adaptive_axis_max_peels_isolated_high_outlier():
    values = np.array([0.0, 1.0, 2.0, 3.0, 100.0])
    assert adaptive_axis_max(values, min_points=5) == pytest.approx(3.0)


def test_adaptive_axis_max_keeps_dense_upper_cluster():
    values = np.array([0.0, 1.0, 2.0, 98.0, 99.0, 100.0])
    assert adaptive_axis_max(values, min_points=5) is None


def test_adaptive_axis_max_keeps_two_point_upper_cluster():
    values = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 50.0, 51.0])
    assert adaptive_axis_max(values, min_points=5) is None


def test_adaptive_axis_max_peels_chained_isolated_outliers():
    values = np.array([0.0, 1.0, 100.0, 200.0])
    assert adaptive_axis_max(values, min_points=4) == pytest.approx(1.0)


def test_adaptive_axis_max_uniform_data_no_trim():
    values = np.linspace(0.0, 10.0, 50)
    assert adaptive_axis_max(values) is None


def test_adaptive_axis_max_too_few_points():
    assert adaptive_axis_max([1.0, 2.0, 100.0], min_points=10) is None


def test_axis_limits_with_adaptive_upper_respects_lower_override():
    values = np.array([0.5, 1.0, 2.0, 3.0, 100.0])
    lim = axis_limits_with_adaptive_upper(values, min_points=5, lower_override=0.0)
    assert lim is not None
    lo, hi = lim
    assert lo == pytest.approx(0.0)
    assert hi > 3.0
    assert hi < 100.0


_POPDEP_DEPTH_SPEC = {
    "tag": "depth",
    "mean": "Depth_Mean",
    "sd": "Depth_SD",
    "var": "Depth_Var",
    "cv": "Depth_CV",
}


def test_percentile_upper_scatter_keeps_low_values():
    n = 2000
    rng = np.random.default_rng(0)
    df = pd.DataFrame(
        {
            "Depth_Mean": rng.uniform(0.0, 10.0, n),
            "Depth_SD": rng.uniform(0.0, 5.0, n),
        }
    )
    trimmed, x_lim, y_lim = _percentile_upper_scatter_limits(
        df,
        "Depth_Mean",
        "Depth_SD",
        _POPDEP_DEPTH_SPEC,
        upper_pct=99.9,
        min_points=10,
    )
    assert trimmed is not None
    assert float(trimmed["Depth_SD"].min()) < float(df["Depth_SD"].quantile(0.01))
    assert x_lim is not None and x_lim[0] == pytest.approx(0.0)
    assert y_lim is not None and y_lim[0] == pytest.approx(0.0)
