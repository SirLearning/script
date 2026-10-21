from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from infra.utils.graph import opaque_legend_handles, plot_joint_regression


def test_opaque_legend_handles_use_full_color(tmp_path: Path):
    fig, ax = plt.subplots()
    coll = ax.scatter([0.0, 1.0], [0.0, 1.0], c="#d95f02", alpha=0.12, s=12)
    proxies = opaque_legend_handles([coll])
    assert coll.get_alpha() == pytest.approx(0.12)
    assert len(proxies) == 1
    assert proxies[0].get_alpha() == pytest.approx(1.0)
    face = np.asarray(proxies[0].get_markerfacecolor(), dtype=float)
    expect = np.asarray(plt.matplotlib.colors.to_rgb("#d95f02"))
    np.testing.assert_allclose(face[:3], expect, atol=1e-6)
    plt.close(fig)
    out = tmp_path / "joint.legend.png"
    df = pd.DataFrame(
        {
            "x": [1.0, 2.0, 3.0, 4.0],
            "y": [0.2, 0.4, 0.3, 0.8],
            "g": ["a", "a", "b", "b"],
        }
    )
    plot_joint_regression(
        df,
        x_col="x",
        y_col="y",
        group_col="g",
        x_label="X",
        y_label="Y",
        filename=str(out),
        draw_regression=False,
        min_points=1,
        bins_x=4,
        bins_y=4,
        hue_order=["a", "b"],
        palette={"a": "#d95f02", "b": "#4c72b0"},
        hue_alpha={"a": 0.12, "b": 0.28},
    )
    assert out.exists()


def test_plot_joint_regression_also_stack_counts(tmp_path: Path):
    df = pd.DataFrame(
        {
            "x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            "y": [0.2, 0.4, 0.3, 0.8, 0.5, 0.9, 0.6, 0.7],
            "g": ["a", "a", "b", "b", "c", "c", "d", "d"],
        }
    )
    out = tmp_path / "joint.stack.png"
    plot_joint_regression(
        df,
        x_col="x",
        y_col="y",
        group_col="g",
        x_label="X",
        y_label="Y",
        filename=str(out),
        draw_regression=False,
        min_points=1,
        bins_x=8,
        bins_y=8,
        marginal_multiple="fill",
        also_stack_counts=True,
    )
    assert out.exists()
    assert out.stat().st_size > 0
