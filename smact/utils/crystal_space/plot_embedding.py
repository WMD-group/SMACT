"""Utility functions for plotting dimension reduced embeddings using plotly."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

if TYPE_CHECKING:
    from pathlib import Path


def update_layout(
    fig: go.Figure,
    title: str,
    num_row: int = 6,
    num_col: int = 3,
    width: float = 1200,
    height: float = 1800,
):
    """Update layout of a plotly figure."""
    # set axis
    for i in range(1, num_row + 1):
        for j in range(1, num_col + 1):
            fig.update_xaxes(
                showticklabels=False,
                linecolor="black",
                showline=True,
                linewidth=1,
                mirror=True,
                row=i,
                col=j,
            )
            fig.update_yaxes(
                showticklabels=False,
                linecolor="black",
                showline=True,
                linewidth=1,
                mirror=True,
                row=i,
                col=j,
            )

    # set layout
    fig.update_layout(
        title=title,
        title_x=0.5,
        title_font_size=30,
        width=width,
        height=height,
        margin=dict(l=10, r=10, t=80, b=50),
        paper_bgcolor="rgba(255,255,255,1)",
        plot_bgcolor="rgba(255,255,255,1)",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            xanchor="center",
            x=0.5,
            y=-0.04,
            font=dict(size=30),
        ),
    )
    return fig


def plot_reducers_embeddings(
    df_label: pd.DataFrame,
    reducers: list[str],
    embedding_names: list[str],
    embedding_dir: Path,
    save_path: Path,
    symbol: str = "circle",
    title: str = "Embedding Visualization",
):
    """Plot dimension reduction plots."""
    fig = make_subplots(
        rows=6,
        cols=3,
        subplot_titles=[f"{reducer} - {embedding_name}" for embedding_name in embedding_names for reducer in reducers],
        vertical_spacing=0.02,
        horizontal_spacing=0.02,
    )

    # updatee the font size of subplot titles
    for i in fig["layout"]["annotations"]:
        i["font"] = dict(size=25)

    legend_colors = {
        "unlikely": "#D9D9D9",
        "interesting": "#22E000",
        "missing": "#FF1201",
        "standard": "#002FFF",
    }

    for i, embedding_name in enumerate(embedding_names):
        for j, reducer in enumerate(reducers):
            print(f"processing {i} {j}...")
            embedding_data = pd.read_pickle(
                embedding_dir / f"{reducer}_{embedding_name}.pkl",
            )
            embedding_data.columns = ["x", "y"]
            df_plot = embedding_data.join(df_label)
            df_plot = df_plot.sample(frac=1, random_state=42)

            fig.add_trace(
                go.Scatter(
                    x=df_plot["x"],
                    y=df_plot["y"],
                    mode="markers",
                    marker=dict(
                        size=8,
                        color=df_plot["label"].map(legend_colors),
                        opacity=0.8,
                        symbol=symbol,
                        line=dict(width=0.5, color="DarkSlateGrey"),
                    ),
                    showlegend=False,
                    text=df_plot["formula"],
                    hovertemplate=("<b>%{text}</b><br><br>"),
                ),
                row=i + 1,
                col=j + 1,
            )

    # add legend
    for label in legend_colors:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(
                    size=8,
                    color=legend_colors[label],
                    opacity=0.8,
                    symbol=symbol,
                    line=dict(width=0.5, color="DarkSlateGrey"),
                ),
                # make only first letter capital
                name=label.capitalize(),
                showlegend=True,
            ),
            row=1,
            col=1,
        )

    # update layout
    fig = update_layout(fig, title=title)

    if save_path is not None:
        if save_path.suffix == ".html":
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, scale=6)
        print(f"Save to {save_path}")

    return fig
