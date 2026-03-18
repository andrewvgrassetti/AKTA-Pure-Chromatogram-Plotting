"""
plotting.py - Module for creating interactive chromatogram plots with Plotly.

Provides functions to create overlay and individual chromatogram plots
matching the functionality of the original R plotChromMultiple function.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# Viridis-inspired color palette (hex values)
_VIRIDIS_COLORS = [
    "#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
    "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
]


def get_viridis_colors(n):
    """Return n evenly spaced colors from the viridis palette.

    Args:
        n: Number of colors needed.

    Returns:
        List of hex color strings.
    """
    if n <= 0:
        return []
    if n == 1:
        return [_VIRIDIS_COLORS[0]]
    indices = np.linspace(0, len(_VIRIDIS_COLORS) - 1, n).astype(int)
    return [_VIRIDIS_COLORS[i] for i in indices]


def lighten_color(hex_color, factor=0.5):
    """Lighten a hex color by mixing with white.

    Args:
        hex_color: A hex color string like '#440154'.
        factor: Blend factor (0 = no change, 1 = white).

    Returns:
        Lightened hex color string.
    """
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    r = r + (1 - r) * factor
    g = g + (1 - g) * factor
    b = b + (1 - b) * factor
    return "#{:02x}{:02x}{:02x}".format(
        int(r * 255), int(g * 255), int(b * 255)
    )


def normalize_data(plot_data_list, min_x, max_x):
    """Shift each channel's data so the minimum in the x-range is zero.

    Args:
        plot_data_list: dict mapping channel name -> DataFrame with 'mL' and 'data'.
        min_x: Minimum x-axis value.
        max_x: Maximum x-axis value.

    Returns:
        Tuple of (modified plot_data_list, max_y in range).
    """
    max_y = 0.0
    for name, pd_df in plot_data_list.items():
        in_range = (pd_df["mL"] >= min_x) & (pd_df["mL"] <= max_x)
        data_in_range = pd_df.loc[in_range, "data"]
        if len(data_in_range) == 0:
            continue
        min_val = np.nanmin(data_in_range.values)
        pd_df["data"] = pd_df["data"] - min_val
        plot_data_list[name] = pd_df
        max_y = max(max_y, np.nanmax(pd_df.loc[in_range, "data"].values))
    return plot_data_list, max_y


def _add_fraction_lines(fig, frac_data, min_x, max_x, max_y, rotate_fracs=False, row=None, col=None):
    """Add fraction boundary lines and labels to a figure.

    Args:
        fig: Plotly figure.
        frac_data: DataFrame with 'mL' and 'data' columns for fractions.
        min_x: Minimum x for filtering.
        max_x: Maximum x for filtering.
        max_y: Maximum y for label placement.
        rotate_fracs: Whether to rotate fraction labels.
        row: Subplot row (None for single plot).
        col: Subplot column (None for single plot).
    """
    if frac_data is None:
        return

    in_range = (
        frac_data["mL"].notna()
        & (frac_data["mL"] >= min_x)
        & (frac_data["mL"] <= max_x)
    )
    positions = frac_data.loc[in_range, "mL"].values
    labels = frac_data.loc[in_range, "data"].values

    for pos, label in zip(positions, labels):
        # Use numeric mL values for fraction positions
        try:
            pos_val = float(pos)
        except (ValueError, TypeError):
            continue

        kwargs = {}
        if row is not None:
            kwargs["row"] = row
            kwargs["col"] = col

        fig.add_vline(
            x=pos_val, line_dash="dot", line_color="gray",
            line_width=1, **kwargs,
        )
        fig.add_annotation(
            x=pos_val,
            y=max_y * 0.95,
            text=str(label),
            showarrow=False,
            textangle=-45 if rotate_fracs else 0,
            font=dict(size=9, color="gray"),
            **kwargs,
        )


def create_overlay_plot(all_sample_data, ax_names=None, plot_fracs=False, rotate_fracs=False):
    """Create an overlay plot with all samples on a single chart.

    Args:
        all_sample_data: List of dicts, each with keys:
            'sample': str, sample name
            'color': str, hex color
            'plot_data': dict of channel -> DataFrame
            'frac_data': DataFrame or None
            'min_x': float
            'max_x': float
            'max_y': float
        ax_names: List of axis labels [x_label, y_label].
        plot_fracs: Whether to show fraction markers.
        rotate_fracs: Whether to rotate fraction labels.

    Returns:
        A Plotly Figure.
    """
    if ax_names is None:
        ax_names = ["Volume (mL)", "mAU"]

    # Compute global ranges
    global_min_x = min(d["min_x"] for d in all_sample_data)
    global_max_x = max(d["max_x"] for d in all_sample_data)
    global_max_y = max(d["max_y"] for d in all_sample_data)

    fig = go.Figure()

    for sample_data in all_sample_data:
        sample_name = sample_data["sample"]
        col_main = sample_data["color"]
        plot_data = sample_data["plot_data"]

        for j, (channel_name, pd_df) in enumerate(plot_data.items()):
            color = col_main if j == 0 else lighten_color(col_main)
            dash = "solid" if j == 0 else "dot"
            # Show legend entry for the first channel of each sample
            legend_name = sample_name if j == 0 else f"{sample_name} ({channel_name})"

            fig.add_trace(go.Scatter(
                x=pd_df["mL"],
                y=pd_df["data"],
                mode="lines",
                name=legend_name,
                line=dict(color=color, width=2.5, dash=dash),
                hovertemplate=f"{sample_name} - {channel_name}<br>"
                              f"Volume: %{{x:.1f}} mL<br>"
                              f"Value: %{{y:.2f}}<extra></extra>",
            ))

        # Add fraction lines for this sample
        if plot_fracs:
            _add_fraction_lines(
                fig, sample_data["frac_data"],
                sample_data["min_x"], sample_data["max_x"],
                global_max_y, rotate_fracs,
            )

    fig.update_layout(
        title="Overlay of All Samples",
        xaxis_title=ax_names[0],
        yaxis_title=ax_names[1],
        xaxis=dict(range=[global_min_x, global_max_x]),
        yaxis=dict(range=[0, global_max_y * 1.05]),
        legend=dict(x=1, y=1, xanchor="right"),
        template="plotly_white",
        height=600,
    )

    return fig


def create_individual_plots(all_sample_data, ax_names=None, plot_fracs=False, rotate_fracs=False):
    """Create individual subplots, one per sample.

    Args:
        all_sample_data: List of dicts (same format as create_overlay_plot).
        ax_names: List of axis labels [x_label, y_label].
        plot_fracs: Whether to show fraction markers.
        rotate_fracs: Whether to rotate fraction labels.

    Returns:
        A Plotly Figure with subplots.
    """
    if ax_names is None:
        ax_names = ["Volume (mL)", "mAU"]

    n = len(all_sample_data)
    fig = make_subplots(
        rows=n, cols=1,
        subplot_titles=[d["sample"] for d in all_sample_data],
        vertical_spacing=0.08,
    )

    for i, sample_data in enumerate(all_sample_data):
        row = i + 1
        col_main = sample_data["color"]
        plot_data = sample_data["plot_data"]

        for j, (channel_name, pd_df) in enumerate(plot_data.items()):
            color = col_main if j == 0 else lighten_color(col_main)
            dash = "solid" if j == 0 else "dot"

            fig.add_trace(
                go.Scatter(
                    x=pd_df["mL"],
                    y=pd_df["data"],
                    mode="lines",
                    name=channel_name,
                    line=dict(color=color, width=2.5, dash=dash),
                    showlegend=(i == 0),
                    hovertemplate=f"{sample_data['sample']} - {channel_name}<br>"
                                  f"Volume: %{{x:.1f}} mL<br>"
                                  f"Value: %{{y:.2f}}<extra></extra>",
                ),
                row=row, col=1,
            )

        fig.update_xaxes(
            title_text=ax_names[0] if row == n else "",
            range=[sample_data["min_x"], sample_data["max_x"]],
            row=row, col=1,
        )
        fig.update_yaxes(
            title_text=ax_names[1],
            range=[0, sample_data["max_y"] * 1.05],
            row=row, col=1,
        )

        if plot_fracs:
            _add_fraction_lines(
                fig, sample_data["frac_data"],
                sample_data["min_x"], sample_data["max_x"],
                sample_data["max_y"], rotate_fracs,
                row=row, col=1,
            )

    fig.update_layout(
        height=400 * n,
        template="plotly_white",
        showlegend=True,
    )

    return fig
