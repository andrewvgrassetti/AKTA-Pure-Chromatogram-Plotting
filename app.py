"""
AKTA Pure Chromatogram Plotting - Streamlit Web Application

Interactive web UI for visualizing AKTA Pure chromatography data.
Upload one or more AKTA CSV exports and configure plotting options to
compare chromatograms across multiple purification runs.
"""

import streamlit as st
import numpy as np

from akta_reader import read_akta_file, extract_plot_data
from plotting import (
    get_viridis_colors,
    create_overlay_plot,
    create_individual_plots,
    normalize_data,
)

st.set_page_config(
    page_title="AKTA Pure Chromatogram Plotter",
    page_icon="📈",
    layout="wide",
)

st.title("📈 AKTA Pure Chromatogram Plotter")
st.markdown(
    "Upload one or more AKTA Pure CSV exports and configure the plotting options below."
)

# ── Sidebar: File Upload & Configuration ──────────────────────────────
with st.sidebar:
    st.header("Upload Files")
    uploaded_files = st.file_uploader(
        "Upload AKTA Pure CSV files",
        type=["csv", "tsv", "txt"],
        accept_multiple_files=True,
        help="Select one or more tab-separated CSV files exported from an AKTA Pure system.",
    )

    st.divider()
    st.header("Plot Settings")

    channel_input = st.text_input(
        "UV channels to plot (comma-separated)",
        value="UV 1_280, UV 2_260",
        help="Enter the column identifiers exactly as they appear in the AKTA export.",
    )
    plot_names = [ch.strip() for ch in channel_input.split(",") if ch.strip()]

    overlay = st.toggle("Overlay all samples", value=True)
    plot_post_inj = st.toggle(
        "Start at sample injection",
        value=True,
        help="Shift the x-axis so that volume 0 corresponds to sample application.",
    )
    plot_fracs = st.toggle("Show fractions", value=False)
    rotate_fracs = st.toggle("Rotate fraction labels", value=False)

    st.subheader("Axis Range")
    use_custom_range = st.toggle("Use custom x-axis range", value=False)
    x_min = st.number_input("X min (mL)", value=0.0, disabled=not use_custom_range)
    x_max = st.number_input("X max (mL)", value=500.0, disabled=not use_custom_range)

    x_label = st.text_input("X-axis label", value="Volume (mL)")
    y_label = st.text_input("Y-axis label", value="mAU")

# ── Main Content ──────────────────────────────────────────────────────
if not uploaded_files:
    st.info("👈 Upload AKTA Pure CSV files using the sidebar to get started.")
    st.stop()

# Allow user to assign sample names
st.subheader("Sample Names")
sample_names = []
cols = st.columns(min(len(uploaded_files), 4))
for idx, f in enumerate(uploaded_files):
    col = cols[idx % len(cols)]
    with col:
        name = st.text_input(
            f"Name for {f.name}",
            value=f.name.replace(".csv", "").replace(".tsv", "").replace(".txt", ""),
            key=f"sample_{idx}",
        )
        sample_names.append(name)

st.divider()

# ── Data Processing ───────────────────────────────────────────────────
colors = get_viridis_colors(len(uploaded_files))
all_sample_data = []
errors = []

progress = st.progress(0, text="Reading files…")

for idx, f in enumerate(uploaded_files):
    progress.progress(
        (idx + 1) / len(uploaded_files),
        text=f"Processing {f.name}…",
    )

    try:
        df = read_akta_file(f)
    except Exception as e:
        errors.append(f"**{f.name}**: {e}")
        continue

    result = extract_plot_data(
        df,
        plot_names=plot_names,
        plot_post_inj=plot_post_inj,
        plot_fracs=plot_fracs,
    )

    plot_data = result["plot_data"]
    frac_data = result["frac_data"]

    if not plot_data:
        errors.append(
            f"**{f.name}**: No matching channels found. "
            f"Looked for: {', '.join(plot_names)}"
        )
        continue

    # Determine x-range
    if use_custom_range:
        min_x, max_x = x_min, x_max
    else:
        min_x = min(
            np.nanmin(d["mL"].values) for d in plot_data.values()
        )
        max_x = max(
            np.nanmax(d["mL"].values) for d in plot_data.values()
        )

    # Normalize (shift baseline to zero)
    plot_data, max_y = normalize_data(plot_data, min_x, max_x)

    all_sample_data.append({
        "sample": sample_names[idx],
        "color": colors[idx],
        "plot_data": plot_data,
        "frac_data": frac_data,
        "min_x": min_x,
        "max_x": max_x,
        "max_y": max_y,
    })

progress.empty()

# Show any errors
for err in errors:
    st.warning(err)

if not all_sample_data:
    st.error("No valid data could be loaded. Please check the uploaded files.")
    st.stop()

# ── Plot ──────────────────────────────────────────────────────────────
ax_names = [x_label, y_label]

if overlay:
    fig = create_overlay_plot(
        all_sample_data,
        ax_names=ax_names,
        plot_fracs=plot_fracs,
        rotate_fracs=rotate_fracs,
    )
else:
    fig = create_individual_plots(
        all_sample_data,
        ax_names=ax_names,
        plot_fracs=plot_fracs,
        rotate_fracs=rotate_fracs,
    )

st.plotly_chart(fig, use_container_width=True)

# ── Data Summary ──────────────────────────────────────────────────────
with st.expander("📊 Data Summary"):
    for sd in all_sample_data:
        st.markdown(f"**{sd['sample']}**")
        for ch_name, ch_df in sd["plot_data"].items():
            valid = ch_df.dropna()
            st.write(
                f"  - {ch_name}: {len(valid)} data points, "
                f"mL range [{valid['mL'].min():.1f}, {valid['mL'].max():.1f}], "
                f"max value {valid['data'].max():.2f}"
            )
