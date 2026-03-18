# AKTA Pure Chromatogram Plotter

An interactive web application for visualizing and comparing chromatograms
exported from AKTA Pure protein purification systems.

Built with **Python**, **Streamlit**, and **Plotly**.

![Python](https://img.shields.io/badge/Python-3.9%2B-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-1.30%2B-FF4B4B)

## Features

- **Upload & Compare** – Load multiple AKTA Pure CSV exports and overlay them on a single chart or view them as individual subplots.
- **Channel Selection** – Choose which UV / detection channels to plot (e.g., UV 1_280, UV 2_260).
- **Baseline Normalization** – Automatically shifts curves so the baseline (minimum value) in the visible range is zero.
- **Injection Start Alignment** – Optionally shift the x-axis so that volume 0 corresponds to the sample application event.
- **Fraction Markers** – Optionally display fraction collection boundaries as vertical lines with labels.
- **Custom Axis Ranges** – Restrict the x-axis to a specific volume window.
- **Interactive Plots** – Powered by Plotly: zoom, pan, hover for data values, and export as PNG.

## Quick Start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Run the app

```bash
streamlit run app.py
```

The app will open in your default web browser. Use the sidebar to upload AKTA
Pure CSV files and configure plotting options.

## Input File Format

The application expects **tab-separated** files exported from an AKTA Pure
system. These files are typically UTF-16LE encoded and contain paired columns
(volume + measurement) for each sensor channel.

## Project Structure

```
├── app.py              # Streamlit web application
├── akta_reader.py      # AKTA CSV file reading and parsing
├── plotting.py         # Plotly chart creation
├── requirements.txt    # Python dependencies
└── README.md
```

## Requirements

- Python 3.9+
- See `requirements.txt` for package dependencies
