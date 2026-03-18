"""
akta_reader.py - Module for reading AKTA Pure chromatography CSV files.

Handles the specific format of AKTA Pure export files:
- UTF-16LE encoded tab-separated values
- Column pairs: identifier column followed by data column
- Row 1: column headers (identifiers), Row 2: units, Row 3+: data
"""

import pandas as pd
import numpy as np


def get_column_data(df, identifier, subtract_start=False, subset_start=0.0, data_is_numeric=True):
    """Extract a pair of columns (mL and data) from an AKTA dataframe by identifier.

    In the AKTA CSV format, each measurement is stored as a pair of columns:
    the first column contains the mL (volume) values with the identifier in row 0,
    and the next column contains the corresponding data values.

    Args:
        df: pandas DataFrame loaded from an AKTA CSV file.
        identifier: The identifier string to search for in the first row.
        subtract_start: If True, subtract subset_start from mL values.
        subset_start: Value to subtract from mL when subtract_start is True.
        data_is_numeric: If True, convert data column to numeric; otherwise keep as string.

    Returns:
        A DataFrame with 'mL' and 'data' columns, or None if identifier not found.
    """
    # Find the column index where the first row matches the identifier
    index = None
    for col_idx, col in enumerate(df.columns):
        first_val = df.iloc[0, col_idx] if len(df) > 0 else None
        if first_val == identifier:
            index = col_idx
            break

    if index is None:
        return None

    # Row 0 is the identifier, row 1 is typically units, row 2+ are data
    ml_series = df.iloc[2:, index]
    data_series = df.iloc[2:, index + 1]

    # Convert mL to numeric
    ml_values = pd.to_numeric(ml_series, errors="coerce").values
    if subtract_start:
        ml_values = ml_values - subset_start

    # Convert data conditionally
    if data_is_numeric:
        data_values = pd.to_numeric(data_series, errors="coerce").values
    else:
        data_values = data_series.values.astype(str)

    result = pd.DataFrame({"mL": ml_values, "data": data_values})
    return result


def read_akta_file(file_obj):
    """Read an AKTA Pure CSV/TSV file.

    AKTA Pure exports files as tab-separated values with UTF-16LE encoding.

    Args:
        file_obj: A file path (str) or file-like object (e.g., from Streamlit upload).

    Returns:
        A pandas DataFrame with all columns as strings.
    """
    return pd.read_csv(
        file_obj,
        sep="\t",
        encoding="utf-16-le",
        header=None,
        dtype=str,
    )


def find_injection_start(df):
    """Find the volume at which sample application begins.

    Args:
        df: pandas DataFrame from read_akta_file.

    Returns:
        The mL value at which 'Sample Application' occurs, or 0.0 if not found.
    """
    run_log = get_column_data(df, "Run Log", data_is_numeric=False)
    if run_log is None:
        return 0.0

    mask = run_log["data"] == "Sample Application"
    indices = run_log.index[mask]
    if len(indices) == 0:
        return 0.0

    # Get the next row's mL value (equivalent to R's +1 index offset)
    next_idx = indices[0] + 1
    if next_idx < len(run_log):
        val = run_log.iloc[next_idx]["mL"]
        if not np.isnan(val):
            return val
    return 0.0


def extract_plot_data(df, plot_names, plot_post_inj=True, plot_fracs=False):
    """Extract all relevant data from an AKTA dataframe for plotting.

    Args:
        df: pandas DataFrame from read_akta_file.
        plot_names: List of channel identifiers (e.g., ['UV 1_280', 'UV 2_260']).
        plot_post_inj: If True, shift mL values so injection start is at 0.
        plot_fracs: If True, also extract fraction data.

    Returns:
        A dict with keys:
            'plot_data': dict mapping channel name -> DataFrame with 'mL' and 'data'
            'frac_data': DataFrame with 'mL' and 'data' for fractions, or None
            'subset_start': The injection start volume used for shifting
    """
    subset_start = 0.0
    if plot_post_inj:
        subset_start = find_injection_start(df)

    plot_data = {}
    for name in plot_names:
        col_data = get_column_data(
            df, name,
            subtract_start=plot_post_inj,
            subset_start=subset_start,
        )
        if col_data is not None:
            plot_data[name] = col_data

    frac_data = None
    if plot_fracs:
        frac_data = get_column_data(
            df, "Fraction",
            subtract_start=plot_post_inj,
            subset_start=subset_start,
            data_is_numeric=False,
        )

    return {
        "plot_data": plot_data,
        "frac_data": frac_data,
        "subset_start": subset_start,
    }
