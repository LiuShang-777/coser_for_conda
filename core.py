from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.preprocessing import MinMaxScaler


def _read_trait_table(
    pheno_file: Union[str, Path],
    sep: str,
    index_col: Optional[Union[int, str]],
    drop_first_trait_col: bool,
    fillna: bool,
) -> pd.DataFrame:
    dat_pheno = pd.read_csv(pheno_file, sep=sep, index_col=index_col)

    if drop_first_trait_col:
        dat_pheno = dat_pheno.iloc[:, 1:]

    dat_pheno = dat_pheno.apply(pd.to_numeric, errors="coerce")

    if fillna:
        dat_pheno = dat_pheno.fillna(dat_pheno.mean(numeric_only=True))

    return dat_pheno


def _prepare_traits(
    dat_pheno: pd.DataFrame,
    x_trait: str,
    y_trait: str,
    scale_data: bool,
) -> pd.DataFrame:
    missing = [trait for trait in (x_trait, y_trait) if trait not in dat_pheno.columns]
    if missing:
        raise ValueError(f"Trait column not found: {', '.join(missing)}")

    if scale_data:
        scaler = MinMaxScaler()
        dat_scaled = scaler.fit_transform(dat_pheno)
        return pd.DataFrame(dat_scaled, index=dat_pheno.index, columns=dat_pheno.columns)

    return dat_pheno.copy()


def _regression_stats(dat_trait: pd.DataFrame, x_trait: str, y_trait: str):
    x = dat_trait[x_trait].values.reshape(-1, 1)
    y = dat_trait[y_trait].values

    model = LinearRegression()
    model.fit(x, y)

    slope = float(model.coef_[0])
    intercept = float(model.intercept_)
    y_pred = model.predict(x)
    r2 = float(r2_score(y, y_pred))
    pearson_r = float(pearsonr(y, y_pred)[0])

    return slope, intercept, r2, pearson_r


def calculate_angle_coupling_score(
    pheno_file: Union[str, Path],
    x_trait: str,
    y_trait: str,
    output_file: Optional[Union[str, Path]] = None,
    sep: str = ",",
    index_col: Optional[Union[int, str]] = 0,
    drop_first_trait_col: bool = False,
    scale_data: bool = True,
    fillna: bool = True,
) -> pd.DataFrame:
    """Calculate angle-based CoSEr score: 1 - sin(theta)."""
    dat_pheno = _read_trait_table(
        pheno_file=pheno_file,
        sep=sep,
        index_col=index_col,
        drop_first_trait_col=drop_first_trait_col,
        fillna=fillna,
    )
    dat_trait = _prepare_traits(dat_pheno, x_trait, y_trait, scale_data)
    slope, intercept, r2, pearson_r = _regression_stats(dat_trait, x_trait, y_trait)

    x_shifted = dat_trait[x_trait].values
    y_shifted = dat_trait[y_trait].values - intercept

    distance_to_line = np.abs(slope * x_shifted - y_shifted) / np.sqrt(slope**2 + 1)
    distance_to_zero = np.sqrt(x_shifted**2 + y_shifted**2)

    sin_theta = np.divide(
        distance_to_line,
        distance_to_zero,
        out=np.zeros_like(distance_to_line),
        where=distance_to_zero != 0,
    )
    sin_theta = np.clip(sin_theta, 0, 1)
    coupling_score = 1 - sin_theta

    result = pd.DataFrame(
        {
            "sample": dat_trait.index,
            x_trait: dat_trait[x_trait].values,
            y_trait: dat_trait[y_trait].values,
            f"{y_trait}_shifted": y_shifted,
            f"{x_trait}_{y_trait}_regression_slope": slope,
            f"{x_trait}_{y_trait}_regression_intercept": intercept,
            f"{x_trait}_{y_trait}_r2": r2,
            f"{x_trait}_{y_trait}_pearson_r": pearson_r,
            f"{x_trait}_{y_trait}_angle_cs": coupling_score,
        },
        index=dat_trait.index,
    ).sort_index()

    if output_file is not None:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_file, index=False)

    return result


def calculate_distance_coupling_score(
    pheno_file: Union[str, Path],
    x_trait: str,
    y_trait: str,
    output_file: Optional[Union[str, Path]] = None,
    sep: str = ",",
    index_col: Optional[Union[int, str]] = 0,
    drop_first_trait_col: bool = False,
    scale_data: bool = True,
    fillna: bool = True,
    normalize_distance: bool = True,
) -> pd.DataFrame:
    """Calculate distance-based CoSEr score from the shifted regression line."""
    dat_pheno = _read_trait_table(
        pheno_file=pheno_file,
        sep=sep,
        index_col=index_col,
        drop_first_trait_col=drop_first_trait_col,
        fillna=fillna,
    )
    dat_trait = _prepare_traits(dat_pheno, x_trait, y_trait, scale_data)
    slope, intercept, r2, pearson_r = _regression_stats(dat_trait, x_trait, y_trait)

    x_shifted = dat_trait[x_trait].values
    y_shifted = dat_trait[y_trait].values - intercept

    distance_to_line = np.abs(slope * x_shifted - y_shifted) / np.sqrt(slope**2 + 1)

    result = pd.DataFrame(
        {
            "sample": dat_trait.index,
            x_trait: dat_trait[x_trait].values,
            y_trait: dat_trait[y_trait].values,
            f"{y_trait}_shifted": y_shifted,
            f"{x_trait}_{y_trait}_regression_slope": slope,
            f"{x_trait}_{y_trait}_regression_intercept": intercept,
            f"{x_trait}_{y_trait}_r2": r2,
            f"{x_trait}_{y_trait}_pearson_r": pearson_r,
            f"{x_trait}_{y_trait}_distance_to_regression_line": distance_to_line,
        },
        index=dat_trait.index,
    )

    if normalize_distance:
        max_dist = float(np.max(distance_to_line))
        if max_dist == 0:
            result[f"{x_trait}_{y_trait}_distance_cs"] = 1
        else:
            result[f"{x_trait}_{y_trait}_distance_cs"] = 1 - distance_to_line / max_dist

    result = result.sort_index()

    if output_file is not None:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_file, index=False)

    return result
