"""Utilities for peptide/protein feature planning and filtering.

This module helps you:
1) load sequence data from Excel,
2) define requested descriptors,
3) exclude descriptors you already computed,
4) keep PAAC with dimension 20 + 2*lambda when requested.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


# Features you said you want to add (from BBPpred-style descriptor pool).
REQUESTED_FEATURES = {
    "AAC",
    "BIT188",
    "BIT12",
    "BIT21",
    "CTD",
    "GDPC",
    "GTPC",
    "IT",
    "OLP",
    "DPC",
    "ASDC",
    "GGAP",
    "PAAC",
    "CTF",
    "BIT20",
    "QSO",
}


# Features you already have from other toolchains.
EXISTING_FEATURES = {
    # Modlamp
    "LENGTH",
    "MW",
    "CHARGE_DENSITY",
    "NET_CHARGE",
    "AROMATICITY",
    "INSTABILITY_INDEX",
    "BOMAN_INDEX",
    "PI",
    "ALIPHATIC_INDEX",
    "HYDROPHOBIC_RATIO",
    # iFeature
    "MOLECULAR_WEIGHT",
    "ISOELECTRIC_POINT",
    "GRAVY",
    "HYDROPHOBICITY_PERCENTAGE",
    "NET_CHARGE_PH7",
    "HYDROPHILIC_RATIO",
    "EXTINCTION_COEFFICIENT_OXIDIZED",
    "EXTINCTION_COEFFICIENT_REDUCED",
    # pFeature
    "AAC",
    "DPC",
    "PAAC",
    "GAAC",
    "CKSAAP",
}


@dataclass(frozen=True)
class FeaturePlan:
    """Computed descriptor plan."""

    to_add: list[str]
    excluded_existing: list[str]
    paac_lambda: int

    @property
    def paac_dimension(self) -> int:
        """PAAC vector length for Chou's pseudo amino acid composition.

        Formula: 20 + 2*lambda
        """
        return 20 + 2 * self.paac_lambda


def load_sequences(file_path: str | Path, sequence_column: str = "sequence") -> pd.DataFrame:
    """Load Excel and uppercase sequence strings."""
    df = pd.read_excel(file_path)
    df[sequence_column] = df[sequence_column].astype(str).str.upper()
    return df


def build_feature_plan(
    requested: Iterable[str] = REQUESTED_FEATURES,
    existing: Iterable[str] = EXISTING_FEATURES,
    paac_lambda: int = 15,
    force_include_paac: bool = True,
) -> FeaturePlan:
    """Create final descriptor plan after excluding duplicates.

    Args:
        requested: candidate features to compute.
        existing: descriptors you already computed.
        paac_lambda: lambda value used in PAAC (dimension=20+2*lambda).
        force_include_paac: if True, keeps PAAC in `to_add` even if existing.
    """
    if paac_lambda < 1:
        raise ValueError("paac_lambda must be >= 1")

    requested_up = {x.upper().strip() for x in requested}
    existing_up = {x.upper().strip() for x in existing}

    excluded_existing = sorted(requested_up & existing_up)
    to_add = sorted(requested_up - existing_up)

    if force_include_paac and "PAAC" in requested_up and "PAAC" not in to_add:
        to_add.append("PAAC")
        to_add = sorted(set(to_add))

    return FeaturePlan(
        to_add=to_add,
        excluded_existing=excluded_existing,
        paac_lambda=paac_lambda,
    )


if __name__ == "__main__":
    # Update this path to your local Excel file.
    file_path = Path("combined.xlsx")

    if file_path.exists():
        data = load_sequences(file_path)
        print(data[["sequence"]].head())
    else:
        print("combined.xlsx not found in current directory; skipping data preview.")

    plan = build_feature_plan(paac_lambda=15, force_include_paac=True)
    print("\\nExcluded because already present:", plan.excluded_existing)
    print("Features to add:", plan.to_add)
    print(f"PAAC dimension (20 + 2λ): 20 + 2*{plan.paac_lambda} = {plan.paac_dimension}")
