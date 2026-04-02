
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import FancyArrowPatch
from pathlib import Path

outdir = Path(".")
years = np.arange(2024, 2101)
years_anchor = np.array([2024, 2026, 2028, 2030, 2035, 2040, 2045, 2050])
wave_anchor = np.array([685.6, 710.9, 733.1, 752.5, 790.0, 814.9, 830.6, 834.3])
tour_anchor = np.array([112.4, 121.2, 130.1, 138.9, 160.9, 183.0, 205.0, 212.5])

wave = np.interp(np.clip(years, 2024, 2050), years_anchor, wave_anchor)
tour = np.interp(np.clip(years, 2024, 2050), years_anchor, tour_anchor)

mask = years > 2050
k = 0.08
wave[mask] = 850.0 - (850.0 - wave_anchor[-1]) * np.exp(-k * (years[mask] - 2050))
tour[mask] = 220.5 - (220.5 - tour_anchor[-1]) * np.exp(-k * (years[mask] - 2050))
direct_total = wave + tour

fishery_total_2100 = 4012.42
mid, steep = 2055, 0.12
cum_fish = 1 / (1 + np.exp(-steep * (years - mid)))
cum_fish = (cum_fish - cum_fish[0]) / (cum_fish[-1] - cum_fish[0]) * fishery_total_2100
annual_fish = np.diff(np.r_[0, cum_fish])

pd.DataFrame({
    "year": years,
    "wave_damage_loss_musd": wave,
    "tourism_loss_musd": tour,
    "direct_total_loss_musd": direct_total,
    "inferred_fishery_annual_loss_musd": annual_fish,
    "inferred_fishery_cumulative_loss_musd": cum_fish,
}).to_csv(outdir / "coral_cascade_losses_2024_2100.csv", index=False)

print("Artifacts recreated.")
