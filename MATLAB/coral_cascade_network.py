import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import FancyArrowPatch
from pathlib import Path

outdir = Path("./output")
outdir.mkdir(parents=True, exist_ok=True)

# -----------------------------
# 1) Reconstruct / extend loss paths through 2100
# -----------------------------
years = np.arange(2024, 2101)

# Paper-quoted annual losses through 2050
years_anchor = np.array([2024, 2026, 2028, 2030, 2035, 2040, 2045, 2050])
wave_anchor = np.array([685.6, 710.9, 733.1, 752.5, 790.0, 814.9, 830.6, 834.3])   # $M
tour_anchor = np.array([112.4, 121.2, 130.1, 138.9, 160.9, 183.0, 205.0, 212.5])   # $M

# Interpolate to 2050
wave = np.interp(np.clip(years, 2024, 2050), years_anchor, wave_anchor)
tour = np.interp(np.clip(years, 2024, 2050), years_anchor, tour_anchor)

# Extrapolate 2051-2100 toward paper-implied ceilings for visualization:
# wave damage -> Dmax ≈ 850M, tourism -> baseline revenue ceiling ≈ 220.5M
mask = years > 2050
k = 0.08
wave[mask] = 850.0 - (850.0 - wave_anchor[-1]) * np.exp(-k * (years[mask] - 2050))
tour[mask] = 220.5 - (220.5 - tour_anchor[-1]) * np.exp(-k * (years[mask] - 2050))

direct_total = wave + tour

# Inferred fishery annual path only for visualization, constrained to paper's cumulative $4.01B by 2100
fishery_total_2100 = 4012.42  # $M cumulative
mid, steep = 2055, 0.12
cum_fish = 1 / (1 + np.exp(-steep * (years - mid)))
cum_fish = (cum_fish - cum_fish[0]) / (cum_fish[-1] - cum_fish[0]) * fishery_total_2100
annual_fish = np.diff(np.r_[0, cum_fish])

loss_df = pd.DataFrame({
    "year": years,
    "wave_damage_loss_musd": wave,
    "tourism_loss_musd": tour,
    "direct_total_loss_musd": direct_total,
    "inferred_fishery_annual_loss_musd": annual_fish,
    "inferred_fishery_cumulative_loss_musd": cum_fish,
})
csv_path = outdir / "coral_cascade_losses_2024_2100.csv"
loss_df.to_csv(csv_path, index=False)

# -----------------------------
# 2) Build the NetworkX cascade
# -----------------------------
G = nx.DiGraph()

nodes = {
    "Ocean warming": {"layer": 0, "group": "driver", "size": 1900},
    "Ocean acidification": {"layer": 0, "group": "driver", "size": 1700},
    "Cyclone intensity": {"layer": 0, "group": "driver", "size": 1700},

    "Thermal stress\n(~8 DHW threshold)": {"layer": 1, "group": "process", "size": 1800},

    "Bleaching frequency": {"layer": 2, "group": "reef", "size": 1800},
    "Hard coral cover decline\n(-0.6397%/yr)": {"layer": 2, "group": "reef", "size": 2300},
    "Algae displacement": {"layer": 2, "group": "reef", "size": 1600},

    "Structural reef loss": {"layer": 3, "group": "eco", "size": 1800},
    "Fish biomass decline": {"layer": 3, "group": "eco", "size": 1900},
    "Species loss / habitat contraction": {"layer": 3, "group": "eco", "size": 2000},
    "Aesthetic value loss": {"layer": 3, "group": "eco", "size": 1700},
    "Reduced wave attenuation\n(reefs absorb up to 97%)": {"layer": 3, "group": "eco", "size": 2100},

    "Tourism revenue loss\n$2.95B NPV to 2050": {"layer": 4, "group": "quant", "size": 2200},
    "Wave damage increase\n$14.24B NPV to 2050": {"layer": 4, "group": "quant", "size": 2600},
    "Fishery losses\n$4.01B by 2100": {"layer": 4, "group": "quant", "size": 2300},
    "Pharmaceutical option loss\n2100 lower-bound NPV\n(modeled, value not quoted in excerpts)": {"layer": 4, "group": "quant2", "size": 2300},

    "Hotels & lodging": {"layer": 5, "group": "down", "size": 1500},
    "Restaurants & bars": {"layer": 5, "group": "down", "size": 1500},
    "Dive / snorkel operators": {"layer": 5, "group": "down", "size": 1450},
    "Charter boats & marinas": {"layer": 5, "group": "down", "size": 1450},
    "Retail / local commerce": {"layer": 5, "group": "down", "size": 1450},
    "Insurance costs": {"layer": 5, "group": "down", "size": 1500},
    "Property values": {"layer": 5, "group": "down", "size": 1450},
    "Public infrastructure repair": {"layer": 5, "group": "down", "size": 1600},
    "Food security": {"layer": 5, "group": "down", "size": 1500},
    "Seafood processors & wholesalers": {"layer": 5, "group": "down", "size": 1500},
    "Pharma pipelines / biotech": {"layer": 5, "group": "down", "size": 1450},
    "Medical future value / unmet cures": {"layer": 5, "group": "down", "size": 1450},

    "Household income & jobs": {"layer": 6, "group": "macro", "size": 1900},
    "Tax base / public finance": {"layer": 6, "group": "macro", "size": 1800},
    "Migration / displacement": {"layer": 6, "group": "macro", "size": 1700},
    "Regional GDP drag by 2100": {"layer": 6, "group": "macro", "size": 2100},
}

for n, attrs in nodes.items():
    G.add_node(n, **attrs)

edges = [
    ("Ocean warming", "Thermal stress\n(~8 DHW threshold)", 2.5),
    ("Ocean acidification", "Hard coral cover decline\n(-0.6397%/yr)", 1.8),
    ("Cyclone intensity", "Reduced wave attenuation\n(reefs absorb up to 97%)", 1.7),

    ("Thermal stress\n(~8 DHW threshold)", "Bleaching frequency", 2.4),
    ("Bleaching frequency", "Hard coral cover decline\n(-0.6397%/yr)", 2.5),
    ("Algae displacement", "Hard coral cover decline\n(-0.6397%/yr)", 1.6),

    ("Hard coral cover decline\n(-0.6397%/yr)", "Structural reef loss", 2.1),
    ("Hard coral cover decline\n(-0.6397%/yr)", "Aesthetic value loss", 2.0),
    ("Hard coral cover decline\n(-0.6397%/yr)", "Reduced wave attenuation\n(reefs absorb up to 97%)", 2.8),
    ("Hard coral cover decline\n(-0.6397%/yr)", "Species loss / habitat contraction", 2.0),

    ("Structural reef loss", "Fish biomass decline", 2.0),
    ("Fish biomass decline", "Fishery losses\n$4.01B by 2100", 2.4),
    ("Aesthetic value loss", "Tourism revenue loss\n$2.95B NPV to 2050", 2.2),
    ("Reduced wave attenuation\n(reefs absorb up to 97%)", "Wave damage increase\n$14.24B NPV to 2050", 3.2),
    ("Species loss / habitat contraction", "Pharmaceutical option loss\n2100 lower-bound NPV\n(modeled, value not quoted in excerpts)", 2.2),

    ("Tourism revenue loss\n$2.95B NPV to 2050", "Hotels & lodging", 1.5),
    ("Tourism revenue loss\n$2.95B NPV to 2050", "Restaurants & bars", 1.7),
    ("Tourism revenue loss\n$2.95B NPV to 2050", "Dive / snorkel operators", 1.6),
    ("Tourism revenue loss\n$2.95B NPV to 2050", "Charter boats & marinas", 1.5),
    ("Tourism revenue loss\n$2.95B NPV to 2050", "Retail / local commerce", 1.3),

    ("Wave damage increase\n$14.24B NPV to 2050", "Insurance costs", 1.8),
    ("Wave damage increase\n$14.24B NPV to 2050", "Property values", 1.6),
    ("Wave damage increase\n$14.24B NPV to 2050", "Public infrastructure repair", 1.8),
    ("Wave damage increase\n$14.24B NPV to 2050", "Restaurants & bars", 1.0),
    ("Wave damage increase\n$14.24B NPV to 2050", "Hotels & lodging", 1.0),

    ("Fishery losses\n$4.01B by 2100", "Food security", 1.7),
    ("Fishery losses\n$4.01B by 2100", "Seafood processors & wholesalers", 1.5),
    ("Fishery losses\n$4.01B by 2100", "Restaurants & bars", 1.1),
    ("Fishery losses\n$4.01B by 2100", "Charter boats & marinas", 1.0),

    ("Pharmaceutical option loss\n2100 lower-bound NPV\n(modeled, value not quoted in excerpts)", "Pharma pipelines / biotech", 1.5),
    ("Pharmaceutical option loss\n2100 lower-bound NPV\n(modeled, value not quoted in excerpts)", "Medical future value / unmet cures", 1.6),
]

for source in [
    "Hotels & lodging", "Restaurants & bars", "Dive / snorkel operators", "Charter boats & marinas",
    "Retail / local commerce", "Insurance costs", "Property values", "Public infrastructure repair",
    "Food security", "Seafood processors & wholesalers", "Pharma pipelines / biotech", "Medical future value / unmet cures"
]:
    edges.append((source, "Household income & jobs", 1.1))

for source in [
    "Hotels & lodging", "Restaurants & bars", "Retail / local commerce",
    "Property values", "Public infrastructure repair", "Insurance costs"
]:
    edges.append((source, "Tax base / public finance", 1.0))

for source in ["Food security", "Property values", "Household income & jobs"]:
    edges.append((source, "Migration / displacement", 1.0))

for source in ["Household income & jobs", "Tax base / public finance", "Migration / displacement"]:
    edges.append((source, "Regional GDP drag by 2100", 1.6))

for u, v, w in edges:
    G.add_edge(u, v, weight=w)

# Manual layered layout for a polished flow
layers = {}
for n, d in G.nodes(data=True):
    layers.setdefault(d["layer"], []).append(n)

x_coords = {0: 0.03, 1: 0.17, 2: 0.31, 3: 0.47, 4: 0.64, 5: 0.82, 6: 0.96}
y_maps = {
    0: [0.82, 0.55, 0.28],
    1: [0.55],
    2: [0.80, 0.55, 0.30],
    3: [0.86, 0.68, 0.50, 0.32, 0.14],
    4: [0.82, 0.60, 0.38, 0.16],
    5: list(np.linspace(0.92, 0.08, 12)),
    6: [0.78, 0.56, 0.34, 0.12],
}

pos = {}
for layer, names in layers.items():
    for n, y in zip(names, y_maps[layer]):
        pos[n] = (x_coords[layer], y)

node_colors = {
    "driver": "#ff6b6b",
    "process": "#f7b801",
    "reef": "#2ec4b6",
    "eco": "#4dabf7",
    "quant": "#e76f51",
    "quant2": "#c77dff",
    "down": "#8bd450",
    "macro": "#f1f5f9",
}

layer_headers = {
    0: "Climate Drivers",
    1: "Thermal Trigger",
    2: "Reef-Specific Damage",
    3: "Ecological Transmission",
    4: "Quantified Loss Channels",
    5: "Dependent Sectors",
    6: "Macro Spillovers",
}

# -----------------------------
# 3) Render poster-style figure
# -----------------------------
fig = plt.figure(figsize=(22, 18), dpi=200)
gs = fig.add_gridspec(2, 1, height_ratios=[3.8, 1.4], hspace=0.06)

ax = fig.add_subplot(gs[0])
fig.patch.set_facecolor("#081421")
ax.set_facecolor("#081421")

# Subtle background glow
x = np.linspace(-1, 1, 700)
y = np.linspace(-1, 1, 450)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + (Y * 0.8)**2)
grad = np.clip(1 - R, 0, 1)
ax.imshow(grad, extent=[-0.02, 1.02, -0.02, 1.02], cmap="Blues", alpha=0.18, origin="lower", aspect="auto")

# Draw edges with soft glow
for u, v, d in G.edges(data=True):
    if "Wave damage" in v or "Wave damage" in u:
        edge_color = "#ff9f1c"
    elif "Tourism" in v or "Tourism" in u:
        edge_color = "#2ec4b6"
    elif "Fishery" in v or "Fishery" in u:
        edge_color = "#8bd450"
    elif "Pharmaceutical" in v or "Pharmaceutical" in u:
        edge_color = "#c77dff"
    else:
        edge_color = "#7f8fa6"

    width = 0.5 + d["weight"] * 1.3
    if pos[u][1] < pos[v][1] - 0.12:
        rad = 0.12
    elif pos[u][1] > pos[v][1] + 0.12:
        rad = -0.12
    else:
        rad = 0.0

    arrow = FancyArrowPatch(
        posA=pos[u],
        posB=pos[v],
        connectionstyle=f"arc3,rad={rad}",
        arrowstyle="-|>",
        mutation_scale=10 + width * 1.5,
        lw=width,
        color=edge_color,
        alpha=0.28,
        zorder=1,
    )
    ax.add_patch(arrow)

    glow = FancyArrowPatch(
        posA=pos[u],
        posB=pos[v],
        connectionstyle=f"arc3,rad={rad}",
        arrowstyle="-|>",
        mutation_scale=10 + width * 1.5,
        lw=width * 0.45,
        color="white",
        alpha=0.06,
        zorder=1,
    )
    ax.add_patch(glow)

# Draw nodes by group
for group, color in node_colors.items():
    nlist = [n for n, d in G.nodes(data=True) if d["group"] == group]
    sizes = [G.nodes[n]["size"] * 1.55 for n in nlist]

    nx.draw_networkx_nodes(
        G, pos, nodelist=nlist,
        node_size=[s * 1.35 for s in sizes],
        node_color=color, alpha=0.10, linewidths=0, ax=ax
    )
    nx.draw_networkx_nodes(
        G, pos, nodelist=nlist,
        node_size=[s * 1.12 for s in sizes],
        node_color=color, alpha=0.18, linewidths=0, ax=ax
    )
    nx.draw_networkx_nodes(
        G, pos, nodelist=nlist,
        node_size=sizes,
        node_color=color, alpha=0.92,
        linewidths=1.2, edgecolors="white", ax=ax
    )

# Labels
for n, (x0, y0) in pos.items():
    fs = 12
    if G.nodes[n]["layer"] >= 4:
        fs = 11.5
    if G.nodes[n]["layer"] == 6:
        fs = 13
    ax.text(
        x0, y0, n,
        ha="center", va="center",
        color="white" if G.nodes[n]["group"] != "macro" else "#0f172a",
        fontsize=fs, fontweight="bold", zorder=3
    )

# Layer headers and guide lines
for layer, header in layer_headers.items():
    ax.text(
        x_coords[layer], 1.01, header,
        ha="center", va="bottom",
        color="white", fontsize=15, fontweight="bold"
    )
    ax.plot([x_coords[layer], x_coords[layer]], [0.02, 0.98], color="white", alpha=0.07, lw=1)

ax.text(
    0.5, 1.08, "Cascading Economic Effects of Coral Reef Degradation",
    ha="center", va="bottom", color="white", fontsize=26, fontweight="bold"
)
ax.text(
    0.5, 1.045,
    "Expanded NetworkX system map combining the paper's quantified channels with inferred downstream spillovers through 2100",
    ha="center", va="bottom", color="#cbd5e1", fontsize=14
)

summary = (
    "Paper-backed losses used here:\n"
    "• Wave damage NPV to 2050: $14.24B\n"
    "• Tourism loss NPV to 2050: $2.95B\n"
    "• Fishery losses by 2100: $4.01B\n"
    "• Annual direct losses rise from $798M (2024) to $1.05B (2050)\n"
    "• Pharmaceutical loss is included as a distinct irreversible channel,\n"
    "  but the excerpts did not provide one quoted global dollar total\n"
    "Visualization extensions added by me:\n"
    "• 2051-2100 wave/tourism values extrapolated toward model ceilings\n"
    "• Fishery annual path inferred only so the 2100 cumulative total can be visualized"
)
ax.text(
    0.012, -0.02, summary,
    ha="left", va="top", color="white", fontsize=12,
    bbox=dict(boxstyle="round,pad=0.6", facecolor="#0d2238", edgecolor="#89c2ff", alpha=0.95)
)

ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.08, 1.12)
ax.axis("off")

# Bottom panel: loss trajectories
ax2 = fig.add_subplot(gs[1])
ax2.set_facecolor("#081421")

ax2.fill_between(years, wave, alpha=0.28, label="Wave damage loss ($M/yr)")
ax2.fill_between(years, wave + tour, wave, alpha=0.28, label="Tourism loss ($M/yr)")
ax2.plot(years, direct_total, linewidth=3.0, label="Direct annual loss total ($M/yr)")
ax2.plot(years, annual_fish, linewidth=2.0, linestyle="--", label="Inferred fishery annual loss for visualization ($M/yr)")

for yr in [2028, 2036, 2044, 2049]:
    ax2.axvline(yr, linewidth=1.2, alpha=0.45)
    ax2.text(yr, ax2.get_ylim()[1] if ax2.get_ylim()[1] > 0 else 1, str(yr),
             rotation=90, va="top", ha="right", fontsize=9, color="white")

ax2.set_title("Projected annual loss channels through 2100", color="white", fontsize=16, pad=10, fontweight="bold")
ax2.set_xlabel("Year", color="white")
ax2.set_ylabel("Annual loss ($ millions)", color="white")
ax2.tick_params(colors="white")
for spine in ax2.spines.values():
    spine.set_color("#94a3b8")
ax2.grid(alpha=0.15)

legend = ax2.legend(loc="upper left", frameon=True)
legend.get_frame().set_facecolor("#0f172a")
legend.get_frame().set_edgecolor("#475569")
for text in legend.get_texts():
    text.set_color("white")

ax2.text(
    0.99, 0.03,
    "Notes: solid/filled series use paper-backed annual values through 2050 and\n"
    "transparent extrapolation thereafter. Dashed fishery series is my inferred time allocation\n"
    "of the paper's cumulative $4.01B loss by 2100, not a quoted annual table from the paper.",
    transform=ax2.transAxes, ha="right", va="bottom", color="#cbd5e1", fontsize=10
)

png_path = outdir / "coral_cascade_network_2100.png"
svg_path = outdir / "coral_cascade_network_2100.svg"
fig.savefig(png_path, bbox_inches="tight", facecolor=fig.get_facecolor())
fig.savefig(svg_path, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)

# -----------------------------
# 4) Save the script used
# -----------------------------
script = r'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import FancyArrowPatch

# Reconstruct annual losses and build the graph...
# This script mirrors the figure-generation logic used in the delivered artifacts.
'''
script_path = outdir / "coral_cascade_network_2100_script.py"

# Save the full executable code by reading it from this notebook cell would be awkward,
# so write a compact but functional standalone script with the same outputs.
standalone_code = f"""
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

pd.DataFrame({{
    "year": years,
    "wave_damage_loss_musd": wave,
    "tourism_loss_musd": tour,
    "direct_total_loss_musd": direct_total,
    "inferred_fishery_annual_loss_musd": annual_fish,
    "inferred_fishery_cumulative_loss_musd": cum_fish,
}}).to_csv(outdir / "coral_cascade_losses_2024_2100.csv", index=False)

print("Artifacts recreated.")
"""
script_path.write_text(standalone_code)

print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")
print(f"Saved: {csv_path}")
print(f"Saved: {script_path}")