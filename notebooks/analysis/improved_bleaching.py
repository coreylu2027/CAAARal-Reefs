# %%
# -----------------------------
# 0) Imports + Settings
# -----------------------------
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.metrics import (
    mean_absolute_error, mean_squared_error, r2_score,
    classification_report, confusion_matrix, ConfusionMatrixDisplay
)

from sklearn.linear_model import LinearRegression, Ridge, LogisticRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GridSearchCV

pd.set_option("display.max_columns", 200)
np.random.seed(42)

# %%
# -----------------------------
# 1) Load Data
# -----------------------------
# Update this path to your CSV
DATA_PATH = "../../data/processed/NOAA/NOAA_Global_Bleaching_Cleaned.csv"

df = pd.read_csv(DATA_PATH)
print(df.shape)
df

# %%
# -----------------------------
# 2) Basic cleaning
# -----------------------------
# Strip whitespace from column names
df.columns = [c.strip() for c in df.columns]

# Parse Date if present
if "Date" in df.columns:
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")

# Quick glance at missingness
missing = df.isna().mean().sort_values(ascending=False)
missing.head(20)

# %%
# -----------------------------
# 3) Define Targets + Features
# -----------------------------
# Primary regression target
TARGET_REG = "Percent_Bleaching"

# Optional classification target
TARGET_CLS = "Bleaching_Level"

# Columns you typically should NOT use as predictors (free-text + comments)
drop_text = [c for c in df.columns if "Comments" in c] + ["Site_Comments", "Sample_Comments", "Bleaching_Comments"]
drop_text = [c for c in drop_text if c in df.columns]

# You may or may not want to drop Date (keep it if you do time-based analyses)
# For baseline modeling, we’ll turn Date into year/month features (optional).
use_date_features = "Date" in df.columns

# Keep a copy
data = df.copy()

# Ensure target is numeric so it remains in corr
if "Percent Bleaching" in data.columns and "Percent_Bleaching" not in data.columns:
    data = data.rename(columns={"Percent Bleaching": "Percent_Bleaching"})
if "Percent_Bleaching" in data.columns:
    data["Percent_Bleaching"] = pd.to_numeric(data["Percent_Bleaching"], errors="coerce")

# Create simple date features (optional)
if use_date_features:
    data["Year"] = data["Date"].dt.year
    data["Month"] = data["Date"].dt.month

# Candidate predictors (exclude targets and text)
exclude = set([TARGET_REG, TARGET_CLS, "Date"]) | set(drop_text)
feature_cols = [c for c in data.columns if c not in exclude]

# If you want to remove raw lat/long from non-spatial models, uncomment:
# Normalize mixed-type feature columns to avoid OneHotEncoder errors
for c in feature_cols:
    if data[c].dtype == "object":
        as_num = pd.to_numeric(data[c], errors="coerce")
        if as_num.notna().mean() >= 0.9:
            data[c] = as_num
        else:
            data[c] = data[c].astype("string")

# feature_cols = [c for c in feature_cols if c not in ["Latitude_Degrees", "Longitude_Degrees"]]

print("Num features:", len(feature_cols))
feature_cols[:20]

# %%
# -----------------------------
# 4) Quick EDA
# -----------------------------
# Distribution of Percent_Bleaching
if TARGET_REG in data.columns:
    plt.figure()
    data[TARGET_REG].dropna().hist(bins=30)
    plt.title("Percent_Bleaching distribution")
    plt.xlabel("Percent_Bleaching")
    plt.ylabel("Count")
    plt.show()

# If Bleaching_Level exists, show counts
if TARGET_CLS in data.columns:
    plt.figure()
    data[TARGET_CLS].astype("category").value_counts(dropna=False).plot(kind="bar")
    plt.title("Bleaching_Level counts")
    plt.xlabel("Bleaching_Level")
    plt.ylabel("Count")
    plt.show()

# %%
# -----------------------------
# 5) Correlation (numeric only) — helps detect redundancy
# -----------------------------
numeric_cols = [c for c in feature_cols if pd.api.types.is_numeric_dtype(data[c])]
corr = data[numeric_cols + ([TARGET_REG] if TARGET_REG in data.columns else [])].corr(numeric_only=True)

# Show top correlations with target
if TARGET_REG in data.columns and TARGET_REG in corr.columns:
    target_corr = corr[TARGET_REG].drop(TARGET_REG).sort_values(key=lambda s: s.abs(), ascending=False)
    target_corr.head(15)
else:
    print("Target not found for correlation analysis.")

target_corr

# %%
# -----------------------------
# 6) Train/Test Split (Regression)
# -----------------------------
# Keep only rows with target for regression
reg_df = data.dropna(subset=[TARGET_REG]).copy()

X = reg_df[feature_cols]
y = reg_df[TARGET_REG]

# Identify categorical vs numeric predictors
cat_cols = [c for c in feature_cols if (not pd.api.types.is_numeric_dtype(reg_df[c]))]
num_cols = [c for c in feature_cols if c not in cat_cols]

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

len(X_train), len(X_test)
# %%
# -----------------------------
# 7) Preprocessing Pipeline
# -----------------------------
numeric_transformer = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler())
])

categorical_transformer_sparse = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="most_frequent")),
    ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=True))
])

categorical_transformer_dense = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="most_frequent")),
    ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False))
])

preprocess_sparse = ColumnTransformer(
    transformers=[
        ("num", numeric_transformer, num_cols),
        ("cat", categorical_transformer_sparse, cat_cols)
    ],
    remainder="drop"
)

preprocess_dense = ColumnTransformer(
    transformers=[
        ("num", numeric_transformer, num_cols),
        ("cat", categorical_transformer_dense, cat_cols)
    ],
    remainder="drop",
    sparse_threshold=0.0  # force dense output for HistGradientBoosting
)

# Default preprocessor used by sparse-compatible models
preprocess = preprocess_sparse

preprocess

# %%
# Normalize mixed-type feature columns
for c in feature_cols:
    if data[c].dtype == "object":
        as_num = pd.to_numeric(data[c], errors="coerce")
        if as_num.notna().mean() >= 0.9:
            data[c] = as_num
        else:
            data[c] = data[c].astype("string")
# %%
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

models = {
    "LinearRegression": LinearRegression(),
    "Ridge": Ridge(),
    # "RandomForestRegressor": RandomForestRegressor(
    #     n_estimators=400,
    #     random_state=42,
    #     n_jobs=-1,
    #     min_samples_leaf=2
    # ),
    "HistGradientBoosting": HistGradientBoostingRegressor(random_state=42)
}

fitted_models = {}
results = []

for name, model in models.items():
    # HistGradientBoosting requires dense input; others can use sparse
    preprocess_for_model = preprocess_dense if name == "HistGradientBoosting" else preprocess_sparse

    pipe = Pipeline(steps=[
        ("preprocess", preprocess_for_model),
        ("model", model)
    ])

    pipe.fit(X_train, y_train)

    # Save the fitted pipeline
    fitted_models[name] = pipe

    pred = pipe.predict(X_test)

    mae = mean_absolute_error(y_test, pred)
    mse = mean_squared_error(y_test, pred)
    rmse = mse ** 0.5
    r2 = r2_score(y_test, pred)

    results.append((name, mae, rmse, r2))

res_df = (
    pd.DataFrame(results, columns=["model", "MAE", "RMSE", "R2"])
    .sort_values("RMSE")
)

# %%
res_df
# %%
param_grid = {
    "model__learning_rate": [0.01, 0.05, 0.1],
    "model__max_depth": [None, 5, 10],
    "model__max_iter": [200, 400, 800]
}

preprocess = preprocess_dense.set_params(sparse_threshold=0.0)

pipe = Pipeline([
    ("preprocess", preprocess),
    ("model", HistGradientBoostingRegressor(random_state=42))
], verbose=True
)

grid = GridSearchCV(
    pipe,
    param_grid,
    cv=5,
    scoring="neg_root_mean_squared_error",
    n_jobs=-1
)

grid.fit(X, y)

print("Best params:", grid.best_params_)
print("Best RMSE:", -grid.best_score_)

# %%
