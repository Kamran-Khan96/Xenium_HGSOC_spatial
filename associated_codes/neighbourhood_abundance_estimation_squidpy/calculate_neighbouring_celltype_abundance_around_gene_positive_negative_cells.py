import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------
# Load and filter edges: tumour source → non-tumour neighbours ## calculated neighbouring celltype counts for the all cells across the TMAs
# ---------------------------------------------------------

tma5_cell_neigh = pd.read_csv("tma5_sqpy_all_cell_cell_touchx_results_fine_labels_cell_id_wise.csv")

tma6_cell_neigh = pd.read_csv("tma6_sqpy_all_cell_cell_touchx_results_fine_labels_cell_id_wise.csv")

df_edges = pd.concat([tma5_cell_neigh, tma6_cell_neigh], ignore_index=True)

df_edges = df_edges[
    (df_edges["src_type"] == "Tumour") &
    (df_edges["neighbor_type"] != "Tumour")
].copy()

# ---------------------------------------------------------
# Load gene status table (cell-id wise gene-positive/negative status)
# ---------------------------------------------------------

df_status = pd.read_csv("tma5_6_cell_id_wise_gene_exp_status.csv")


df_status_long = df_status.melt(
    id_vars="cell_id",
    var_name="Gene_Name",
    value_name="Gene_Status"
).dropna()

# ---------------------------------------------------------
# Merge edges × gene status (keeps +ve and -ve)
# ---------------------------------------------------------
df_merged = df_edges.merge(df_status_long, on="cell_id", how="inner")

# ---------------------------------------------------------
# Aggregate counts & unstack into a wide matrix
# ---------------------------------------------------------
df_wide = (
    df_merged
    .groupby(["radius_um", "PatientID", "Gene_Name", "Gene_Status", "neighbor_type"])["tgt_count"]
    .sum()
    .unstack(fill_value=0)
    .reset_index()
)

# ---------------------------------------------------------
# Normalise to 100% fractions
# ---------------------------------------------------------
print("Normalizing to 100% non-malignant microenvironment fractions...")
neighbor_cols = sorted(
    c for c in df_wide.columns
    if c not in ["radius_um", "PatientID", "Gene_Name", "Gene_Status"]
)

row_totals = df_wide[neighbor_cols].sum(axis=1).replace(0, 1)
df_wide[neighbor_cols] = df_wide[neighbor_cols].div(row_totals, axis=0) * 100

df_wide = df_wide.rename(columns={"radius_um": "Radius_um"})
