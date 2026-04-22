#%%
import sys 
import pandas as pd

meta = sys.argv[1]
contrast = sys.argv[2]
contrast_num = int(sys.argv[3]) - 1
condition = sys.argv[4]
expr = sys.argv[5]
out_condition = sys.argv[6]
out_expr = sys.argv[7]

# meta = "meta.csv"
# contrast = "contrast.de.csv"
# contrast_num = 2 - 1
# condition = "original_condition.csv"
# expr = "original_EXPR.csv"
# out_condition = "condition.csv"
# out_expr = "EXPR.csv"

#%%
# Read meta data and contrast definitions
meta_df = pd.read_csv(meta)
contrast_df = pd.read_csv(contrast)
condition_df = pd.read_csv(condition)
expr_df = pd.read_csv(expr)
# Parse treatment (first line) and control (second line)
# Split by ';' to allow multiple groups
treatment_groups = set(contrast_df.iloc[0, contrast_num].split(';'))
control_groups = set(contrast_df.iloc[1, contrast_num].split(';'))
# if the two groups overlap, warning and exit
if len(treatment_groups & control_groups) > 0:
    raise ValueError("Treatment and control groups overlap!")

# Function to map group to condition
def assign_condition(group):
    if group in control_groups:
        return "control"
    elif group in treatment_groups:
        return "treatment"
    else:
        return "NA"

meta_df["condition"] = meta_df["GROUP_LABEL"].apply(assign_condition)
# Save only the two requested columns excluding rows where condition is "NA"
meta_df = meta_df[["SAMPLE_LABEL", "condition"]].rename(columns={"SAMPLE_LABEL": "SampleID"})
meta_df = meta_df[meta_df["condition"] != "NA"] 
meta_df.to_csv(out_condition, index=False)

# Subset expr_df to only keep and reorder columns based on meta_df's SampleID while always keeping the first column
expr_df = expr_df[["TE"] + meta_df["SampleID"].tolist()]
expr_df.to_csv(out_expr, index=False)
# %%
