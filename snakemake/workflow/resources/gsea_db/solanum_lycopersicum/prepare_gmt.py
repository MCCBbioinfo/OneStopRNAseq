# Tomato gene sets are retrieve from: https://systemsbiology.cau.edu.cn/PlantGSEAv2/download.php

#%% 
import csv
import os
from collections import defaultdict

#%%
# Input file
input_file = "Sly_ALL.txt"   # change to your filename
output_dir = "gsea_db"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Data structure:
# category -> gene_set_name -> set(genes)
data = defaultdict(lambda: defaultdict(set))

#%%
# Read file (tab-delimited)
with open(input_file, "r", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    if not reader.fieldnames:
        raise ValueError(f"No header found in input file: {input_file}")

    normalized_to_actual = {}
    for name in reader.fieldnames:
        if name is None:
            continue
        normalized_to_actual[name.strip().lower()] = name

    required_columns = {
        "gene": "gene",
        "gene_set": "gene_set",
        "category": "category",
    }

    resolved_columns = {}
    missing = []
    for normalized, display in required_columns.items():
        actual = normalized_to_actual.get(normalized)
        if actual is None:
            missing.append(display)
        else:
            resolved_columns[normalized] = actual

    if missing:
        raise KeyError(
            "Missing required column(s) in input file header: "
            + ", ".join(missing)
            + f". Found columns: {reader.fieldnames}"
        )

    for row in reader:
        gene = (row.get(resolved_columns["gene"]) or "").strip()
        gene_set = (row.get(resolved_columns["gene_set"]) or "").strip()
        category = (row.get(resolved_columns["category"]) or "").strip()
        if not gene or not gene_set or not category:
            continue
        data[category][gene_set].add(gene)

#%%
# Write GMT files (one per Category)
for category, gene_sets in data.items():

    out_file = os.path.join(output_dir, f"{category}.gmt")

    with open(out_file, "w") as out:

        for gene_set, genes in gene_sets.items():

            # GMT line:
            # name  description  gene1 gene2 ...
            line = [
                gene_set,
                category,
                *sorted(genes)
            ]

            out.write("\t".join(line) + "\n")

    print(f"Written: {out_file}")

# %%
