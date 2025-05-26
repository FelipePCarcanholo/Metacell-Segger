
#how to use it
#conda environment cellxgene
#will have to change this pathways: single_cell_path, metacell_output_path, name_column_cell_type
#the single cell must have quality control before running here
#the gene list must be in SYMBOL not ENSEMBL
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import argparse
import anndata as ad             # For reading/writing AnnData files
import matplotlib.pyplot as plt  # For plotting
import metacells as mc           # The Metacells package
import numpy as np               # For array/matrix operations
import pandas as pd              # For data frames
import os                        # For filesystem operations
import seaborn as sb             # For plotting
import scipy.sparse as sp        # For sparse matrices
import shutil                    # for filesystem operations
from math import hypot           # For plotting

# Parse command line arguments
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-sc", "--single_cell_path", help="Path to single cell anndata")
parser.add_argument("-mtc", "--metacell_output_path", help="Path to metacell output")
parser.add_argument("-cn", "--name_column_cell_type", default="cell_type", help="Name of the collumn to pass the information, defoult: cell_type")
parser.add_argument("-sz", "--metacell_size", default=48, type=int, help="Name of the collumn to pass the information, defoult: cell_type")

args = vars(parser.parse_args())
 
# Set up parameters
single_cell_path = args["single_cell_path"]
metacell_output_path = args["metacell_output_path"]
name_column_cell_type = args["name_column_cell_type"]
metacell_size = args["metacell_size"]


full = ad.read_h5ad(single_cell_path)
full.var.index = full.var.index.astype(str)
#full.raw.var.index = full.raw.var.index.astype(str)
#put in the rownames
full.var["ENSEMBL"] = full.var.index
id_to_symbol = dict(zip(full.var['ENSEMBL'], full.var['feature_name']))

# full = gbm_sc[:,(result_df['Ensembl_ID'].tolist())] #leaving just the genes that has a convertion to gene symbol, without converting yet
# Update the var index using the dictionary and filter out missing genes
full.var.index = full.var.index.map(id_to_symbol)

#full.raw.var["ENSEMBL"] = full.raw.var.index
#id_to_symbol = dict(zip(full.raw.var['ENSEMBL'], full.raw.var['feature_name']))
#full.raw.var.index = full.raw.var.index.map(id_to_symbol)
full.var_names_make_unique()
mc.ut.top_level(full)
mc.ut.set_name(full, "sc_BRCA_atlas")
PROPERLY_SAMPLED_MIN_CELL_TOTAL = 200
PROPERLY_SAMPLED_MAX_CELL_TOTAL = 20000
total_umis_per_cell = mc.ut.get_o_numpy(full, "__x__", sum=True)
plot = sb.displot(total_umis_per_cell, log_scale=(10, None))
plot.set(xlabel="UMIs", ylabel="Density", yticks=[])

plot.refline(x=PROPERLY_SAMPLED_MIN_CELL_TOTAL, color="darkgreen")
plot.refline(x=PROPERLY_SAMPLED_MAX_CELL_TOTAL, color="crimson")


too_small_cells_count = np.sum(total_umis_per_cell < PROPERLY_SAMPLED_MIN_CELL_TOTAL)
too_large_cells_count = np.sum(total_umis_per_cell > PROPERLY_SAMPLED_MAX_CELL_TOTAL)

total_umis_per_cell = mc.ut.get_o_numpy(full, name="__x__", sum=True)
too_small_cells_percent = 100.0 * too_small_cells_count / full.n_obs
too_large_cells_percent = 100.0 * too_large_cells_count / full.n_obs

print(
    f"Will exclude {too_small_cells_count} ({too_small_cells_percent:.2f}%%) cells"
    f" with less than {PROPERLY_SAMPLED_MIN_CELL_TOTAL} UMIs"
)
print(
    f"Will exclude {too_large_cells_count} ({too_large_cells_percent:.2f}%%) cells"
    f" with more than {PROPERLY_SAMPLED_MAX_CELL_TOTAL} UMIs"
)
EXCLUDED_GENE_NAMES = [
    "XIST", "MALAT1",   # Sex-specific genes.
    "NEAT1"             # Non-coding.
]
EXCLUDED_GENE_PATTERNS = ["MT-.*"]  # Mytochondrial.
if full.X.dtype != np.float32:
    print(f"Converting data matrix from {full.X.dtype} to float32...")
    full.X = full.X.astype(np.float32)
print(type(full.X), full.X.dtype)
mc.pl.exclude_genes(
    full,
    excluded_gene_names=EXCLUDED_GENE_NAMES, 
    excluded_gene_patterns=EXCLUDED_GENE_PATTERNS,
    random_seed=123456,
)
mc.tl.compute_excluded_gene_umis(full)
PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION = 0.25
excluded_umis_fraction_regularization = 1e-3  # Avoid 0 values in log scale plot.
excluded_umis_per_cell = mc.ut.get_o_numpy(full, "excluded_umis")
excluded_umis_fraction_per_cell = excluded_umis_per_cell / total_umis_per_cell

excluded_umis_fraction_per_cell += excluded_umis_fraction_regularization
plot = sb.displot(excluded_umis_fraction_per_cell, log_scale=(10, None))
excluded_umis_fraction_per_cell -= excluded_umis_fraction_regularization

plot.set(xlabel="Fraction of excluded gene UMIs", ylabel="Density", yticks=[])
plot.refline(x=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION, color="crimson")

too_excluded_cells_count = np.sum(
    excluded_umis_fraction_per_cell > PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION
)
too_excluded_cells_fraction = too_excluded_cells_count / full.n_obs

print(
    f"Will exclude {too_excluded_cells_count} ({100 * too_excluded_cells_fraction:.2f}%) cells"
    f" with more than {100 * PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION:.2f}% excluded gene UMIs"
)


mc.pl.exclude_cells(
    full,
    properly_sampled_min_cell_total=PROPERLY_SAMPLED_MIN_CELL_TOTAL,
    properly_sampled_max_cell_total=PROPERLY_SAMPLED_MAX_CELL_TOTAL,
    properly_sampled_max_excluded_genes_fraction=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION,
    #additional_cells_masks=["|doublet_cell"]
)


clean = mc.pl.extract_clean_data(full, name="sc_BRCA_atlas.clean")
mc.ut.top_level(clean)
cells = clean
LATERAL_GENE_NAMES = [
    "ACSM3", "ANP32B", "APOE", "AURKA", "B2M", "BIRC5", "BTG2", "CALM1", "CD63", "CD69", "CDK4",
    "CENPF", "CENPU", "CENPW", "CH17-373J23.1", "CKS1B", "CKS2", "COX4I1", "CXCR4", "DNAJB1",
    "DONSON", "DUSP1", "DUT", "EEF1A1", "EEF1B2", "EIF3E", "EMP3", "FKBP4", "FOS", "FOSB", "FTH1",
    "G0S2", "GGH", "GLTSCR2", "GMNN", "GNB2L1", "GPR183", "H2AFZ", "H3F3B", "HBM", "HIST1H1C",
    "HIST1H2AC", "HIST1H2BG", "HIST1H4C", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB",
    "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-E", "HLA-F", "HMGA1",
    "HMGB1", "HMGB2", "HMGB3", "HMGN2", "HNRNPAB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B",
    "HSPA6", "HSPD1", "HSPE1", "HSPH1", "ID2", "IER2", "IGHA1", "IGHA2", "IGHD", "IGHG1", "IGHG2",
    "IGHG3", "IGHG4", "IGHM", "IGKC", "IGKV1-12", "IGKV1-39", "IGKV1-5", "IGKV3-15", "IGKV4-1",
    "IGLC2", "IGLC3", "IGLC6", "IGLC7", "IGLL1", "IGLL5", "IGLV2-34", "JUN", "JUNB", "KIAA0101",
    "LEPROTL1", "LGALS1", "LINC01206", "LTB", "MCM3", "MCM4", "MCM7", "MKI67", "MT2A", "MYL12A",
    "MYL6", "NASP", "NFKBIA", "NUSAP1", "PA2G4", "PCNA", "PDLIM1", "PLK3", "PPP1R15A", "PTMA",
    "PTTG1", "RAN", "RANBP1", "RGCC", "RGS1", "RGS2", "RGS3", "RP11-1143G9.4", "RP11-160E2.6",
    "RP11-53B5.1", "RP11-620J15.3", "RP5-1025A1.3", "RP5-1171I10.5", "RPS10", "RPS10-NUDT3", "RPS11",
    "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS19BP1",
    "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L",
    "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6", "RPS6KA1",
    "RPS6KA2", "RPS6KA2-AS1", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2",
    "RPS6KC1", "RPS6KL1", "RPS7", "RPS8", "RPS9", "RPSA", "RRM2", "SMC4", "SRGN", "SRSF7", "STMN1",
    "TK1", "TMSB4X", "TOP2A", "TPX2", "TSC22D3", "TUBA1A", "TUBA1B", "TUBB", "TUBB4B", "TXN", "TYMS",
    "UBA52", "UBC", "UBE2C", "UHRF1", "YBX1", "YPEL5", "ZFP36", "ZWINT"
]
LATERAL_GENE_PATTERNS = ["RP[LS].*"]  # Ribosomal
# This will mark as "lateral_gene" any genes that match the above, if they exist in the clean dataset.
mc.pl.mark_lateral_genes(
    cells,
    lateral_gene_names=LATERAL_GENE_NAMES,
    lateral_gene_patterns=LATERAL_GENE_PATTERNS,
)

lateral_gene_mask = mc.ut.get_v_numpy(cells, "lateral_gene")
lateral_gene_names = set(cells.var_names[lateral_gene_mask])
print(sorted([
    name for name in lateral_gene_names
    if not name.startswith("RPL") and not name.startswith("RPS")
]))
print(f"""and {len([
    name for name in lateral_gene_names if name.startswith("RPL") or name.startswith("RPS")
])} RP[LS].* genes""")
NOISY_GENE_NAMES = [
    "CCL3", "CCL4", "CCL5", "CXCL8", "DUSP1", "FOS", "G0S2", "HBB", "HIST1H4C", "IER2", "IGKC",
    "IGLC2", "JUN", "JUNB", "KLRB1", "MT2A", "RPS26", "RPS4Y1", "TRBC1", "TUBA1B", "TUBB"
]
mc.pl.mark_noisy_genes(cells, noisy_gene_names=NOISY_GENE_NAMES)
# Either use the guesstimator:
max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
# Or, if running out of memory manually override:
# max_paralle_piles = ...
print(max_parallel_piles)
mc.pl.set_max_parallel_piles(max_parallel_piles)
with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(cells, random_seed=123456, target_metacell_size=metacell_size)

metacells = \
    mc.pl.collect_metacells(cells, name="sc_BRCA_atlas.metacells", random_seed=123456)
print(f"Preliminary: {metacells.n_obs} metacells, {metacells.n_vars} genes")
mc.tl.convey_obs_to_group(
    adata=cells, gdata=metacells,
    property_name=name_column_cell_type, to_property_name="cell_type",
    method=mc.ut.most_frequent  # This is the default, for categorical data
)
cells_meta = cells.obs
cells_meta = cells_meta[[name_column_cell_type, "metacell_name"]]
counts = cells_meta.groupby(["metacell_name", name_column_cell_type]).size()
most_cell_type = counts.unstack(fill_value=0)
most_cell_type['total_cells'] = most_cell_type.sum(axis=1)
#getting the % of cells
normalized_counts = most_cell_type.drop(columns='total_cells')
normalized_counts = normalized_counts.div(most_cell_type['total_cells'], axis=0)
most_cell_type[normalized_counts.columns] = normalized_counts
path_most_cell_type = os.path.join(metacell_output_path, "most_cell_type.csv")
most_cell_type.to_csv(path_most_cell_type)
metacells.obs["porct_cell_type"] = None
metacells.obs["total_cells"] = None
for line in range(len(metacells.obs)):
    cell_type = metacells.obs["cell_type"].iloc[line]
    metacell_name = metacells.obs.index[line]
    metacells.obs["porct_cell_type"].iloc[line] = most_cell_type.loc[metacell_name, cell_type]
    metacells.obs["total_cells"].iloc[line] = most_cell_type["total_cells"].loc[metacell_name]

if 'porct_cell_type' in metacells.obs.columns:
    metacells.obs['porct_cell_type'] = metacells.obs['porct_cell_type'].astype(str)
metacells.obs['porct_cell_type'] = metacells.obs['porct_cell_type'].fillna("NaN").astype(str)
if 'total_cells' in metacells.obs.columns:
    metacells.obs['total_cells'] = metacells.obs['total_cells'].astype(str)
metacells.obs['total_cells'] = metacells.obs['total_cells'].fillna("NaN").astype(str)
path_metacell = os.path.join(metacell_output_path, "metacell_output.h5ad")
metacells.write_h5ad(path_metacell)
path_cells = os.path.join(metacell_output_path, "sc_with_mtc_annotation.h5ad")
cells.write_h5ad(path_cells)




