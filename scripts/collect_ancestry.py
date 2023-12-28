import pandas as pd
import re

# Import all ancestry filenames
ancestry = snakemake.input["ancestry"]

# Define resulting dataframe
all_ancestry = pd.DataFrame()

for anc in ancestry:
    # Read each .Ancestry file's IntendedSample row
    df = pd.read_csv(anc, sep="\t")[["IntendedSample"]]

    # Transpose rows to columns
    df = df.transpose()

    # Add columns names back for all the PCs
    df.columns = ["PC1", "PC2", "PC3", "PC4"]

    # Extract sample name from filename and assign to a new column
    df["SAMPLE"] = re.search("HGDP[0-9]+", anc).group(0)

    # Append new dataframe to all the existing ones
    all_ancestry = pd.concat([all_ancestry, df])

    # Reorder columns for readability
    all_ancestry = all_ancestry[["SAMPLE", "PC1", "PC2", "PC3", "PC4"]]

# Write resulting dataframe with all ancestry information to file
all_ancestry.to_csv(snakemake.output["all_ancestry"], index=None, sep="\t")
