import matplotlib.pyplot as plt
import pandas as pd

reference_pc = snakemake.input["reference_pc"]
reference_pc = pd.read_csv(reference_pc, sep="\t", header=None)
reference_pc.columns = ["SAMPLE", "PC1", "PC2", "PC3", "PC4"]

fig, ax = plt.subplots()
ax.scatter(reference_pc["PC1"], reference_pc["PC2"])
plt.savefig(snakemake.output["pc_1_2_plot"])
