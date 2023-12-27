import matplotlib.pyplot as plt
import pandas as pd

reference_pc = pd.read_csv(snakemake.input["reference_pc"], sep="\t", header=None)
reference_pc.columns = ["SAMPLE", "PC1", "PC2", "PC3", "PC4"]


def two_axis_pc_plot(first: str, second: str, out_file: str):
    fig, ax = plt.subplots()
    ax.scatter(reference_pc[first], reference_pc[second])
    plt.savefig(out_file)


two_axis_pc_plot("PC1", "PC2", snakemake.output["pc_1_2_plot"])
two_axis_pc_plot("PC2", "PC3", snakemake.output["pc_2_3_plot"])
two_axis_pc_plot("PC3", "PC4", snakemake.output["pc_3_4_plot"])
