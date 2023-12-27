import matplotlib.pyplot as plt
import pandas as pd

reference_pc = pd.read_csv(snakemake.input["reference_pc"], sep="\t", header=None)
reference_pc.columns = ["SAMPLE", "PC1", "PC2", "PC3", "PC4"]

samples_populations = pd.read_csv(
    snakemake.input["samples_populations"], sep="\t", header=None
)
samples_populations.columns = ["SAMPLE", "POPULATION"]

reference_pc_populations = reference_pc.merge(samples_populations, on="SAMPLE")


def two_axis_pc_plot(first: str, second: str, out_file: str):
    fig, ax = plt.subplots()
    ax.scatter(
        x=reference_pc_populations[first],
        y=reference_pc_populations[second],
        c=pd.factorize(reference_pc_populations["POPULATION"])[0],
        alpha=0.3,
        linewidths=0,
    )
    plt.xlabel(first)
    plt.ylabel(second)
    plt.savefig(out_file)


two_axis_pc_plot("PC1", "PC2", snakemake.output["pc_1_2_plot"])
two_axis_pc_plot("PC2", "PC3", snakemake.output["pc_2_3_plot"])
two_axis_pc_plot("PC3", "PC4", snakemake.output["pc_3_4_plot"])
