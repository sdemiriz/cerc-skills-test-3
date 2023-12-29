import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

reference_pc = pd.read_csv(snakemake.input["reference_pc"], sep="\t", header=None)
reference_pc.columns = ["SAMPLE", "PC1", "PC2", "PC3", "PC4"]

reference_populations = pd.read_csv(
    snakemake.input["thousandG_reference_populations"], sep="\t", header=None
)
reference_populations.columns = ["SAMPLE", "POPULATION"]

reference_pc_populations = reference_pc.merge(reference_populations, on="SAMPLE")
study_pc_populations = pd.read_csv(snakemake.input["ancestry"], sep="\t")

color_map = {
    "EUR": "red",
    "SAS": "blue",
    "EAS": "green",
    "AFR": "magenta",
    "AMR": "turquoise",
}


def two_axis_pc_plot(first: str, second: str, out_file: str):
    fig, ax = plt.subplots()
    ax.scatter(
        x=reference_pc_populations[first],
        y=reference_pc_populations[second],
        c=reference_pc_populations["POPULATION"].map(color_map),
        alpha=0.3,
        linewidths=0,
    )
    ax.scatter(
        x=study_pc_populations[first],
        y=study_pc_populations[second],
        c="black",
        marker="^",
    )
    for i, sample_number in enumerate(study_pc_populations["SAMPLE"]):
        ax.annotate(
            sample_number,
            (study_pc_populations[first][i], study_pc_populations[second][i]),
        )
    custom = [
        Line2D([], [], marker="o", color=clr, linestyle="None")
        for clr in color_map.values()
    ] + [Line2D([], [], marker="^", color="black", linestyle="None")]
    plt.legend(handles=custom, labels=[*color_map.keys(), "Study"])
    plt.xlabel(first)
    plt.ylabel(second)
    plt.savefig(out_file)


two_axis_pc_plot("PC1", "PC2", snakemake.output["pc12_plot"])
two_axis_pc_plot("PC2", "PC3", snakemake.output["pc23_plot"])
two_axis_pc_plot("PC3", "PC4", snakemake.output["pc34_plot"])

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.scatter(
    xs=reference_pc_populations["PC1"],
    ys=reference_pc_populations["PC2"],
    zs=reference_pc_populations["PC3"],
    c=pd.factorize(reference_pc_populations["POPULATION"])[0],
)
ax.scatter(
    xs=study_pc_populations["PC1"],
    ys=study_pc_populations["PC2"],
    zs=study_pc_populations["PC3"],
    marker="^",
    c="black",
)
for x, y, z, sample_number in zip(
    study_pc_populations["PC1"],
    study_pc_populations["PC2"],
    study_pc_populations["PC3"],
    study_pc_populations["SAMPLE"],
):
    ax.text(x, y, z, f"{sample_number}", zdir=None)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")
plt.savefig(snakemake.output["pc123_plot"])
