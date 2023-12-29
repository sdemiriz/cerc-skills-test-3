import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

# Import reference and study populations, add PC column when necessary
reference_pc = pd.read_csv(snakemake.input["reference_pc"], sep="\t", header=None)
reference_pc.columns = ["SAMPLE", "PC1", "PC2", "PC3", "PC4"]

reference_populations = pd.read_csv(
    snakemake.input["thousandG_reference_populations"], sep="\t", header=None
)
reference_populations.columns = ["SAMPLE", "POPULATION"]

reference_pc_populations = reference_pc.merge(reference_populations, on="SAMPLE")
study_pc_populations = pd.read_csv(snakemake.input["ancestry"], sep="\t")

# Specify colors of each population on the plots
color_map = {
    "EUR": "red",
    "SAS": "blue",
    "EAS": "green",
    "AFR": "magenta",
    "AMR": "turquoise",
}


# Function to create the 2D PC plot due to needing three plots of this kind
def two_axis_pc_plot(first: str, second: str, out_file: str):
    # Define plot
    fig, ax = plt.subplots()

    # Add reference population scatterpoints
    ax.scatter(
        x=reference_pc_populations[first],
        y=reference_pc_populations[second],
        c=reference_pc_populations["POPULATION"].map(color_map),
        alpha=0.3,
    )

    # Add study samples scatterpoints with standout visual
    ax.scatter(
        x=study_pc_populations[first],
        y=study_pc_populations[second],
        c="black",
        marker="^",
    )

    # Annotate study samples with sample numbers
    for i, sample_number in enumerate(study_pc_populations["SAMPLE"]):
        ax.annotate(
            sample_number,
            (study_pc_populations[first][i], study_pc_populations[second][i]),
        )

    # Define custom legend
    custom = [
        Line2D([], [], marker="o", color=clr, linestyle="None")
        for clr in color_map.values()
    ] + [Line2D([], [], marker="^", color="black", linestyle="None")]

    # Add axis labels and custom legend
    plt.xlabel(first)
    plt.ylabel(second)
    plt.legend(handles=custom, labels=[*color_map.keys(), "Study"])

    # Save to file
    plt.savefig(out_file)


# Use defined function
two_axis_pc_plot("PC1", "PC2", snakemake.output["pc12_plot"])
two_axis_pc_plot("PC2", "PC3", snakemake.output["pc23_plot"])
two_axis_pc_plot("PC3", "PC4", snakemake.output["pc34_plot"])

# Define 3D figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# Add reference population scatterpoints
ax.scatter(
    xs=reference_pc_populations["PC1"],
    ys=reference_pc_populations["PC2"],
    zs=reference_pc_populations["PC3"],
    c=reference_pc_populations["POPULATION"].map(color_map),
    alpha=0.3,
)

# Add study samples scatterpoints with standout visual
ax.scatter(
    xs=study_pc_populations["PC1"],
    ys=study_pc_populations["PC2"],
    zs=study_pc_populations["PC3"],
    marker="^",
    c="black",
)

# Annotate study samples with sample numbers
for x, y, z, sample_number in zip(
    study_pc_populations["PC1"],
    study_pc_populations["PC2"],
    study_pc_populations["PC3"],
    study_pc_populations["SAMPLE"],
):
    ax.text(x, y, z, f"{sample_number}", zdir=None)

# Define custom legend
custom = [
    Line2D([], [], marker="o", color=clr, linestyle="None")
    for clr in color_map.values()
] + [Line2D([], [], marker="^", color="black", linestyle="None")]

# Add axis labels and custom legend
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")
plt.legend(
    handles=custom,
    labels=[*color_map.keys(), "Study"],
    bbox_to_anchor=(1.35, 1.15),
)
plt.savefig(snakemake.output["pc123_plot"])
