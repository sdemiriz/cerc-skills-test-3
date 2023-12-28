from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
import pandas as pd

# Import and set column names for reference samples, study samples and 1000G population files
ancestry = pd.read_csv(snakemake.input["ancestry"], sep="\t")
reference_pc = pd.read_csv(snakemake.input["reference_pc"], sep="\t", header=None)
reference_pc.columns = ancestry.columns

thousandG_reference_populations = pd.read_csv(
    snakemake.input["thousandG_reference_populations"], sep="\t", header=None
)
thousandG_reference_populations.columns = ["SAMPLE", "POPULATION"]

# Join 1000G populations with reference panels
reference_pc_populations = reference_pc.merge(thousandG_reference_populations)

# Set up Nearest Neighbors Classifier
kn_classifier = KNeighborsClassifier(
    n_neighbors=snakemake.params["n_neighbors"],
)

# Set up train-test split for validation purposes
X = reference_pc_populations[["PC1", "PC2", "PC3", "PC4"]]
y = reference_pc_populations["POPULATION"]
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=42)

# Train classifier on training set
kn_classifier.fit(X_train, y_train)

# Good results: 0.9952076677316294 on my machine for n=5 neighbors
print(kn_classifier.score(X_test, y_test))

# Classify study samples using the trained model
ancestry["POPULATION"] = kn_classifier.predict(ancestry[["PC1", "PC2", "PC3", "PC4"]])

# Write classified study samples with predicted population labels to file
ancestry.to_csv(snakemake.output["classified_samples"], sep="\t", index=None)
