# CERC Skills Test 3
by: Sedat Demiriz

## Introduction
Find my submission for CERC Skills Test 3, a Snakemake workflow. The workflow is written in Snakemake 
as the workflow manager, uses Bash and Python for data transformation tasks and plotting.

## Running the workflow
To run the workflow, please have `conda` installed on your system. After downloading and unpacking the contents 
of the submission, `cd` into the unpacked directory and run the following command to create 
the `conda` environment to be able to run the workflow.
```
conda env create -n st3 -f envs/environment.yml
```
* All necessary tools, packages etc. for running the tools as well as the `conda` channels they can be found
in can be found in the `envs/environment.yml` file. I did not include them here because `conda` should be 
able to automatically fetch and install everything necessary for this workflow to run.

Activate the new `st3` environment by running after all packages and dependencies have been fetched and 
installed:
```
conda activate st3
```

Finally, execute the workflow by running the following command while in the same directory,
specifying the number of cores to use, `N`, as any number between 1-10 inclusive, based on your hardware. 
As the workflow runs on a maximum of 10 samples at a time, any greater number will not significantly 
improve runtime running on fewer than 10 cores will increase overall runtime due to batching of jobs. 
```
snakemake -cN
```

## Walk-through
During execution, a number of downloads will first occur. These include CRAM files from the Skills Test 3 
repository, reference panel files from the `verifybamid2` tool's repository as well as the human genome
and its index from the URLs mentioned in the Skills Test repository README.

Next, following the instructions for Skills Test 3, `verifybamid2` is run on the downloaded CRAM and reference
panel files, resulting in contamination `.selfSM` and ancestry `.Ancestry` files being generated epr sample
(total of 10 files each).

Next, the results from both sets of 10 files are concatenated using Bash or Python and transformed
according Skills Test 3 instructions.

PCA plots are generated using `matplotlib` in Python and populations are assigned to the 10 samples using a 
`sklearn` K Nearest Neighbors Classifier model, trained on the reference samples and their PC coordinates 
and population labels.

## Results
After the workflow completes, copies of the requested text files and plots should become available in the 
`deliverables/` directory for convenience.

