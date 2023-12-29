# CERC Skills Test 3
by: Sedat Demiriz

## Introduction
Find my submission, a Snakemake workflow, to Skills Test 3 for my application to CERC, for Fall 2024 intake.
The workflow is written in Snakemake for workflow management, uses Bash for simple data transformation tasks,
and Python for more involved data transformation tasks and plotting.

## Running the workflow
To run the workflow, please have `conda` installed on your system. After downloading and unpacking the contents 
of the submission, `cd` into the unpacked directory and run the following command to create 
the `conda` environment to be able to run the workflow.
```
conda env create -n st3 -f envs/environment.yml
```

Activate the new `st3` environment by running:
```
conda activate st3
```

After the installations complete, execute the workflow by running the following command while in the same directory,
where N can be any number between 1-10 inclusive, due to the workflow running on at most 10 samples at the same time.
Any larger number will not significantly improve runtime and smaller numbers will work through the 10 samples in
batches, slowing the execution:
```
snakemake -cN
```

## Results
After the workflow completes, copies of the requested text files and plots should become available in the 
`deliverables/` directory for convenience.
