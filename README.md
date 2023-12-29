# CERC Skills Test 3
## by: Sedat Demiriz

Find my submission, a Snakemake workflow, to Skills Test 3 for my application to CERC, for Fall 2024 intake.

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

After the downloads complete, execute the workflow by running while in the same directory, where
N can be any number between 1-10 inclusive, due to the workflow running on at most 10 samples at
the same time:
```
snakemake -cN
```

After workflow completes, the desired tab-separated files and plots should become available in the 
`deliverables/` directory.
