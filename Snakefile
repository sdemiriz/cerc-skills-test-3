configfile: "config/config.yaml"

rule all:
  input: 
    "results/HGDP00082.Ancestry"

rule download_human_genome_and_index:
  """
  Download GRCh38 reference genome and the corresponding index
  """
  params:
    human_genome_url = config["human_genome_url"],
    human_genome_index_url = config["human_genome_index_url"],
    
  output:
    genome = 'inputs/GRCh38.fa',
    genome_index = 'inputs/GRCh38.fa.fai',

  shell:
    """
    curl {params.human_genome_url} -o {output.genome}
    curl {params.human_genome_index_url} -o {output.genome_index}
    """

rule download_verifybamid_resources:
  """
  The verifybamid2 tool does provide the necessary reference panel when installed, 
  but a copy can also be downloaded fromthe tool's Github repository so we will do 
  that for future reproducibility
  """
  params:
    resources_url = config["verify_bam_id_url"],

  output:
    verify_bam_id_resources = expand("inputs/VerifyBamID_resource/1000g.phase3.100k.b38.vcf.gz.dat.{ext}", ext=["UD", "V", "bed", "mu"]),
  
  shell:
    """
    git clone {params.resources_url}
    mv VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat* inputs/VerifyBamID_resource/
    rm -rf VerifyBamID
    """

rule download_crams:
  """
  Download the CRAM files and their corresponding index files from the skills_test_3 repository using a flag file
  since the pipeline is not aware of what the files are named beforehand and downloading each pair individually
  would be redundant
  """
  params:
    inputs_url = config["skills_test_url"],
    temp_dir = "skills_test",

  output:
    samples_downloaded_flag = "results/flags/samples_downloaded_flag",

  shell:
    """
    git clone {params.inputs_url} {params.temp_dir}
    mv {params.temp_dir}/input/* inputs/crams/
    rm -rf skills_test
    touch {output.samples_downloaded_flag}
    """

# Quick sanity check to list all downloaded file names
import os
sample_names = [file[:9] for file in os.listdir("inputs/crams/") if file.startswith("HGDP")]
print(sample_names)

rule verify_bam_id:
  """
  Run the verifybamid2 tool on the downloaded files using wildcards to extend its functionality
  """
  input:
    reference = "inputs/GRCh38.fa",
    bam_file = "inputs/crams/{sample_number}.GRCh38.low_coverage.cram",
    samples_downloaded_flag = "results/flags/samples_downloaded_flag",
     
  output:
    #"results/verifybamid/{sample_number}.Ancestry",
    #"results/verifybamid/{sample_number}.selfSM",
    "results/flags/verifybamid_ran_flag",

  params:
    svd_prefix = "inputs/VerifyBamID_resource/1000g.phase3.100k.b38.vcf.gz.dat",
    num_pc = 4,

  shell:
    """
    verifybamid2 \
      --SVDPrefix {params.svd_prefix} \
      --Reference {input.reference} \
      --BamFile {input.bam_file} \
      --NumPC {params.num_pc} \
      --Output results/verifybamid/{wildcards.sample_number}
    """

rule collect_contamination:
  """
  Concatanate all .selfSM files into one, leaving only the two required columns:
  1. Concatanate all sample contamination samples into one
  2. Remove all but the first header in the concatanated file
  3. Select the two desired columns
  4. Rename #SEQ_ID field to SAMPLE
  """
  input:
    verifybamid_ran_flag = "results/flags/verifybamid_ran_flag",

  output:
    all_samples_contamination = "results/verifybamid/all.selfSM",

  params:
    all_samples = "results/verifybamid/*.selfSM",

  shell:
    """
    cat {params.all_samples} | \
    sed -e '1p' -e '/#SEQ_ID/d' | \
    cut -f 1,7 | \
    sed -e 's/#SEQ_ID/SAMPLE/' > \
    {output.all_samples_contamination}
    """

rule generate_pc_plots:
  """
  Generate the plots with the requested Principal Component matchups
  """
  input:
    reference_pc = "inputs/VerifyBamID_resource/1000g.phase3.100k.b38.vcf.gz.dat.V",
    samples_populations = "inputs/crams/1000G_reference_populations.txt",

  output:
    pc_1_2_plot = "results/PC1_PC2.png",
    pc_2_3_plot = "results/PC2_PC3.png",
    pc_3_4_plot = "results/PC3_PC4.png",
    pc_1_2_3_plot = "results/PC1_PC2_PC3.png",

  script:
    "scripts/generate_pc_plots.py"
