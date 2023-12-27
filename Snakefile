configfile: "config/config.yaml"
sample_numbers = config["sample_numbers"]

rule download_human_genome_and_index:
  """
  Download GRCh38 reference genome and the corresponding index
  """
  params:
    human_genome_url = config["human_genome_url"],
    human_genome_index_url = config["human_genome_index_url"],
    
  output:
    genome = 'inputs/genome/GRCh38.fa',
    genome_index = 'inputs/genome/GRCh38.fa.fai',

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
    resources_url = config["verifybamid2_url"],

  output:
    verify_bam_id_resources = expand("inputs/verifybamid_resources/1000g.phase3.100k.b38.vcf.gz.dat.{ext}", ext=["UD", "V", "bed", "mu"]),
  
  shell:
    """
    git clone {params.resources_url}
    mv verifybamid/resources/1000g.phase3.100k.b38.vcf.gz.dat* inputs/verifybamid_resources/
    rm -rf verifybamid
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
    cram_files = expand("inputs/crams/{sample_numbers}.GRCh38.low_coverage.cram", sample_numbers=sample_numbers),
    cram_index_files = expand("inputs/crams/{sample_numbers}.GRCh38.low_coverage.cram.fai", sample_numbers=sample_numbers),
    thousandG_reference_populations = "inputs/crams/1000G_reference_populations.txt",

  shell:
    """
    git clone {params.inputs_url} {params.temp_dir}
    mv {params.temp_dir}/input/* inputs/crams/
    rm -rf skills_test
    """

rule verifybamid:
  """
  Run the verifybamid2 tool on the downloaded files using wildcards to extend its functionality
  """
  input:
    reference = "inputs/genome/GRCh38.fa",
    bam_file = expand("inputs/crams/{sample_numbers}.GRCh38.low_coverage.cram", sample_numbers=sample_numbers),
     
  output:
    ancestry = expand("results/verifybamid/{sample_numbers}.Ancestry", sample_numbers=sample_numbers),
    selfSM = expand("results/verifybamid/{sample_numbers}.selfSM", sample_numbers=sample_numbers),

  params:
    svd_prefix = "inputs/verifybamid_resources/1000g.phase3.100k.b38.vcf.gz.dat",
    num_pc = 4,
    output_prefix = lambda wildcards: "results/verifybamid/{wildcards.sample_numbers}",

  shell:
    """
    verifybamid2 \
      --SVDPrefix {params.svd_prefix} \
      --Reference {input.reference} \
      --BamFile {input.bam_file} \
      --NumPC {params.num_pc} \
      --Output {params.output_prefix}
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
    selfSM = expand("results/verifybamid/{sample_numbers}.selfSM", sample_numbers=sample_numbers),

  output:
    all_samples_contamination = "results/verifybamid/all.selfSM",

  params:
    all_samples_glob = "results/verifybamid/*.selfSM",

  shell:
    """
    cat {params.all_samples_glob} | \
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
    ancestry = expand("results/verifybamid/{sample_numbers}.Ancestry", sample_numbers=sample_numbers),
    reference_pc = "inputs/verifybamid_resources/1000g.phase3.100k.b38.vcf.gz.dat.V",
    thousandG_reference_populations = "inputs/crams/1000G_reference_populations.txt",

  output:
    pc12_plot = "results/PC1_PC2.png",
    pc23_plot = "results/PC2_PC3.png",
    pc34_plot = "results/PC3_PC4.png",
    pc123_plot = "results/PC1_PC2_PC3.png",

  params:
    ancestry_filenames = lambda wildcards: sample_numbers,

  script:
    "scripts/generate_pc_plots.py"
