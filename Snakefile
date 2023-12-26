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
    "results/verifybamid/{sample_number}.Ancestry",
    "results/verifybamid/{sample_number}.selfSM",
    "results/flags/verifybamid_ran_flag"

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
  """
  input:
    verifybamid_ran_flag = "results/flags/verifybamid_ran_flag",
  output:
    all_samples_contamination = "results/verifybamid/all.selfSM",
  shell:
  """
  cat results/verifybamid/*.selfSM | sed -e '1p' -e '/#SEQ_ID/d' > {output.all_samples_contamination}
  """
