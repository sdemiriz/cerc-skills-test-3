configfile: "config/config.yaml"

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
  params:
    resources_url = config["verify_bam_id_url"],

  output:
    verify_bam_id_resources = expand("inputs/VerifyBamID_resource/1000g.phase3.100k.b38.vcf.gz.dat.{ext}", ext=["UD", "V", "bed", "mu"])
  
  shell:
    """
    git clone {params.resources_url}
    mv VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat* inputs/VerifyBamID_resource/
    rm -rf VerifyBamID
    """

rule download_crams:
  params:
    inputs_url = config["skills_test_url"],
    temp_dir = "skills_test"

  output:
    "inputs/crams/1000G_reference_populations.txt",

  shell:
    """
    git clone {params.inputs_url} {params.temp_dir}
    mv {params.temp_dir}/input/* inputs/crams/
    rm -rf skills_test
    """
