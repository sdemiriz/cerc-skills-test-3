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
    git clone https://github.com/Griffan/VerifyBamID.git
    mv VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat* inputs/VerifyBamID_resource/
    rm -rf VerifyBamID
    """

rule download_skills_test_3_inputs:
  params:
    inputs_url = config["skills_test_url"],

  output:
    "inputs/skills_test_3_inputs/HGDP{sample_number}.GRCh38.low_coverage.cram",
    "inputs/skills_test_3_inputs/HGDP{sample_number}.GRCh38.low_coverage.cram.crai",

  shell:
    """
    git clone https://github.com/CERC-Genomic-Medicine/skills_test_3.git
    mv skills_test_3/input/* inputs/skills_test_3_inputs/
    rm -rf skills_test_3
    """
