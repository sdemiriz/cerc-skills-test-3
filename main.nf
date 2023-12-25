process downloadGenomeAndIndex {
  input:
    val human_genome_url
    val human_genome_index_url

  output:
    path 'GRCh38.fa'
    path 'GRCh38.fa.fai'

  script:
    """
    curl $human_genome_url -o GRCh38.fa
    curl $human_genome_index_url -o GRCh38.fa.fai
    """
}

process downloadInputs {
  debug true
  input:
    val skills_test_url

  output:
    path 'input/'

  script:
    """
    git clone $skills_test_url skills_test_3
    mv skills_test_3/input .
    """
}

process downloadVerifyBamIdResource {
  debug true
  input:
    val verify_bam_id_url

  output:
    path '1000g.phase3.100k.b38.vcf.gz.dat*'

  script:
    """
    git clone $verify_bam_id_url verify_bam_id
    mv verify_bam_id/resource/1000g.phase3.100k.b38.vcf.gz.dat* .
    """
}

workflow {
  downloadGenomeAndIndex(params.human_genome_url, params.human_genome_index_url)
  downloadInputs(params.skills_test_url)
  downloadVerifyBamIdResource(params.verify_bam_id_url)
}
