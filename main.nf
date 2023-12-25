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

workflow {
  downloadGenomeAndIndex(params.human_genome_url, params.human_genome_index_url)
}
