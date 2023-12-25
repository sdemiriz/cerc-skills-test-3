process helloWorld {
  debug true
  script:
  """
  echo 'Hello world!'
  """
}

workflow {
  helloWorld()
}
