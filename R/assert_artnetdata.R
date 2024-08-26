assert_artnetdata <- function() {
  if (system.file(package = "ARTnetData") == "") {
    warning(
      "This function requires the `ARTnetData` package to be installed.\n",
      "Follow the instructions at the link below to get access to it.\n",
      "https://github.com/EpiModel/ARTnet/tree/main?tab=readme-ov-file#artnetdata-dependency\n"
    )
    return(FALSE)
  }
  return(TRUE)
}
