assert_artnetdata <- function() {
  if (!system.file(package = "ARTnetData")) {
    stop(
      "This function requires the `ARTnetData` package to be installed.",
      "Follow the instructions at the link below to get access to it.",
      "https://github.com/EpiModel/ARTnet/tree/main?tab=readme-ov-file#artnetdata-dependency"
    )
  }
}
