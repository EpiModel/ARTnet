## REFERENCE COPY (2026-04-19) of EpiModelHIV-Template/R/A-networks/initialize.R
## DO NOT EDIT — this exists to pin the downstream consumer contract.
## If the upstream file changes, refresh this copy and update
## `inst/validation/netstats_contract.md`.

## Initialize the ARTnet data objects and the networks to be fitted
##
## This script should not be run directly. But `sourced` by `1-estimation.R`

if (system.file(package = "ARTnetData") == "") {
  message(
    "=================================================================\n",
    "You are currently using the example population provided by ARTnet\n",
    "Install ARTnetData to get all the features.\n",
    "Follow the instructions at the link below to get access to it.\n",
    "https://github.com/EpiModel/ARTnet/tree/main?tab=readme-ov-file#artnetdata-dependency\n",
    "=================================================================\n"
  )

  epistats <- readRDS(system.file("epistats-example.rds", package = "ARTnet"))
  netstats <- readRDS(system.file("netstats-example.rds", package = "ARTnet"))
} else {
  epistats <- build_epistats(
    geog.lvl = "city",
    geog.cat = "Atlanta",
    init.hiv.prev = c(0.33, 0.137, 0.084),
    race = TRUE,
    time.unit = time_unit
  )

  netparams <- build_netparams(
    epistats = epistats,
    smooth.main.dur = TRUE
  )

  netstats <- build_netstats(
    epistats,
    netparams,
    expect.mort = 0.000478213,
    network.size = networks_size
  )
}


nw <- EpiModel::network_initialize(netstats$demog$num)
nw_main <- EpiModel::set_vertex_attribute(
  nw,
  names(netstats$attr),
  netstats$attr
)

nw_casl <- nw_main
nw_inst <- nw_main
