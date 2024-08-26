#' @import utils

utils::globalVariables(c( "age",
                         "race.cat3",
                         "p_race.cat3",
                         "AMIS_ID",
                         "survey_year",
                         "prep",
                         "ptype",
                         "duration",
                         "comb.age",
                         "geogYN",
                         "race.combo",
                         "RAI",
                         "IAI",
                         "hiv.concord.pos",
                         "anal.acts.week",
                         "anal.acts.week.cp",
                         "RECUAI",
                         "INSUAI",
                         "hiv2",
                         "ongoing",
                         "index.age.grp",
                         "geog",
                         "same.age.grp",
                         "race.dist.city",
                         "race.dist.state",
                         "race.dist.census.division",
                         "race.dist.census.region",
                         "survey.year",
                         "ongoing2"))

missing_data_msg <- paste0(
  "This function requires the `ARTnetData` package to be installed.\n",
  "Follow the instructions at the link below to get access to it.\n",
  "https://github.com/EpiModel/ARTnet/tree/main?tab=readme-ov-file#artnetdata-dependency\n"
)
