
#' Epidemic and Network Model Parameterization with the ARTnet Dataset
#'
#' \tabular{ll}{
#'    Package: \tab ARTnet\cr
#'    Type: \tab Package\cr
#'    Version: \tab 1.0.0\cr
#'    Date: \tab 2019-05-23\cr
#'    License: \tab GPL-3\cr
#'    LazyLoad: \tab yes\cr
#' }
#'
#' @name ARTnet-package
#' @aliases ARTnet
#' @import ARTnetData EpiModel dplyr
#' @importFrom stats binomial glm lm median poisson predict rbinom runif
#' @importFrom utils head
#'
#' @docType package
#' @keywords package
#'
#' @details
#' The ARTnet package provides a suite of functions for the parameterization
#' of epidemic and network models using data from the ARTnet cross-sectional survery of men
#' who have sex with men (MSM) in the United States. Three functions make up this workflow
#'  and are described below:
#'
#' @section EpiStats Is This Title Long Enough Is That The Reason:
#' \code{build_epistats} builds epidemic models governing act rates and probability of condom use
#' among main, casual and one-of sexual partnerships. It takes as parameters the following variables:
#' \itemize{
#' \item{\code{geog.lvl}}: level of geographic stratification desired. Acceptable values are "city",
#'  "state", "region", "division", and "all" corresponding to cirty, state, census region, census
#'  division and complete geographic area respectively.
#' \item{\code{geog.cat}}: given a geographic level above, `geog.cat` gives the desired feature of
#' interest. Acceptable values are based on the chosen geographic level:
#'   \itemize{}
#' \item{\code{city}}: Atlanta, Boston, Chicago, Dallas, Denver, Detroit, Houston, Los Angeles, Miami,
#' New York City, Philadelphia, San Diego, San Franciso, Seattle, Washington DC
#' \item{\code{state}}: AK, AL, AR, AZ, CA, CO, CT, DC, DE, FL, GA, HI, IA, ID, IL, IN, KS, KY, LA, MA,
#'  MD, ME, MI, MN, MO, MS, MT, NC, ND, NE, NH, NJ, NM, NV, NY, OH, OK, OR, PA, RI, SC, SD, TN, TX, UT,
#'   VA, VT, WA, WI, WV, WY
#' \item{\code{division}}: 1 (New England), 2 (Middle Atlantic), 3 (East North Central),
#' 4 (West North Central), 5 (South Atlantic), 6 (East South Central), 7 (West South Central), 8 (Mountain) 9 (Pacific)
#' \item{\code{region}}: 1 (Northeast), 2 (Midwest), 3 (South), 4 (North)
#' \item{\code{all}}:  No input required if `geog.lvl` is set to "all".
#'  }
#' \item{\code{race}}: whether to introduce modeling by racial stratification. TRUE or FALSE.
#' FALSE by default.
#' \item{\code{age.lim}}: a vector giving the lower and upper limit for the age of interest. Set
#' to `c(15, 65)` by default.
#' \item{\code{age.bks}}: a vector giving the upper age breaks to categorize data by age. Must be
#' within the bounds specified by `age.lim`
#' }
#'
#' @section NetParams
#' builds network formation models. It takes as parameters the folliwing:
#'
#' \itemize
#' \item{\code{epistats}}: an object created by `build_epistats`.
#' \item{\code{smooth.main.dur.old}}: whether to smooth partnership durations of the oldest age group
#' by averaging over current oldest and second oldest age groups.
#' \end{itemize}
#'
#' @section NetStats
#' builds network attributes as well as target statsitics for said network. It takes as parameters the following:
#'
#' \itemize
#' \item{\code{epistats}}: an object created by \code{build_epistats}.
#' \item{\code{netparams}}: an objects created by \code{build_netparams}.
#' \item{\code{expect.mort}}: the expected mortality rate of population of interest. Numeric valued and > 0.
#' \item{\code{network.size}}: desired network size. For two group networks, this is the size of each
#' group. Integer valued and > 0.
#' \end{itemize}
NULL
