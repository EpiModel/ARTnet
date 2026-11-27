out <- readRDS("targets.rds"); co <- readRDS("coefs.rds")
tab <- do.call(rbind, lapply(out, `[[`, "rows"))
vars <- c("main","casl","inst","conc_main","conc_casl","dur_main","dur_casl","dur_main_max","nmr_main_share","nmr_casl_share","nmr_inst_share",
          "nmr_main_B","nmr_main_H","nmr_main_W","nmr_casl_B","nmr_casl_H","nmr_casl_W","nmr_inst_B","nmr_inst_H","nmr_inst_W")
cur <- tab[tab$scen == "current", ]; sh <- tab[tab$scen == "nested_dd", ]; fl <- tab[tab$scen == "flat_dd", ]
options(width = 250)
for (v in vars) {
  cat(sprintf("%-15s", v)); for (i in seq_len(nrow(cur))) cat(sprintf(" %s %8.1f -> %8.1f (%+5.0f%%, flat %+5.0f%%)", cur$city[i], cur[[v]][i] * ifelse(grepl("share", v), 100, 1),
      sh[[v]][i] * ifelse(grepl("share", v), 100, 1), 100 * (sh[[v]][i] / cur[[v]][i] - 1), 100 * (fl[[v]][i] / cur[[v]][i] - 1))); cat("\n")
}
# CIs on current city shifts (rate/odds/duration ratios vs rest of ARTnet)
k <- co[co$city %in% c("det","mia") & co$model %in% c("joint_model","joint_concurrent_model","joint_duration_model","joint_nm_race_model"), ]
se <- ifelse(is.na(k$se_quasi), k$se_cur, k$se_quasi)
k$ratio <- exp(k$b_cur); k$lo <- exp(k$b_cur - 1.96 * se); k$hi <- exp(k$b_cur + 1.96 * se); k$shrunk <- exp(k$b_nest)
k$nb <- exp(k$b_nbfe); k$nb_lo <- exp(k$b_nbfe - 1.96 * k$se_nbfe); k$nb_hi <- exp(k$b_nbfe + 1.96 * k$se_nbfe)
print(cbind(k[, c("city","layer","model","n_city")], round(k[, c("ratio","lo","hi","nb","nb_lo","nb_hi","shrunk")], 2)), row.names = FALSE)
