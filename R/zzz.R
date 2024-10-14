.onAttach <- function(...) {
  # needed <- core[!is_attached(core)]
  # if (length(needed) == 0)
  #   return()
  #
  crayon::num_colors(TRUE)
  massdatabase_attach()

  # if (!"package:conflicted" %in% search()) {
  #   x <- massdatabase_conflicts()
  #   msg(massdatabase_conflict_message(x), startup = TRUE)
  # }
  packageStartupMessage(paste0("massdatabase ", massdatabase_version, " (", update_date, ')'))
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}
