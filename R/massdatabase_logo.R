#' @title Show the base information of massdatabase pacakge
#' @description Show the base information of massdatabase pacakge.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @return A ASCII log of massdatabase
#' @importFrom magrittr %>%
#' @importFrom dplyr case_when everything select filter
#' @importFrom purrr map map2 walk
#' @importFrom crayon green
#' @importFrom utils download.file head read.table
#' @export
#' @examples
#' massdatabase_logo()

massdatabase_logo <- function() {
  cat(crayon::green("Thank you for using massdatabase!\n"))
  message(crayon::green("Version", massdatabase_version, "(", update_date, ')\n'))
  cat(crayon::green("More information: google tidymass massdatabase.\n"))
  cat(crayon::green(
    c(
      "                _    _____  ___ ",
      " _ __ ___   ___| |_  \\_   \\/   \\",
      "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /",
      "| | | | | |  __/ |_/\\/ /_/ /_// ",
      "|_| |_| |_|\\___|\\__\\____/___,'  ",
      "                                "
    )

  ), sep = "\n")
}

massdatabase_version <- packageVersion(pkg = "massdatabase")
update_date <- as.character(Sys.time())


# library(cowsay)
# # https://onlineasciitools.com/convert-text-to-ascii-art
# # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("logo.txt")
# dput(art)
# massdatabase_logo <-
#   c("                _    _____  ___ ", " _ __ ___   ___| |_  \\_   \\/   \\",
#     "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /", "| | | | | |  __/ |_/\\/ /_/ /_// ",
#     "|_| |_| |_|\\___|\\__\\____/___,'  ", "                                "
#   )
# cat(massdatabase_logo, sep = "\n")
