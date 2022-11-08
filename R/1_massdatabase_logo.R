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
#' @importFrom ggplot2 aes ggplot
#' @importFrom Rdisop getMass getMolecule
#' @importFrom masstools name_duplicated convert_precursor_mz2accurate_mass
#' @importFrom masstools show_progresser sum_formula
#' @export
#' @examples
#' massdatabase_logo()

massdatabase_logo <-
  function() {
    message("Thank you for using massdatabase!")
    message("Version ", massdatabase_version, " (", update_date, ')')
    message("More information: massdatabase.tidymass.org")
    cat(
      c(
        "                          _____        _        _",
        "                         |  __ \\      | |      | |",
        "  _ __ ___   __ _ ___ ___| |  | | __ _| |_ __ _| |__   __ _ ___  ___",
        " | '_ ` _ \\ / _` / __/ __| |  | |/ _` | __/ _` | '_ \\ / _` / __|/ _ \\",
        " | | | | | | (_| \\__ \\__ \\ |__| | (_| | || (_| | |_) | (_| \\__ \\  __/",
        " |_| |_| |_|\\__,_|___/___/_____/ \\__,_|\\__\\__,_|_.__/ \\__,_|___/\\___|",
        ""
      ),
      sep = "\n"
    )
  }

massdatabase_version <-
  as.character(packageVersion(pkg = "massdatabase"))
update_date <- as.character(Sys.time())


# library(cowsay)
# # https://onlineasciitools.com/convert-text-to-ascii-art
# # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("logo.txt")
# dput(art)
# massdatabase_logo <-
#   c("                          _____        _        _", "                         |  __ \\      | |      | |",
#     "  _ __ ___   __ _ ___ ___| |  | | __ _| |_ __ _| |__   __ _ ___  ___",
#     " | '_ ` _ \\ / _` / __/ __| |  | |/ _` | __/ _` | '_ \\ / _` / __|/ _ \\",
#     " | | | | | | (_| \\__ \\__ \\ |__| | (_| | || (_| | |_) | (_| \\__ \\  __/",
#     " |_| |_| |_|\\__,_|___/___/_____/ \\__,_|\\__\\__,_|_.__/ \\__,_|___/\\___|",
#     "")
# cat(massdatabase_logo, sep = "\n")
#
#
#
