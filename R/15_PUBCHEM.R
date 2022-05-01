#'
#'
#' #' @title Request the information of one compound
#' #' @description Request the information of one compound
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param compound_id compound_id name
#' #' @return A data frame.
#' #' @importFrom rvest read_html html_table html_element
#' #' @importFrom dplyr filter
#' #' @importFrom stringr str_replace_all
#' #' @importFrom magrittr %>%
#' #' @export
#' #' @examples
#' #' request_wikipedia_compound(compound_id = "Microcystin-LR")
#' #' request_wikipedia_compound(compound_id = "Glucose")
#'
#' request_wikipedia_compound <-
#'   function(compound_id = "Microcystin-LR") {
#'     url <-
#'       paste0("https://en.wikipedia.org/wiki/",
#'              stringr::str_replace_all(compound_id, " ", "_"))
#'
#'
#'     result <-
#'       read_html("https://pubchem.ncbi.nlm.nih.gov/compound/9307")
#'
#'
#'     result %>%
#'       html_elements("p")
#'
#'
#'     result <-
#'       tryCatch(
#'         rvest::read_html(url),
#'         error = function(e)
#'           NULL
#'       )
#'
#'     if (is.null(result)) {
#'       return(NA)
#'     }
#'
#'     result <-
#'       tryCatch(result %>%
#'                  rvest::html_element(".ib-chembox") %>%
#'                  rvest::html_table() %>%
#'                  tibble::as_tibble(),
#'                error = function(e) NULL)
#'
#'     if(is.null(result)){
#'       return(NA)
#'     }
#'
#'     result <-
#'       result %>%
#'       dplyr::filter(!is.na(X1)) %>%
#'       dplyr::filter(X1 != "")
#'
#'     colnames(result) <-
#'       c("name", "value")
#'
#'     idx <- which(result$name == "Names")
#'
#'     result <-
#'       result[idx:nrow(result),]
#'
#'     # result$name[grep("IUPAC name", result$name)] <-
#'     #   result$name[grep("IUPAC name", result$name)] %>%
#'     #   stringr::str_split(pattern = "\n") %>%
#'     #   `[[`(1) %>%
#'     #   `[`(1)
#'     #
#'     # result$value[grep("IUPAC name", result$name)] <-
#'     #   result$value[grep("IUPAC name", result$name)] %>%
#'     #   stringr::str_split(pattern = "\n") %>%
#'     #   `[[`(1) %>%
#'     #   `[`(2)
#'     #
#'     # result$name[grep("Other names", result$name)] <-
#'     #   result$name[grep("Other names", result$name)] %>%
#'     #   stringr::str_split(pattern = "\n") %>%
#'     #   `[[`(1) %>%
#'     #   `[`(1)
#'     #
#'     # result$value[grep("Other names", result$name)] <-
#'     #   result$value[grep("Other names", result$name)] %>%
#'     #   stringr::str_split(pattern = "\n") %>%
#'     #   `[[`(1) %>%
#'     #   `[`(2)
#'     #
#'     # result$name[grep("InChI", result$name)] <-
#'     #   result$name[grep("InChI", result$name)] %>%
#'     #   stringr::str_split(pattern = "\n") %>%
#'     #   `[[`(1) %>%
#'     #   `[`(1)
#'     #
#'     # result$value[grep("InChI", result$name)] <-
#'     #   result$value[grep("InChI", result$name)] %>%
#'     #   stringr::str_split(pattern = "\n") %>%
#'     #   `[[`(1) %>%
#'     #   `[`(2)
#' #
#' #     idx <-
#' #       which(result$name == result$value)
#' #
#' #     idx <-
#' #     data.frame(idx1 = idx,
#' #                idx2 = c(idx-1, nrow(result))[-1])
#' #
#' #     seq_len(nrow(idx)) %>%
#' #       purrr::map(function(i){
#' #         result[idx$idx1[i] : idx$idx2[i], ]
#' #       })
#' #
#' #
#' #     result$name <-
#' #       result$name %>%
#' #       stringr::str_split("\n", "")
#'     result
#'   }
#'
