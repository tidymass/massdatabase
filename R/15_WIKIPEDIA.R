#' @title Request the scientific classification of one species from wikipedia
#' @description Request the scientific classification of one species from wikipedia
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param species_id species_id
#' @return A data frame.
#' @importFrom rvest read_html html_table html_element
#' @importFrom dplyr filter
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_wikipedia_scientific_classification(species_id = "Aaptos ciliata")
#' x

request_wikipedia_scientific_classification <-
  function(species_id = "Aaptos ciliata") {
    url <-
      paste0("https://en.wikipedia.org/wiki/",
             stringr::str_replace_all(species_id, " ", "_"))

    result <-
      tryCatch(
        rvest::read_html(url),
        error = function(e)
          NULL
      )

    if (is.null(result)) {
      new_species_id <-
        stringr::str_replace_all(species_id, " ", "_")
      if (length(stringr::str_extract_all(new_species_id, "_")[[1]]) == 2) {
        new_species_id <-
          head(stringr::str_split(new_species_id, "_")[[1]], 2) %>%
          paste0(collapse = "_")
        result <-
          tryCatch(
            rvest::read_html(
              paste0("https://en.wikipedia.org/wiki/",
                     new_species_id)
            ),
            error = function(e)
              NULL
          )

        if(is.null(result)){
          return(NA)
        }

      }else{
        return(NA)
      }
    }

    result <-
      tryCatch(result %>%
                 rvest::html_element(".biota") %>%
                 rvest::html_table() %>%
                 as.data.frame(), error = function(e) NULL)

    if(is.null(result)){
      return(NA)
    }

    colnames(result) <-
      c("name", "value")
    result <-
      result %>%
      dplyr::filter(name != value)

    result$name <-
      result$name %>%
      stringr::str_replace("\\:", "")
    result
  }









#' @title Request the information of one compound
#' @description Request the information of one compound
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param compound_id compound_id name
#' @return A data frame.
#' @importFrom rvest read_html html_table html_element
#' @importFrom dplyr filter
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
#' @examples
#' request_wikipedia_compound(compound_id = "Microcystin-LR")
#' request_wikipedia_compound(compound_id = "Glucose")

request_wikipedia_compound <-
  function(compound_id = "Microcystin-LR") {
    url <-
      paste0("https://en.wikipedia.org/wiki/",
             stringr::str_replace_all(compound_id, " ", "_"))

    result <-
      tryCatch(
        rvest::read_html(url),
        error = function(e)
          NULL
      )

    if (is.null(result)) {
      return(NA)
    }

    result <-
      tryCatch(result %>%
                 rvest::html_element(".ib-chembox") %>%
                 rvest::html_table() %>%
                 tibble::as_tibble(),
               error = function(e) NULL)

    if(is.null(result)){
      return(NA)
    }

    result <-
      result %>%
      dplyr::filter(!is.na(X1)) %>%
      dplyr::filter(X1 != "")

    colnames(result) <-
      c("name", "value")

    idx <- which(result$name == "Names")

    result <-
      result[idx:nrow(result),]

    # result$name[grep("IUPAC name", result$name)] <-
    #   result$name[grep("IUPAC name", result$name)] %>%
    #   stringr::str_split(pattern = "\n") %>%
    #   `[[`(1) %>%
    #   `[`(1)
    #
    # result$value[grep("IUPAC name", result$name)] <-
    #   result$value[grep("IUPAC name", result$name)] %>%
    #   stringr::str_split(pattern = "\n") %>%
    #   `[[`(1) %>%
    #   `[`(2)
    #
    # result$name[grep("Other names", result$name)] <-
    #   result$name[grep("Other names", result$name)] %>%
    #   stringr::str_split(pattern = "\n") %>%
    #   `[[`(1) %>%
    #   `[`(1)
    #
    # result$value[grep("Other names", result$name)] <-
    #   result$value[grep("Other names", result$name)] %>%
    #   stringr::str_split(pattern = "\n") %>%
    #   `[[`(1) %>%
    #   `[`(2)
    #
    # result$name[grep("InChI", result$name)] <-
    #   result$name[grep("InChI", result$name)] %>%
    #   stringr::str_split(pattern = "\n") %>%
    #   `[[`(1) %>%
    #   `[`(1)
    #
    # result$value[grep("InChI", result$name)] <-
    #   result$value[grep("InChI", result$name)] %>%
    #   stringr::str_split(pattern = "\n") %>%
    #   `[[`(1) %>%
    #   `[`(2)
#
#     idx <-
#       which(result$name == result$value)
#
#     idx <-
#     data.frame(idx1 = idx,
#                idx2 = c(idx-1, nrow(result))[-1])
#
#     seq_len(nrow(idx)) %>%
#       purrr::map(function(i){
#         result[idx$idx1[i] : idx$idx2[i], ]
#       })
#
#
#     result$name <-
#       result$name %>%
#       stringr::str_split("\n", "")
    result
  }

