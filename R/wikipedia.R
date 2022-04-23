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
