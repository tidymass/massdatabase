#' @title Get the similarity between two words
#' @description Get the similarity between two words
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param word1 word1.
#' @param word2 word2.
#' @param ... Other arguments for stringdist::stringsim()
#' @return A simility score.
#' @importFrom stringr str_extract_all str_to_lower
#' @importFrom magrittr %>%
#' @importFrom stringdist stringsim
#' @export
#' @examples
#' get_words_similarity(word1 = "(2S)-2-Amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid",
#'word2 = "(2S)-2-Amino-3-(1-methyl-1H-imidazol-4-yl)propanoate")


get_words_similarity <-
  function(word1 = "(2S)-2-Amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid",
           word2 = "(2S)-2-Amino-3-(1-methyl-1H-imidazol-4-yl)propanoate",
           ...) {
    new_word1 <-
      word1 %>%
      lapply(function(x) {
        x %>%
          stringr::str_extract_all("[0-9a-zA-Z]{1,100}") %>%
          `[[`(1) %>%
          paste(collapse = "") %>%
          stringr::str_to_lower()
      }) %>%
      unlist()

    new_word2 <-
      word2 %>%
      lapply(function(x) {
        x %>%
          stringr::str_extract_all("[0-9a-zA-Z]{1,100}") %>%
          `[[`(1) %>%
          paste(collapse = "") %>%
          stringr::str_to_lower()
      }) %>%
      unlist()

    stringdist::stringsim(a = new_word1,
                          b = new_word2,
                          ...)

  }
