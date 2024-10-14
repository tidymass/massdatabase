#' @title Request one specific reaction information in Reactome
#' @description Request one specific reaction information in Reactome
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param reaction_id reaction id. For example, reaction_id
#' @return Alist.
#' @importFrom curl curl_download
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_reactome_reaction(reaction_id = "R-HSA-8876188")
#' x$reactants
#' x$products

request_reactome_reaction <-
  function(reaction_id = "R-HSA-8876188") {
    url <-
      paste0("https://reactome.org/ContentService/exporter/event/",
             reaction_id,
             ".sbml")

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    curl::curl_download(url = url,
                        destfile = file.path(temp_file, "file.sbml"))

    result <-
      tryCatch(
        parse_reactome_reaction(file_name = file.path(temp_file, "file.sbml")),
        error = function(e) {
          return(NULL)
        }
      )

    invisible(result)
  }




#' @title Parse the sbml reaction data from Reactome
#' @description Parse the sbml reaction data from Reactome
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file_name file name of the data.
#' @return A data frame or list.
#' @importFrom dplyr filter
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
parse_reactome_reaction <-
  function(file_name) {
    result <-
      readLines(file_name)

    result <-
      XML::xmlToList(result)

    species <-
      result$model$listOfSpecies %>%
      lapply(function(x) {
        x$.attrs
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame() %>%
      dplyr::select(
        -c(
          boundaryCondition,
          compartment,
          constant,
          hasOnlySubstanceUnits,
          metaid,
          sboTerm
        )
      )

    reactants <-
      result$model$listOfReactions$reaction$listOfReactants %>%
      dplyr::bind_rows() %>%
      as.data.frame() %>%
      dplyr::select(-c(constant, id, sboTerm, stoichiometry))

    products <-
      result$model$listOfReactions$reaction$listOfProducts %>%
      dplyr::bind_rows() %>%
      as.data.frame() %>%
      dplyr::select(-c(constant, id, sboTerm, stoichiometry))

    reactants <-
      reactants %>%
      dplyr::left_join(species,
                       by = c("species" = "id"))

    products <-
      products %>%
      dplyr::left_join(species,
                       by = c("species" = "id"))

    result <- list(reactants = reactants,
                   products = products)
  }
