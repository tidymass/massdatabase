#' @title Request all the compound information in FoodB
#' @description Request all the compound information in FoodB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://foodb.ca/compounds".
#' @param sleep Default is 1 second.
#' @param pages default is from 1:2838
#' @return A data frame.
#' @importFrom rvest read_html html_table
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_foodb_compound_info(pages = 1)
#' head(x)

request_foodb_compound_info <-
  function(url = "https://foodb.ca/compounds",
           sleep = 1,
           pages = c(1:2838)) {
    result <-
      purrr::map(pages, function(idx) {
        cat(idx, " ")
        Sys.sleep(sleep)
        new_url <-
          paste0(url, "?page=", idx)
        x <-
          rvest::read_html(x = new_url)
        
        x <-
          x %>%
          rvest::html_table(fill = TRUE) %>%
          `[[`(1)
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    invisible(result)
  }



#' @title Request one specific the compound information in FoodB
#' @description Request one specific the compound information in FoodB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://foodb.ca/compounds".
#' @param compound_id compound id. For example, FDB000004.
#' @param return_form data.frame or list.
#' @return A data frame or list.
#' @importFrom XML xmlTreeParse xmlToList
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_foodb_compound( compound_id = "FDB000004", return_form = "list")
#' x
#' y = request_foodb_compound( compound_id = "FDB000004", return_form = "data.frame")
#'

request_foodb_compound <-
  function(url = "https://foodb.ca/compounds",
           compound_id = "FDB000004",
           return_form = c("list" , "data.frame")) {
    return_form <- match.arg(return_form)
    result <-
      readLines(paste0(url, "/", compound_id, ".xml"))
    result <-
      XML::xmlTreeParse(file = result, asText = TRUE)
    result <-
      XML::xmlToList(result)
    
    # names(result)
    
    if (return_form == "list") {
      result <- result
    } else{
      version <- result$version
      foods <- paste(unlist(result$foods[1,]), collapse = "{}")
      result <-
        data.frame(
          version = version,
          creation_date = result$creation_date,
          update_date = result$update_date,
          accession = result$accession,
          name = result$name,
          description = result$description,
          synonyms = paste(unlist(result$synonyms), collapse = "{}"),
          chemical_formula = result$chemical_formula,
          average_molecular_weight = as.numeric(result$average_molecular_weight),
          monisotopic_moleculate_weight = as.numeric(result$monisotopic_moleculate_weight),
          iupac_name = ifelse(is.null(result$iupac_name), NA, result$iupac_name),
          traditional_iupac = ifelse(
            is.null(result$traditional_iupac),
            NA,
            result$traditional_iupac
          ),
          cas_registry_number = ifelse(
            is.null(result$cas_registry_number),
            NA,
            result$cas_registry_number
          ),
          smiles = ifelse(is.null(result$smiles), NA, result$smiles),
          inchi = ifelse(is.null(result$inchi), NA, result$inchi),
          inchikey = ifelse(is.null(result$inchikey), NA, result$inchikey),
          state = ifelse(is.null(result$state), NA, is.null(result$state)),
          pathways = ifelse(is.null(result$pathways), NA, is.null(result$pathways)),
          hmdb_id = ifelse(is.null(result$hmdb_id), NA, is.null(result$hmdb_id)),
          pubchem_compound_id = ifelse(
            is.null(result$pubchem_compound_id),
            NA,
            is.null(result$pubchem_compound_id)
          ),
          chemspider_id = ifelse(
            is.null(result$chemspider_id),
            NA,
            is.null(result$chemspider_id)
          ),
          kegg_id = ifelse(is.null(result$kegg_id), NA, is.null(result$kegg_id)),
          chebi_id = ifelse(is.null(result$chebi_id), NA, is.null(result$chebi_id)),
          biocyc_id = ifelse(is.null(result$biocyc_id), NA, is.null(result$biocyc_id)),
          het_id = ifelse(is.null(result$het_id), NA, is.null(result$het_id)),
          wikipidia = ifelse(is.null(result$wikipidia), NA, is.null(result$wikipidia)),
          vmh_id = ifelse(is.null(result$vmh_id), NA, is.null(result$vmh_id)),
          fbonto_id = ifelse(is.null(result$fbonto_id), NA, is.null(result$fbonto_id)),
          foodb_id = ifelse(is.null(result$foodb_id), NA, is.null(result$foodb_id)),
          general_references = ifelse(
            is.null(result$general_references),
            NA,
            is.null(result$general_references)
          ),
          foods = foods,
          flavors = ifelse(is.null(result$flavors), NA, is.null(result$flavors)),
          enzymes = ifelse(is.null(result$enzymes), NA, is.null(result$enzymes)),
          health_effects = ifelse(
            is.null(result$health_effects),
            NA,
            is.null(result$health_effects)
          )
        )
    }
    invisible(result)
  }



#' @title Request all the food information in FoodB
#' @description Request all the food information in FoodB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://foodb.ca/foods".
#' @param sleep Default is 1 second.
#' @param pages default is from 1:2838
#' @return A data frame.
#' @importFrom rvest read_html html_table
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_foodb_food_info(pages = 1)
#' head(x)

request_foodb_food_info <-
  function(url = "https://foodb.ca/foods",
           sleep = 1,
           pages = c(1:32)) {
    result <-
      purrr::map(pages, function(idx) {
        cat(idx, " ")
        Sys.sleep(sleep)
        new_url <-
          paste0(url, "?page=", idx)
        x <-
          rvest::read_html(x = new_url)
        
        x <-
          x %>%
          rvest::html_table(fill = TRUE) %>%
          `[[`(1)
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    invisible(result)
  }



#' #' @title Request one specific the food information in FoodB
#' #' @description Request one specific the food information in FoodB
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param url Default is "https://foodb.ca/foods".
#' #' @param food_id food id. For example, FOOD00971
#' #' @param return_form data.frame or list.
#' #' @return A data frame or list.
#' #' @importFrom XML xmlTreeParse xmlToList
#' #' @importFrom magrittr %>%
#' #' @export
#' #' @examples
#' #' x = request_foodb_food( food_id = "FOOD00971", return_form = "list")
#' #' x
#' #' y = request_foodb_food( food_id = "FOOD00971", return_form = "data.frame")
#' #'
#'
#' request_foodb_food <-
#'   function(url = "https://foodb.ca/foods",
#'            food_id = "FOOD00971",
#'            return_form = c("list" ,"data.frame")) {
#'     return_form <- match.arg(return_form)
#'     result <-
#'       readLines(paste0(url, "/", food_id, ".xml"))
#'     result <-
#'       XML::xmlTreeParse(file = result, asText = TRUE)
#'     result <-
#'       XML::xmlToList(result)
#'
#'     # names(result)
#'
#'     if (return_form == "list") {
#'       result <- result
#'     } else{
#'       version <- result$version
#'       foods <- paste(unlist(result$foods[1, ]), collapse = "{}")
#'       result <-
#'         data.frame(
#'           version = version,
#'           creation_date = result$creation_date,
#'           update_date = result$update_date,
#'           accession = result$accession,
#'           name = result$name,
#'           description = result$description,
#'           synonyms = paste(unlist(result$synonyms), collapse = "{}"),
#'           chemical_formula = result$chemical_formula,
#'           average_molecular_weight = as.numeric(result$average_molecular_weight),
#'           monisotopic_moleculate_weight = as.numeric(result$monisotopic_moleculate_weight),
#'           iupac_name = ifelse(is.null(result$iupac_name), NA, result$iupac_name),
#'           traditional_iupac = ifelse(
#'             is.null(result$traditional_iupac),
#'             NA,
#'             result$traditional_iupac
#'           ),
#'           cas_registry_number = ifelse(
#'             is.null(result$cas_registry_number),
#'             NA,
#'             result$cas_registry_number
#'           ),
#'           smiles = ifelse(is.null(result$smiles), NA, result$smiles),
#'           inchi = ifelse(is.null(result$inchi), NA, result$inchi),
#'           inchikey = ifelse(is.null(result$inchikey), NA, result$inchikey),
#'           state = ifelse(is.null(result$state), NA, is.null(result$state)),
#'           pathways = ifelse(is.null(result$pathways), NA, is.null(result$pathways)),
#'           hmdb_id = ifelse(is.null(result$hmdb_id), NA, is.null(result$hmdb_id)),
#'           pubchem_food_id = ifelse(
#'             is.null(result$pubchem_food_id),
#'             NA,
#'             is.null(result$pubchem_food_id)
#'           ),
#'           chemspider_id = ifelse(
#'             is.null(result$chemspider_id),
#'             NA,
#'             is.null(result$chemspider_id)
#'           ),
#'           kegg_id = ifelse(is.null(result$kegg_id), NA, is.null(result$kegg_id)),
#'           chebi_id = ifelse(is.null(result$chebi_id), NA, is.null(result$chebi_id)),
#'           biocyc_id = ifelse(is.null(result$biocyc_id), NA, is.null(result$biocyc_id)),
#'           het_id = ifelse(is.null(result$het_id), NA, is.null(result$het_id)),
#'           wikipidia = ifelse(is.null(result$wikipidia), NA, is.null(result$wikipidia)),
#'           vmh_id = ifelse(is.null(result$vmh_id), NA, is.null(result$vmh_id)),
#'           fbonto_id = ifelse(is.null(result$fbonto_id), NA, is.null(result$fbonto_id)),
#'           foodb_id = ifelse(is.null(result$foodb_id), NA, is.null(result$foodb_id)),
#'           general_references = ifelse(
#'             is.null(result$general_references),
#'             NA,
#'             is.null(result$general_references)
#'           ),
#'           foods = foods,
#'           flavors = ifelse(is.null(result$flavors), NA, is.null(result$flavors)),
#'           enzymes = ifelse(is.null(result$enzymes), NA, is.null(result$enzymes)),
#'           health_effects = ifelse(
#'             is.null(result$health_effects),
#'             NA,
#'             is.null(result$health_effects)
#'           )
#'         )
#'     }
#'     invisible(result)
#'   }
