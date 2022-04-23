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
#' x[1:2]
#' y = request_foodb_compound(compound_id = "FDB000004", return_form = "data.frame")
#' head(y)

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
          state = ifelse(is.null(result$state), NA, result$state),
          pathways = ifelse(is.null(result$pathways), NA, result$pathways),
          hmdb_id = ifelse(is.null(result$hmdb_id), NA, result$hmdb_id),
          pubchem_compound_id = ifelse(
            is.null(result$pubchem_compound_id),
            NA,
            result$pubchem_compound_id
          ),
          chemspider_id = ifelse(
            is.null(result$chemspider_id),
            NA,
            result$chemspider_id
          ),
          kegg_id = ifelse(is.null(result$kegg_id), NA, result$kegg_id),
          chebi_id = ifelse(is.null(result$chebi_id), NA, result$chebi_id),
          biocyc_id = ifelse(is.null(result$biocyc_id), NA, result$biocyc_id),
          het_id = ifelse(is.null(result$het_id), NA, result$het_id),
          wikipidia = ifelse(is.null(result$wikipidia), NA, result$wikipidia),
          vmh_id = ifelse(is.null(result$vmh_id), NA, result$vmh_id),
          fbonto_id = ifelse(is.null(result$fbonto_id), NA, result$fbonto_id),
          foodb_id = ifelse(is.null(result$foodb_id), NA, result$foodb_id),
          general_references = ifelse(
            is.null(result$general_references),
            NA,
           result$general_references
          ),
          foods = foods,
          flavors = ifelse(is.null(result$flavors), NA, result$flavors),
          enzymes = ifelse(is.null(result$enzymes), NA, result$enzymes),
          health_effects = ifelse(
            is.null(result$health_effects),
            NA,
            result$health_effects
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



#' @title Request MS2 spectra of one compound in FoodB
#' @description Request MS2 spectra of one compound in FoodB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param compound_id compound id. For example, FDB000004.
#' @return A data frame
#' @importFrom XML xmlTreeParse xmlToList
#' @importFrom magrittr %>%
#' @importFrom rvest read_html html_element html_attr
#' @export
#' @examples
#' x = request_foodb_compound_ms2(compound_id = "FDB000004")

request_foodb_compound_ms2 <-
  function(compound_id = "FDB000013") {
    url <- paste0("https://foodb.ca/compounds/", compound_id)

    result <-
      readLines(paste0(url, ".xml"))

    result <-
      XML::xmlTreeParse(file = result, asText = TRUE)

    result <-
      XML::xmlToList(result)

    spectra <-
      result$spectra %>%
      as.data.frame()

    idx <-
      which(unname(unlist(spectra[1,])) == "Specdb::MsMs")

    if (length(idx) == 0) {
      message('No MS/MS.')
      return(NA)
    }

    ms2_id <-
      unname(unlist(spectra[2, idx, drop = TRUE]))

    ms2_url <- "https://foodb.ca/spectra/ms_ms/"

    ms2_spectra <-
      lapply(ms2_id, function(temp_id) {
        # cat(temp_id, " ")
        html_document <-
          rvest::read_html(paste0(ms2_url, temp_id))
        link <-
          html_document %>%
          rvest::html_element("tr:nth-child(1) a") %>%
          rvest::html_attr("href")
        if (is.na(link)) {
          return(NULL)
        }
        ms2 <-
          read.table(link, header = FALSE)
        colnames(ms2) <-
          c("mz", "intensity")

        ms1_info <-
          html_document %>%
          rvest::html_table()

        ms1_info <-
          rbind(ms1_info[[1]],
                ms1_info[[2]])

        colnames(ms1_info) <-
          c("name", "value")

        list(ms1_info = ms1_info,
             ms2 = ms2)

      })

    names(ms2_spectra) <- ms2_id
    ms2_spectra

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
