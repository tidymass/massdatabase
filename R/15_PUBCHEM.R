#' @title Request the information of one compound from PubChem
#' @description Request the information of one compound from PubChem
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param compound_id compound_id name
#' @return A list.
#' @importFrom rvest read_html html_table html_element
#' @importFrom dplyr filter
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x <- request_pubchem_compound(compound_id = "5288826")
#' names(x)
#' x$Computed_Descriptors
#' x$Molecular_Formula
#' x$Other_Identifiers
#' x$Synonyms

request_pubchem_compound <-
  function(compound_id = "5288826") {
    url <-
      paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/",
        compound_id,
        "/JSON/?response_type=save&response_basename=compound_CID_",
        compound_id
      )

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    curl::curl_download(url = url,
                        destfile = file.path(temp_file, "file.json"))

    result <-
      tryCatch(
        parse_pubchem_compound(file_name = file.path(temp_file, "file.json")),
        error = function(e) {
          return(NULL)
        }
      )

    if (is.null(result)) {
      message('Please check your PubChemID: ', compound_id, " online.")
      return(NULL)
    }

    unlink(x = temp_file, recursive = TRUE)
    return(result)
  }




#' @title Request the information of one compound from PubChem
#' @description Request the information of one compound from PubChem
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param part_index Which part you want to download. See here:
#' https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/
#' @param path default is .
#' @return xml data
#' @importFrom rvest read_html html_table html_element
#' @importFrom utils download.file
#' @export

download_pubchem_compound <-
  function(part_index = as.character(1:329),
           path = ".") {
    part_index <- as.character(part_index)
    path <- file.path(path, "data")
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    part_index <-
      match.arg(part_index)
    idx1 <- as.numeric(part_index) * 500000 + 1
    idx1 <-
      format(idx1, scientific = FALSE)
    idx2 <- (as.numeric(part_index) + 1) * 500000
    idx2 <-
      format(idx2, scientific = FALSE)
    idx1 <-
      paste0(paste(rep(0, 9 - nchar(idx1)), collapse = ""), idx1) %>%
      stringr::str_trim()
    idx2 <-
      paste0(paste0(rep(0, 9 - nchar(idx2)), collapse = ""), idx2) %>%
      stringr::str_trim()
    url <-
      paste0(
        "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/Compound_",
        idx1,
        "_",
        idx2,
        ".xml.gz"
      )

    file_name <-
      paste0("Compound_",
             idx1,
             "_",
             idx2,
             ".xml.gz")

    utils::download.file(url = url, destfile = file.path(path, file_name))
    message('Done')
  }



#'
#' #' @title Read the xml database from download_pubchem_compound function
#' #' @description Read the xml database from download_pubchem_compound function
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param file should be xml format
#' #' @param path Default is .. Should be same with download_pubchem_compound function.
#' #' @return A list
#' #' @importFrom magrittr %>%
#' #' @importFrom plyr dlply .
#' #' @importFrom readr read_delim
#' #' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' #' @importFrom tidyr pivot_wider
#' #' @importFrom purrr map
#' #' @importFrom XML xmlParse
#' #' @importFrom R.utils gunzip isGzipped
#' #' @importFrom utils untar
#' #' @importFrom xml2 read_xml
#' #' @export
#' read_pubchem_xml <-
#'   function(file,
#'            path = ".") {
#'     if (R.utils::isGzipped(file.path(path, "data", file))) {
#'       message("Uncompressing data...")
#'       R.utils::gunzip(file.path(path, "data", file))
#'       message("Done")
#'     }
#'
#'     message("Reading data, it may take a while...")
#'     result <-
#'       xml2::read_xml(stringr::str_replace(file.path(path, "data", file), "\\.gz", ""))
#'     message("Done")
#'
#'     message("Parsing data, it may take a while...")
#'     result <-
#'       XML::xmlParse(result)
#'     message("Done")
#'
#'     result <-
#'       XML::xmlToList(result)
#'
#'     message("Organizing...")
#'     pb <- progress::progress_bar$new(total = length(lipidmaps))
#'
#'     lipidmaps_result <-
#'       seq_len(length(lipidmaps)) %>%
#'       purrr::map(function(i) {
#'         # cat(i, " ")
#'         pb$tick()
#'         x <- lipidmaps[[i]]
#'         result <-
#'           tryCatch(
#'             matrix(x[[4]], nrow = 1) %>%
#'               as.data.frame(),
#'             error = NULL
#'           )
#'         if (is.null(result)) {
#'           return(NULL)
#'         }
#'         colnames(result) <- names(x[[4]])
#'         result
#'       })
#'     message("Done.")
#'     return(lipidmaps_result)
#'   }
#'

#' @title Parse the json or xml compound data from PubChem
#' @description Parse the json or xml compound data from PubChem
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file_name file name of the data.
#' @return A data frame.
#' @importFrom dplyr filter
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
parse_pubchem_compound <-
  function(file_name) {
    if(requireNamespace("rjson", quietly = TRUE)){
      result <-
        rjson::fromJSON(file = file_name)
    }else{
      stop("Please install rjson package first.")
    }

    RecordType <-
      result$Record$RecordType

    RecordNumber <-
      result$Record$RecordNumber

    RecordTitle <-
      result$Record$RecordTitle

    section <-
      result$Record$Section

    section_tocheading <-
      section %>%
      lapply(function(x) {
        x$TOCHeading
      }) %>%
      unlist()

    section_description <-
      section %>%
      lapply(function(x) {
        x$Description
      })

    ###Names and Identifiers
    names <-
      section[[which(section_tocheading == "Names and Identifiers")]]$Section

    ###Computed Descriptors
    names_section_tocheading <-
      names %>%
      lapply(function(x) {
        x$TOCHeading
      }) %>%
      unlist()

    computed_descriptors <-
      names[[which(names_section_tocheading == "Computed Descriptors")]]

    id_class <-
      computed_descriptors$Section %>%
      lapply(function(x) {
        x$TOCHeading
      }) %>%
      unlist()

    id_value <-
      computed_descriptors$Section %>%
      lapply(function(x) {
        x$Information[[1]]$Value$StringWithMarkup[[1]]$String
      }) %>%
      unlist()

    names(id_value) <- id_class


    ####molecular formula
    molecular_formula <-
      names[[which(names_section_tocheading == "Molecular Formula")]]
    molecular_formula <-
      molecular_formula$Information[[1]]$Value$StringWithMarkup[[1]]$String

    ####Other Identifiers
    other_identifiers <-
      names[[which(names_section_tocheading == "Other Identifiers")]]

    other_identifiers_id_name <-
      other_identifiers$Section %>%
      lapply(function(x) {
        x$TOCHeading
      }) %>%
      unlist()

    other_identifiers_id_value <-
      other_identifiers$Section %>%
      lapply(function(x) {
        x$Information[[1]]$Value$StringWithMarkup[[1]]$String
      }) %>%
      unlist()

    names(other_identifiers_id_value) <-
      other_identifiers_id_name

    ##Synonyms
    ####Other Identifiers
    synonyms <-
      names[[which(names_section_tocheading == "Synonyms")]]

    synonyms_id_name <-
      synonyms$Section %>%
      lapply(function(x) {
        x$TOCHeading
      }) %>%
      unlist()

    synonyms <-
      synonyms$Section[which(synonyms_id_name == "Depositor-Supplied Synonyms")][[1]]$Information[[1]]$Value

    synonyms <-
      synonyms[[1]] %>%
      lapply(function(y)
        y$String) %>%
      unlist() %>%
      paste(collapse = "{}")

    result <- list(
      Computed_Descriptors = id_value,
      Molecular_Formula = molecular_formula,
      Other_Identifiers = other_identifiers_id_value,
      Synonyms = synonyms
    )
  }
