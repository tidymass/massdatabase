#' @title Request all the compound information in FooDB based on web crawler
#' @description Request all the compound information in FooDB based on web crawler
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
#' x = request_foodb_compound_info_crawler(pages = 1)
#' head(x)

request_foodb_compound_info_crawler <-
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
      foods <- paste(unlist(result$foods[1, ]), collapse = "{}")
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
          chemspider_id = ifelse(is.null(result$chemspider_id),
                                 NA,
                                 result$chemspider_id),
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
#' @param url Default is "https://foodb.ca/downloads".
#' @param sleep sleep default is 1 second.
#' @param pages default is from 1:2838
#' @return A data frame.
#' @importFrom rvest read_html html_table
#' @importFrom magrittr %>%
#' @importFrom utils download.file untar
#' @export
#' @examples
#' x = request_foodb_compound_info()
#' head(x)

request_foodb_compound_info <-
  function(url = c(
    "https://raw.githubusercontent.com/jaspershen/databases/main/data/foodb_compound_info.rda",
    "https://foodb.ca/downloads",
    "https://foodb.ca/compounds"
  ),
  sleep = 1,
  pages = 1:2838) {
    url <- match.arg(url)
    if (url == "https://raw.githubusercontent.com/jaspershen/databases/main/data/foodb_compound_info.rda") {
      temp_file <- tempfile()
      dir.create(temp_file, showWarnings = FALSE)
      options(timeout = 10000)
      utils::download.file(url = url,
                           destfile = file.path(temp_file, "foodb_compound_info.rda"))

      load(file.path(temp_file, "foodb_compound_info.rda"))

      unlink(temp_file)
      x <-
        food_compound_info %>%
        dplyr::select(public_id, name)
      return(x)
    }

    if (url ==  "https://foodb.ca/downloads") {
      temp_file <- tempfile()
      dir.create(temp_file, showWarnings = FALSE)
      options(timeout = 10000)
      utils::download.file(url = "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz",
                           destfile = file.path(temp_file, "foodb_2020_4_7_csv.tar"))

      utils::untar(tarfile = file.path(temp_file, "foodb_2020_4_7_csv.tar"),
                   exdir = temp_file)

      food_compound_info <-
        readr::read_csv(file.path(temp_file,
                                  "foodb_2020_04_07_csv/Compound.csv"),
                        show_col_types = FALSE)

      unlink(temp_file)
      x <-
        food_compound_info %>%
        dplyr::select(public_id, name)
      return(x)
    }

    if (url == "https://foodb.ca/compounds") {
      request_foodb_compound_info_crawler(sleep = sleep,
                                          pages = pages)
    }

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
      which(unname(unlist(spectra[1, ])) == "Specdb::MsMs")

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








#' @title Download FOODB compound data
#' @description Download FOODB compound data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param compound_id "all": download all compounds, or a vector of compound IDs.
#' @param path Default is .
#' @return FOODB compound database, rda format.
#' @export

download_foodb_compound <-
  function(compound_id = "all",
           path = ".") {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    foodb_id <-
      request_foodb_compound_info(url =
                                    "https://raw.githubusercontent.com/jaspershen/databases/main/data/foodb_compound_info.rda")

    if (all(compound_id != "all")) {
      foodb_id <-
        foodb_id %>%
        dplyr::filter(public_id %in% compound_id)
    }

    foodb_id <-
      foodb_id$public_id

    pb <- progress::progress_bar$new(total = length(foodb_id))

    foodb_compound_database <-
      seq_along(foodb_id) %>%
      purrr::map(function(i) {
        pb$tick()
        result <-
          tryCatch(
            request_foodb_compound(compound_id = foodb_id[i]),
            error = function(e)
              NULL
          )
        if (is.null(result)) {
          return(NULL)
        } else{
          return(result)
        }
      })

    save(foodb_compound_database,
         file = file.path(path, "foodb_compound_database"))
  }






#' @title Read the FOODB compound database from download_foodb_compound function
#' @description Read the FOODB compound database from download_foodb_compound function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with download_foodb_compound function.
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @export
read_foodb_compound <-
  function(path = ".") {
    load(file.path(path, "foodb_compound_database"))

    pb <-
      progress::progress_bar$new(total = length(foodb_compound_database))

    foodb_result <-
      seq_along(foodb_compound_database) %>%
      purrr::map(function(i) {
        # cat(x$accession, " ")
        pb$tick()
        x <- foodb_compound_database[[i]]
        if (is.null(x)) {
          return(NULL)
        }
        Kingdom <-
          ifelse(is.null(x$taxonomy$kingdom), NA, x$taxonomy$kingdom)
        Super_class <-
          ifelse(is.null(x$taxonomy$super_class),
                 NA,
                 x$taxonomy$super_class)
        Class <-
          ifelse(is.null(x$taxonomy$class), NA, x$taxonomy$class)
        Sub_class <-
          ifelse(is.null(x$taxonomy$sub_class), NA, x$taxonomy$sub_class)
        State <- ifelse(is.null(x$state), NA, x$state)
        HMDB.ID <- ifelse(is.null(x$hmdb_id), NA, x$hmdb_id)
        PUBCHEM.ID <-
          ifelse(is.null(x$pubchem_compound_id),
                 NA,
                 x$pubchem_compound_id)
        CHEMSPIDER.ID <-
          ifelse(is.null(x$chemspider_id), NA, x$chemspider_id)
        KEGG.ID <- ifelse(is.null(x$kegg_id), NA, x$kegg_id)
        CHEBI.ID <- ifelse(is.null(x$chebi_id), NA, x$chebi_id)
        BIOCYC.ID <- ifelse(is.null(x$biocyc_id), NA, x$biocyc_id)
        HET.ID <- ifelse(is.null(x$het_id), NA, x$het_id)
        WIKIPEDIA.ID <-
          ifelse(is.null(x$wikipidia), NA, x$wikipidia)
        VMH.ID <- ifelse(is.null(x$vmh_id), NA, x$vmh_id)
        Synonyms = ifelse(is.null(x$synonyms), NA, paste(unname(unlist(x$synonyms)), collapse = "{}"))
        foods <- x$foods
        if (is.null(foods)) {
          Food_name <- NA
          Food_type <- NA
          Food_category <- NA
          Food_scientific_name <- NA
        } else{
          if (class(foods)[1] == "list") {
            foods <-
              foods %>%
              lapply(function(y) {
                data.frame(
                  name = ifelse(is.null(y$name), NA, y$name),
                  food_type = ifelse(is.null(y$food_type), NA, y$food_type),
                  category = ifelse(is.null(y$category), NA, y$category),
                  name_scientific = ifelse(is.null(y$name_scientific), NA, y$name_scientific)
                )
              }) %>%
              dplyr::bind_rows() %>%
              as.data.frame()

            Food_name <- as.character(foods[, "name"]) %>%
              paste(collapse = "{}")
            Food_type <- as.character(foods[, "food_type"]) %>%
              paste(collapse = "{}")
            Food_category <- as.character(foods[, "category"]) %>%
              paste(collapse = "{}")
            Food_scientific_name <-
              as.character(foods[, "name_scientific"]) %>%
              paste(collapse = "{}")

          } else{
            colnames(foods) <- paste("V", 1:ncol(foods), sep = "")

            foods <-
              foods[1:4, , drop = FALSE] %>%
              as.data.frame() %>%
              apply(2, function(y) {
                y <-
                  y %>%
                  lapply(function(z) {
                    if (is.null(z)) {
                      return(NA)
                    } else{
                      return(z)
                    }
                  }) %>%
                  unlist()
                y
              })
            Food_name <- as.character(foods["name",]) %>%
              paste(collapse = "{}")
            Food_type <- as.character(foods["food_type",]) %>%
              paste(collapse = "{}")
            Food_category <- as.character(foods["category",]) %>%
              paste(collapse = "{}")
            Food_scientific_name <-
              as.character(foods["name_scientific",]) %>%
              paste(collapse = "{}")
          }
        }

        data.frame(
          Lab.ID = x$accession,
          Create_date = ifelse(is.null(x$creation_date), NA, x$creation_date),
          Updated_date = ifelse(is.null(x$update_date), NA, x$update_date),
          Compound.name = ifelse(is.null(x$name), NA, x$name),
          Description = ifelse(is.null(x$description), NA, x$description),
          Formula = ifelse(is.null(x$chemical_formula), NA, x$chemical_formula),
          Synonyms = Synonyms,
          Average.mass =  ifelse(
            is.null(x$average_molecular_weight),
            NA,
            x$average_molecular_weight
          ),
          mz = ifelse(
            is.null(x$monisotopic_moleculate_weight),
            NA,
            x$monisotopic_moleculate_weight
          ),
          IUPAC_name = ifelse(is.null(x$iupac_name), NA, x$iupac_name),
          Traditional_IUPAC_name = ifelse(is.null(x$traditional_iupac), NA, x$traditional_iupac),
          CAS.ID = ifelse(
            is.null(x$cas_registry_number),
            NA,
            x$cas_registry_number
          ),
          SMILES.ID = ifelse(is.null(x$smiles), NA, x$smiles),
          INCHI.ID = ifelse(is.null(x$inchi), NA, x$inchi),
          INCHIKEY.ID = ifelse(is.null(x$inchikey), NA, x$inchikey),
          Kingdom = Kingdom,
          Super_class = Super_class,
          Class = Class,
          Sub_class = Sub_class,
          State = State,
          FOODB.ID = ifelse(is.null(x$accession), NA, x$accession),
          HMDB.ID = HMDB.ID,
          PUBCHEM.ID = PUBCHEM.ID,
          CHEMSPIDER.ID = CHEMSPIDER.ID,
          KEGG.ID = KEGG.ID,
          CHEBI.ID = CHEBI.ID,
          BIOCYC.ID = BIOCYC.ID,
          HET.ID = HET.ID,
          WIKIPEDIA.ID = WIKIPEDIA.ID,
          VMH.ID = VMH.ID,
          Food_name,
          Food_type,
          Food_category,
          Food_scientific_name
        )
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    foodb_result$mz <-
      as.numeric(foodb_result$mz)

    foodb_result$RT <- NA
    foodb_result$From_food <- TRUE
    foodb_result$mz.pos = NA
    foodb_result$mz.neg = NA
    foodb_result$Submitter = "FOODB"

    foodb_result <-
      foodb_result %>%
      dplyr::filter(!is.na(mz) & !is.na(Formula))

    return(foodb_result)
  }



#' @title Convert FOODB compound data (list,
#' from download_foodb_compound function)
#' to metID format database
#' @description Convert FOODB compound data (list,
#' from download_foodb_compound function)
#' to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data data.frame, from read_foodb_compound function.
#' @param path Default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_foodb2metid <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = data,
                     file = file.path(temp_file, "data.csv"))

    foodb_ms1 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "data.csv",
        source = "FOODB",
        link = "https://foodb.ca/",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "data.csv"))
    unlink(temp_file)

    save(foodb_ms1, file = file.path(path, "foodb_ms1"))
    invisible(foodb_ms1)
  }
