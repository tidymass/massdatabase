#' @title Request all compound information in KEGG
#' @description Request all compound information in KEGG
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @return A data frame.
#' @importFrom KEGGREST keggGet keggList
#' @importFrom magrittr %>%
#' @export
#' @examples
#' head(request_kegg_compound_info, 3)

request_kegg_compound_info <-
  function() {
    compound_id <-
      KEGGREST::keggList(database = "compound")
    result <-
      data.frame(KEGG.ID = names(compound_id),
                 Synonyms = unname(compound_id)) %>%
      dplyr::mutate(KEGG.ID = stringr::str_replace(KEGG.ID, "cpd\\:", ""))

    return(result)
  }

#' @title Request all drug information in KEGG
#' @description Request all drug information in KEGG
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @return A data frame.
#' @importFrom KEGGREST keggGet keggList
#' @importFrom magrittr %>%
#' @export
#' @examples
#' head(request_kegg_drug_info, 3)

request_kegg_drug_info <-
  function() {
    drug_id <-
      KEGGREST::keggList(database = "drug")
    result <-
      data.frame(KEGG.ID = names(drug_id),
                 Synonyms = unname(drug_id)) %>%
      dplyr::mutate(KEGG.ID = stringr::str_replace(KEGG.ID, "dr\\:", ""))

    return(result)
  }





#' @title Request one specific compound information in KEGG
#' @description Request one specific compound information in KEGG
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param compound_id compound id. For example, C02886
#' @param return_form data.frame or list.
#' @return A data frame or list.
#' @importFrom KEGGREST keggGet
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_kegg_compound(compound_id = "C02886", return_form = "list")
#' x[1:2]
#' y = request_kegg_compound(compound_id = "C02886", return_form = "data.frame")
#' head(y)

request_kegg_compound <-
  function(compound_id = "C02886",
           return_form = c("list" , "data.frame")) {
    return_form <- match.arg(return_form)

    x <-
      KEGGREST::keggGet(dbentries = compound_id)[[1]]

    if (return_form == "list") {
      result <- x
    } else{
      KEGG.ID = x$ENTRY
      x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
      Compound.name = paste(x$NAME, collapse = "{}")
      Formula = x$FORMULA
      if (is.null(x$FORMULA)) {
        Formula = NA
      }
      mz = as.numeric(x$EXACT_MASS)
      if (is.null(x$EXACT_MASS)) {
        mz = NA
      }

      CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
        stringr::str_trim(side = "both")

      PUBCHEM.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
        stringr::str_trim(side = "both")

      CHEBI.ID <-
        stringr::str_replace(grep("ChEBI", x$DBLINKS, value = TRUE), "ChEBI: ", "") %>%
        stringr::str_trim(side = "both")

      CHEMBL.ID <-
        stringr::str_replace(grep("ChEMBL", x$DBLINKS, value = TRUE), "ChEMBL: ", "") %>%
        stringr::str_trim(side = "both")

      LIPIDMAPS.ID <-
        stringr::str_replace(grep("LIPIDMAPS", x$DBLINKS, value = TRUE), "LIPIDMAPS: ", "") %>%
        stringr::str_trim(side = "both")

      LIPIDBANK.ID <-
        stringr::str_replace(grep("LipidBank", x$DBLINKS, value = TRUE), "LipidBank: ", "") %>%
        stringr::str_trim(side = "both")

      DRUGBANK.ID <-
        stringr::str_replace(grep("DrugBank", x$DBLINKS, value = TRUE), "DrugBank: ", "") %>%
        stringr::str_trim(side = "both")

      if (length(CAS.ID) == 0) {
        CAS.ID = NA
      }

      if (length(PUBCHEM.ID) == 0) {
        PUBCHEM.ID = NA
      }

      if (length(CHEBI.ID) == 0) {
        CHEBI.ID = NA
      }

      if (length(CHEMBL.ID) == 0) {
        CHEMBL.ID = NA
      }

      if (length(LIPIDMAPS.ID) == 0) {
        LIPIDMAPS.ID = NA
      }

      if (length(LIPIDBANK.ID) == 0) {
        LIPIDBANK.ID = NA
      }

      if (length(DRUGBANK.ID) == 0) {
        DRUGBANK.ID = NA
      }

      From_human = "Yes"
      REMARK <- x$REMARK
      if (is.null(REMARK)) {
        From_drug <- "No"
        KEGG_DRUG.ID <- NA
      } else{
        KEGG_DRUG.ID <-
          paste(stringr::str_extract_all(REMARK, "D[0-9]{5,6}")[[1]],
                collapse = "{}")

        if (length(KEGG_DRUG.ID) == 0) {
          From_drug <- "No"
          KEGG_DRUG.ID <- NA
        } else{
          From_drug <- "Yes"
          KEGG_DRUG.ID <- KEGG_DRUG.ID
        }
      }

      result <-
        data.frame(
          Lab.ID = KEGG.ID,
          Compound.name,
          Formula,
          mz,
          CAS.ID,
          HMDB.ID = NA,
          KEGG.ID,
          PUBCHEM.ID,
          CHEBI.ID,
          CHEMBL.ID,
          LIPIDMAPS.ID,
          LIPIDBANK.ID,
          DRUGBANK.ID = DRUGBANK.ID,
          From_human = From_human,
          From_drug = From_drug,
          KEGG_DRUG.ID = KEGG_DRUG.ID
        ) %>%
        dplyr::mutate(Synonyms = Compound.name)
    }
    return(result)
  }



#' @title Request one specific drug information in KEGG
#' @description Request one specific drug information in KEGG
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param drug_id drug id. For example, D00001
#' @param return_form data.frame or list.
#' @return A data frame or list.
#' @importFrom KEGGREST keggGet
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_kegg_drug(drug_id = "D00001", return_form = "list")
#' x[1:2]
#' y = request_kegg_drug(drug_id = "D00001", return_form = "data.frame")
#' head(y)

request_kegg_drug <-
  function(drug_id = "D00001",
           return_form = c("list" , "data.frame")) {
    return_form <- match.arg(return_form)
    request_kegg_compound(compound_id = drug_id,
                          return_form = return_form)
  }

#' @title Download KEGG compound data
#' @description Download KEGG compound data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .
#' @param sleep Default is 1 second.
#' @return KEGG compound database, rda format.
#' @export

download_kegg_compound <-
  function(path = ".",
           sleep = 1) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    kegg_id <-
      request_kegg_compound_info()

    pb <- progress::progress_bar$new(total = nrow(kegg_id))

    kegg_compound_database <-
      seq_along(kegg_id$KEGG.ID) %>%
      purrr::map(function(i) {
        pb$tick()
        Sys.sleep(time = sleep)
        KEGGREST::keggGet(dbentries = kegg_id$KEGG.ID[i])[[1]]
      })
    save(kegg_compound_database,
         file = file.path(path, "kegg_compound_database"))
  }


#' @title Download KEGG drug data
#' @description Download KEGG drug data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .
#' @param sleep Default is 1 second.
#' @return KEGG drug database, rda format.
#' @export

download_kegg_drug <-
  function(path = ".",
           sleep = 1) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    kegg_id <-
      request_kegg_drug_info()

    pb <- progress::progress_bar$new(total = nrow(kegg_id))

    kegg_drug_database <-
      seq_along(kegg_id$KEGG.ID) %>%
      purrr::map(function(i) {
        pb$tick()
        Sys.sleep(time = sleep)
        KEGGREST::keggGet(dbentries = kegg_id$KEGG.ID[i])[[1]]
      })
    save(kegg_drug_database,
         file = file.path(path, "kegg_drug_database"))
  }


#' @title Read the KEGG compound database from download_kegg_compound function
#' @description Read the KEGG compound database from download_kegg_compound function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with download_kegg_compound function.
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @export
read_kegg_compound <-
  function(path = ".") {
    load(file.path(path, "kegg_compound_database"))

    pb <-
      progress::progress_bar$new(total = length(kegg_compound_database))

    kegg_metabolite =
      seq_len(length(kegg_compound_database)) %>%
      purrr::map(function(i) {
        pb$tick()
        x <- kegg_compound_database[[i]]
        KEGG.ID = x$ENTRY
        x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
        Compound.name = paste(x$NAME, collapse = "{}")
        Formula = x$FORMULA
        if (is.null(x$FORMULA)) {
          Formula = NA
        }
        mz = as.numeric(x$EXACT_MASS)
        if (is.null(x$EXACT_MASS)) {
          mz = NA
        }

        CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
          stringr::str_trim(side = "both")

        PUBCHEM.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
          stringr::str_trim(side = "both")

        CHEBI.ID <-
          stringr::str_replace(grep("ChEBI", x$DBLINKS, value = TRUE), "ChEBI: ", "") %>%
          stringr::str_trim(side = "both")

        CHEMBL.ID <-
          stringr::str_replace(grep("ChEMBL", x$DBLINKS, value = TRUE), "ChEMBL: ", "") %>%
          stringr::str_trim(side = "both")

        LIPIDMAPS.ID <-
          stringr::str_replace(grep("LIPIDMAPS", x$DBLINKS, value = TRUE), "LIPIDMAPS: ", "") %>%
          stringr::str_trim(side = "both")

        LIPIDBANK.ID <-
          stringr::str_replace(grep("LipidBank", x$DBLINKS, value = TRUE), "LipidBank: ", "") %>%
          stringr::str_trim(side = "both")

        DRUGBANK.ID <-
          stringr::str_replace(grep("DrugBank", x$DBLINKS, value = TRUE), "DrugBank: ", "") %>%
          stringr::str_trim(side = "both")

        if (length(CAS.ID) == 0) {
          CAS.ID = NA
        }

        if (length(PUBCHEM.ID) == 0) {
          PUBCHEM.ID = NA
        }

        if (length(CHEBI.ID) == 0) {
          CHEBI.ID = NA
        }

        if (length(CHEMBL.ID) == 0) {
          CHEMBL.ID = NA
        }

        if (length(LIPIDMAPS.ID) == 0) {
          LIPIDMAPS.ID = NA
        }

        if (length(LIPIDBANK.ID) == 0) {
          LIPIDBANK.ID = NA
        }

        if (length(DRUGBANK.ID) == 0) {
          DRUGBANK.ID = NA
        }

        From_human = "Yes"
        REMARK <- x$REMARK
        if (is.null(REMARK)) {
          From_drug <- "No"
          KEGG_DRUG.ID <- NA
        } else{
          KEGG_DRUG.ID <-
            paste(stringr::str_extract_all(REMARK, "D[0-9]{5,6}")[[1]],
                  collapse = "{}")

          if (length(KEGG_DRUG.ID) == 0) {
            From_drug <- "No"
            KEGG_DRUG.ID <- NA
          } else{
            From_drug <- "Yes"
            KEGG_DRUG.ID <- KEGG_DRUG.ID
          }
        }

        data.frame(
          Lab.ID = KEGG.ID,
          Compound.name,
          Formula,
          mz,
          CAS.ID,
          HMDB.ID = NA,
          KEGG.ID,
          PUBCHEM.ID,
          CHEBI.ID,
          CHEMBL.ID,
          LIPIDMAPS.ID,
          LIPIDBANK.ID,
          DRUGBANK.ID = DRUGBANK.ID,
          From_human = From_human,
          From_drug = From_drug,
          KEGG_DRUG.ID = KEGG_DRUG.ID
        )
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    return(kegg_metabolite)
  }


#' @title Read the KEGG drug database from download_kegg_drug function
#' @description Read the KEGG drug database from download_kegg_drug function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with download_kegg_drug function.
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @export
read_kegg_drug <-
  function(path = ".") {
    load(file.path(path, "kegg_drug_database"))

    pb <-
      progress::progress_bar$new(total = length(kegg_drug_database))

    kegg_metabolite =
      seq_len(length(kegg_drug_database)) %>%
      purrr::map(function(i) {
        pb$tick()
        x <- kegg_drug_database[[i]]
        KEGG.ID = x$ENTRY
        x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
        Compound.name = paste(x$NAME, collapse = "{}")
        Formula = x$FORMULA
        if (is.null(x$FORMULA)) {
          Formula = NA
        }
        mz = as.numeric(x$EXACT_MASS)
        if (is.null(x$EXACT_MASS)) {
          mz = NA
        }

        CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
          stringr::str_trim(side = "both")

        PUBCHEM.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
          stringr::str_trim(side = "both")

        CHEBI.ID <-
          stringr::str_replace(grep("ChEBI", x$DBLINKS, value = TRUE), "ChEBI: ", "") %>%
          stringr::str_trim(side = "both")

        CHEMBL.ID <-
          stringr::str_replace(grep("ChEMBL", x$DBLINKS, value = TRUE), "ChEMBL: ", "") %>%
          stringr::str_trim(side = "both")

        LIPIDMAPS.ID <-
          stringr::str_replace(grep("LIPIDMAPS", x$DBLINKS, value = TRUE), "LIPIDMAPS: ", "") %>%
          stringr::str_trim(side = "both")

        LIPIDBANK.ID <-
          stringr::str_replace(grep("LipidBank", x$DBLINKS, value = TRUE), "LipidBank: ", "") %>%
          stringr::str_trim(side = "both")

        DRUGBANK.ID <-
          stringr::str_replace(grep("DrugBank", x$DBLINKS, value = TRUE), "DrugBank: ", "") %>%
          stringr::str_trim(side = "both")

        if (length(CAS.ID) == 0) {
          CAS.ID = NA
        }

        if (length(PUBCHEM.ID) == 0) {
          PUBCHEM.ID = NA
        }

        if (length(CHEBI.ID) == 0) {
          CHEBI.ID = NA
        }

        if (length(CHEMBL.ID) == 0) {
          CHEMBL.ID = NA
        }

        if (length(LIPIDMAPS.ID) == 0) {
          LIPIDMAPS.ID = NA
        }

        if (length(LIPIDBANK.ID) == 0) {
          LIPIDBANK.ID = NA
        }

        if (length(DRUGBANK.ID) == 0) {
          DRUGBANK.ID = NA
        }

        From_drug = "Yes"
        REMARK <- x$REMARK
        if (is.null(REMARK)) {
          From_human <- "No"
          KEGG.ID <- NA
        } else{
          KEGG.ID <-
            paste(stringr::str_extract_all(REMARK, "C[0-9]{4,6}")[[1]],
                  collapse = "{}")

          if (length(KEGG.ID) == 0) {
            From_human <- "No"
            KEGG.ID <- NA
          } else{
            From_human <- "Yes"
            KEGG.ID <- KEGG_DRUG.ID
          }
        }

        data.frame(
          Lab.ID = KEGG.ID,
          Compound.name,
          Formula,
          mz,
          CAS.ID,
          HMDB.ID = NA,
          KEGG.ID,
          PUBCHEM.ID,
          CHEBI.ID,
          CHEMBL.ID,
          LIPIDMAPS.ID,
          LIPIDBANK.ID,
          DRUGBANK.ID = DRUGBANK.ID,
          From_human = From_human,
          From_drug = From_drug,
          KEGG_DRUG.ID = KEGG_DRUG.ID
        )
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    return(kegg_metabolite)
  }


#' @title Convert KEGG compound/drug data (data.frame,
#' from read_kegg_compound or read_kegg_drug function)
#' to metID format database
#' @description  KEGG compound/drug data (data.frame,
#' from read_kegg_compound or read_kegg_drug function)
#' to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data data.frame, from read_kegg_compound or read_kegg_drug function
#' @param path Default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_kegg2metid <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    data <-
      data %>%
      dplyr::filter(!is.na(mz)) %>%
      dplyr::mutate(Synonyms = Compound.name) %>%
      dplyr::mutate(
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "KEGG"
      ) %>%
      dplyr::select(
        Lab.ID,
        Compound.name,
        mz,
        RT,
        CAS.ID,
        HMDB.ID,
        KEGG.ID,
        Formula,
        mz.pos,
        mz.neg,
        Submitter,
        everything()
      )

    data$Compound.name =
      data$Compound.name %>%
      stringr::str_split(pattern = "\\{\\}") %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist() %>%
      stringr::str_replace(pattern = ";", "")

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = data,
                     file = file.path(temp_file, "data.csv"))

    kegg_ms1 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "data.csv",
        source = "KEGG",
        link = "https://www.genome.jp/kegg",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "data.csv"))
    unlink(temp_file)

    save(kegg_ms1, file = file.path(path, "kegg_ms1"))
    invisible(kegg_ms1)
  }
