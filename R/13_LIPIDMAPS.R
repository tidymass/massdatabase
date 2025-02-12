#' @title Download LIPIDMAPS data
#' @description Download LIPIDMAPS data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .
#' @return version.
#' @export

download_lipidmaps_lipid <-
  function(path = ".") {
    download.file(url = "https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip", destfile = file.path(path, "LMSD.sdf.zip"))
  }


#' @title Read the sdf database from download_lipidmaps_lipid function
#' @description Read the sdf database from download_lipidmaps_lipid function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file should be sdf format
#' @param path Default is .. Should be same with download_lipidmaps_lipid function.
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @export
read_sdf_data_lipidmaps <-
  function(file, path = ".") {
    message("Reading data, it may take a while...")
    if (requireNamespace("ChemmineR", quietly = TRUE)) {
      lipidmaps <-
        ChemmineR::read.SDFset(sdfstr = file.path(path, file), skipErrors = TRUE)
    } else{
      stop('Please install ChemmineR package first.')
    }
    message("Done.")
    message("Organizing...")
    pb <- progress::progress_bar$new(total = length(lipidmaps))

    lipidmaps_result <-
      seq_len(length(lipidmaps)) %>%
      purrr::map(function(i) {
        # cat(i, " ")
        pb$tick()
        x <- lipidmaps[[i]]
        result <-
          tryCatch(
            matrix(x[[4]], nrow = 1) %>%
              as.data.frame(),
            error = NULL
          )
        if (is.null(result)) {
          return(NULL)
        }
        colnames(result) <- names(x[[4]])
        result
      })
    message("Done.")
    return(lipidmaps_result)
  }


#' @title Request one lipid information in LIPIDMAPS
#' @description Request one lipid information in LIPIDMAPS
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://www.lipidmaps.org/databases/lmsd".
#' @param lipid_id lipid id. For example, LMFA01030001
#' @return A data frame.
#' @importFrom magrittr %>%
#' @importFrom rvest html_node html_text2 html_nodes
#' @importFrom utils tail
#' @export
#' @examples
#' request_lipidmaps_lipid(lipid_id = "LMFA01030001")
#' request_lipidmaps_lipid(lipid_id = "LMGL00000140")




request_lipidmaps_lipid <-
  function(url = "https://www.lipidmaps.org/databases/lmsd", lipid_id = "LMFA01030001") {
    result <-
      tryCatch(
        rvest::read_html(paste0(url, "/", lipid_id)),
        error = function(e) {
          NULL
        }
      )

    if (is.null(result)) {
      message("Check this url: ", paste0(url, "/", lipid_id))
      return(NULL)
    }


    ###base information
    # idx <- grep("Show lipids differing", info)
    # info <- info[1:idx - 1]
    base_info <-
      result %>%
      rvest::html_node(".block .px-5") %>%
      # html_nodes("div") %>%
      html_text2() %>%
      stringr::str_split("\n") %>%
      `[[`(1)

    base_info <- base_info[-grep("Select m/z", base_info)]

    base_info <-
      base_info[!base_info %in% c("Calculate m/z")]

    base_info <-
      base_info[!base_info %in% ""]

    idx1 <-
      which(base_info == "Synonyms")
    idx2 <-
      which(base_info == "LM ID")

    if (idx2 - idx1 == 1) {
      Synonyms <- NA
    } else{
      Synonyms <-
        paste(base_info[(idx1 + 1):(idx2 - 1)], collapse = "{}")
    }

    Common_Name <-
      base_info[which(base_info == "Common Name") + 1]
    Common_Name <-
      ifelse(length(Common_Name) == 0, NA, Common_Name)

    Systematic_Name <-
      base_info[which(base_info == "Systematic Name") + 1]

    Systematic_Name <-
      ifelse(length(Systematic_Name) == 0, NA, Systematic_Name)

    LM_ID <-
      base_info[which(base_info == "LM ID") + 1]
    LM_ID <-
      ifelse(length(LM_ID) == 0, NA, LM_ID)

    Status <-
      base_info[which(base_info == "Status") + 1]
    Status <-
      ifelse(length(Status) == 0, NA, Status)

    Exact_Mass <-
      base_info[which(base_info == "Exact Mass") + 1] %>%
      as.numeric()
    Exact_Mass <-
      ifelse(length(Exact_Mass) == 0, NA, Exact_Mass)

    Formula <-
      base_info[which(base_info == "Formula") + 1]
    Formula <-
      ifelse(length(Formula) == 0, NA, Formula)

    Abbrev <-
      base_info[which(base_info == "Abbrev") + 1]
    Abbrev <-
      ifelse(length(Abbrev) == 0, NA, Abbrev)

    all_other_results <-
      1:20 %>%
      purrr::map(function(i) {
        node <-
          paste(".rounded-lg:nth-child(", i, ") .px-5", sep = "")
        result %>%
          html_node(node) %>%
          html_nodes("div") %>%
          html_text2() %>%
          trimws()
      })

    remove_index <-
      all_other_results %>%
      lapply(function(x) {
        length(x) == 0
      }) %>%
      unlist() %>%
      which()

    if (length(remove_index) > 0) {
      all_other_results <- all_other_results[-remove_index]
    }

    ####classification
    classification_index <-
      lapply(all_other_results, function(x) {
        any(stringr::str_detect(x, "Category"))
      }) %>%
      unlist() %>%
      which()

    if (length(classification_index) == 0) {
      classification <-
        data.frame(
          Category = NA,
          Main_Class = NA,
          Sub_Class = NA,
          stringsAsFactors = FALSE
        )
    } else{
      classification <-
        all_other_results[[classification_index]]

      classification <- data.frame(
        Category = classification[5],
        # First value after labels
        Main_Class = classification[6],
        Sub_Class = classification[7],
        stringsAsFactors = FALSE
      )
    }


    ####Biological Context
    biological_context_index <-
      lapply(all_other_results, function(x) {
        any(stringr::str_detect(x, "Biological Context"))
      }) %>%
      unlist() %>%
      which()

    if (length(biological_context_index) == 0) {
      biological_context <-
        data.frame(Biological_Context = NA,
                   stringsAsFactors = FALSE)
    } else{
      biological_context <-
        all_other_results[[biological_context_index]]

      biological_context <- data.frame(Biological_Context = biological_context[1],
                                       stringsAsFactors = FALSE)
    }


    ##References
    taxonomy_information <-
      result %>%
      html_node(".gap-4") %>%
      html_nodes("div") %>%
      html_text2()

    if (length(taxonomy_information) > 0) {
      taxonomy_information <-
        taxonomy_information[!taxonomy_information %in% c("Reference")]
      taxonomy_information <-
        matrix(taxonomy_information[1:4], nrow = 2, byrow = TRUE) %>%
        as.data.frame()

      colnames(taxonomy_information) <-
        as.character(taxonomy_information[1, ])

      taxonomy_information <-
        taxonomy_information[-1, , drop = FALSE]
    } else{
      taxonomy_information <-
        data.frame(
          "Curated from" = NA,
          "NCBI taxonomy class" = NA,
          check.names = FALSE
        )
    }



    ###String Representations
    string_representations_index <-
      lapply(all_other_results, function(x) {
        any(stringr::str_detect(x, "InChiKey \\(Click to copy\\)"))
      }) %>%
      unlist() %>%
      which()

    if (length(string_representations_index) == 0) {
      InChiKey = NA
      InChi = NA
      SMILES = NA
    } else{
      string_representations <-
        all_other_results[[string_representations_index]]

      InChiKey <-
        string_representations[which(string_representations == "InChiKey (Click to copy)") + 1]

      InChi <-
        string_representations[which(string_representations == "InChi (Click to copy)") + 1]

      SMILES <-
        string_representations[which(string_representations == "SMILES (Click to copy)") + 1]

    }




    ##Other databases
    other_database_index <-
      lapply(all_other_results, function(x) {
        any(
          stringr::str_detect(
            x,
            "Wikipedia|KEGG ID|HMDB ID|CHEBI ID|LIPIDBANK ID|PubChem CID|PlantFA ID|SwissLipids ID|Cayman ID|PDB ID|GuidePharm ID"
          )
        )
      }) %>%
      unlist() %>%
      which()


    if (length(other_database_index) == 0) {
      other_databases <-
        data.frame(
          Wikipedia = NA,
          `KEGG ID` = NA,
          `HMDB ID` = NA,
          `CHEBI ID` = NA,
          `LIPIDBANK ID` = NA,
          `PubChem CID` = NA,
          `PlantFA ID` = NA,
          `SwissLipids ID` = NA,
          `Cayman ID` = NA,
          `PDB ID` = NA,
          `GuidePharm ID` = NA,
          check.names = FALSE
        )
    } else{
      other_databases <-
        all_other_results[[other_database_index]]

      other_databases <-
        other_databases[-1]

      other_databases <-
        other_databases %>%
        matrix(nrow = 2, byrow = FALSE) %>%
        as.data.frame()

      colnames(other_databases) <-
        as.character(other_databases[1, ])

      other_databases <-
        other_databases[-1, , drop = FALSE]
    }



    return_result <-
      data.frame(
        Synonyms,
        Common_Name,
        Systematic_Name,
        LM_ID,
        Status,
        Exact_Mass,
        Formula,
        Abbrev,
        classification,
        InChiKey,
        InChi,
        SMILES,
        taxonomy_information,
        other_databases
      )

    return_result[which(return_result == "-", arr.ind = TRUE)] <- NA
    return_result[which(return_result == "", arr.ind = TRUE)] <- NA
    return_result
  }



#' @title Convert LIPIDMAPS compound data (list,
#' from read_sdf_data_lipidmaps function)
#' to metID format database
#' @description Convert LIPIDMAPS compound data (list,
#' from read_sdf_data_lipidmaps function)
#' to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_sdf_data_lipidmaps function
#' @param path Default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_lipidmaps2metid <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    all_column_name <-
      data %>%
      lapply(colnames) %>%
      unlist() %>%
      unique()

    message("Organizing...")
    pb <-
      progress::progress_bar$new(total = length(data))

    lipidmaps_ms1 <-
      data <-
      seq_len(length(data)) %>%
      purrr::map(function(i) {
        # cat(i, " ")
        pb$tick()
        x <- data[[i]]

        diff_names <- setdiff(all_column_name, colnames(x))
        if (length(diff_names) == 0) {
          return(x[, all_column_name])
        }
        add_info <-
          matrix(NA, nrow = 1, ncol = length(diff_names)) %>%
          as.data.frame()
        colnames(add_info) <- diff_names
        cbind(x, add_info) %>%
          dplyr::select(all_column_name)
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    lipidmaps_ms1 <-
      lipidmaps_ms1 %>%
      dplyr::rename(
        Lab.ID = LM_ID,
        PUBCHEM.ID = PUBCHEM_CID,
        Compound.name = NAME,
        Systematic.name = SYSTEMATIC_NAME,
        Main_class_lipidmaps = MAIN_CLASS,
        Sub_class_lipidmaps = SUB_CLASS,
        Abbreviation = ABBREVIATION,
        Class_level4 = CLASS_LEVEL4,
        HMDB.ID = HMDB_ID,
        Category = CATEGORY,
        Formula = FORMULA,
        mz = EXACT_MASS,
        KEGG.ID = KEGG_ID,
        CHEBI.ID = CHEBI_ID,
        SWISSLIPIDS.ID = SWISSLIPIDS_ID,
        LIPIDBANK.ID = LIPIDBANK_ID,
        PLANTFA.ID = PLANTFA_ID,
        SMILES.ID = SMILES,
        INCHI.ID = INCHI,
        INCHIKEY.ID = INCHI_KEY,
        Synonyms = SYNONYMS
      ) %>%
      dplyr::mutate(
        LIPIDMAPS.ID = Lab.ID,
        CAS.ID = NA,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "LIPIDMAPS"
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

    lipidmaps_ms1$mz <-
      as.numeric(lipidmaps_ms1$mz)

    lipidmaps_ms1$Synonyms <-
      lipidmaps_ms1$Synonyms %>%
      stringr::str_replace_all('; ', "{}")

    lipidmaps_ms1[which(data == "", arr.ind = TRUE)] <-
      NA

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = lipidmaps_ms1,
                     file = file.path(temp_file, "lipidmaps_ms1.csv"))

    lipidmaps_ms1 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "lipidmaps_ms1.csv",
        source = "LIPIDMAPS",
        link = "https://www.lipidmaps.org/",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "lipidmaps_ms1.csv"))
    unlink(temp_file)

    save(lipidmaps_ms1, file = file.path(path, "lipidmaps_ms1"))
    invisible(lipidmaps_ms1)
  }
