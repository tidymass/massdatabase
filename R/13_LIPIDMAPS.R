#' @title Download LIPIDMAPS data
#' @description Download LIPIDMAPS data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .
#' @return version.
#' @export

download_lipidmaps_lipid <-
  function(path = ".") {
    url <-
      paste0('http://bigg.ucsd.edu/static/models/',
             model_id,
             '.json')
    download.file(url = "https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip",
                  destfile = file.path(path, "LMSD.sdf.zip"))
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
#' @importFrom ChemmineR read.SDFset
#' @export
read_sdf_data_lipidmaps <-
  function(file,
           path = ".") {
    message("Reading data, it may take a while...")
    lipidmaps <-
      ChemmineR::read.SDFset(sdfstr = file, skipErrors = TRUE)
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
  function(url = "https://www.lipidmaps.org/databases/lmsd",
           lipid_id = "LMFA01030001") {
    result <-
      tryCatch(rvest::read_html(paste0(url, "/", lipid_id)), error = function(e){
        NULL
      })

    if(is.null(result)){
      message("Check this url: ", paste0(url, "/", lipid_id))
      return(NULL)
    }

    ###base information
    info <-
      result %>%
      rvest::html_node(".block .px-5") %>%
      html_text2() %>%
      stringr::str_split("\n") %>%
      `[[`(1)

    idx <- grep("Show lipids differing", info)
    info <- info[1:idx - 1]
    info <- info[-grep("Select m/z", info)]

    info <-
      info[!info %in% c("Calculate m/z")]

    idx1 <-
      which(info == "Synonyms")
    idx2 <-
      which(info == "LM ID")

    if (idx2 - idx1 == 1) {
      Synonyms <- NA
    } else{
      Synonyms <-
        paste(info[(idx1 + 1):(idx2 - 1)], collapse = "{}")
    }

    Common_Name <-
      info[which(info == "Common Name") + 1]
    Common_Name <-
      ifelse(length(Common_Name) == 0, NA, Common_Name)

    Systematic_Name <-
      info[which(info == "Systematic Name") + 1]

    Systematic_Name <-
      ifelse(length(Systematic_Name) == 0, NA, Systematic_Name)

    LM_ID <-
      info[which(info == "LM ID") + 1]
    LM_ID <-
      ifelse(length(LM_ID) == 0, NA, LM_ID)

    Status <-
      info[which(info == "Status") + 1]
    Status <-
      ifelse(length(Status) == 0, NA, Status)

    Exact_Mass <-
      info[which(info == "Exact Mass") + 1] %>%
      as.numeric()
    Exact_Mass <-
      ifelse(length(Exact_Mass) == 0, NA, Exact_Mass)

    Formula <-
      info[which(info == "Formula") + 1]
    Formula <-
      ifelse(length(Formula) == 0, NA, Formula)

    Abbrev <-
      info[which(info == "Abbrev") + 1]
    Abbrev <-
      ifelse(length(Abbrev) == 0, NA, Abbrev)

    ##main
    main <-
      result %>%
      html_node("#main+ .card .px-5") %>%
      html_nodes("div") %>%
      html_text2()

    main <-
      main[-1]

    class <-
      matrix(main[1:6], nrow = 2, byrow = TRUE) %>%
      as.data.frame()

    colnames(class) <- as.character(class[1, ])
    class <- class[-1, , drop = FALSE]

    InChiKey <-
      main[which(main == "InChiKey (Click to copy)") + 1]

    InChi <-
      main[which(main == "InChi (Click to copy)") + 1]

    SMILES <-
      main[which(main == "SMILES (Click to copy)") + 1]

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
        as.character(taxonomy_information[1,])

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


    other_databases <-
      result %>%
      html_node(".mb-2 .gap-0") %>%
      html_nodes("div") %>%
      html_text2()

    other_databases <-
      other_databases %>%
      matrix(nrow = 2, byrow = FALSE) %>%
      as.data.frame()

    colnames(other_databases) <-
      as.character(other_databases[1,])

    other_databases <-
      other_databases[-1, , drop = FALSE]

    Admin <-
      result %>%
      html_node(".mb-2:nth-child(12) .px-5") %>%
      html_text2() %>%
      stringr::str_split("\n") %>%
      `[[`(1)

    Admin <-
      Admin %>%
      matrix(nrow = 2, byrow = FALSE) %>%
      as.data.frame()

    colnames(Admin) <-
      as.character(Admin[1,])

    Admin <-
      Admin[-1, , drop = FALSE]

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
        class,
        InChiKey,
        InChi,
        SMILES,
        taxonomy_information,
        other_databases,
        Admin
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
