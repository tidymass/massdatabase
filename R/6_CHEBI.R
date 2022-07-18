#' @title Download the CHEBI compound database
#' @description Download the CHEBI compound database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/".
#' @param path Default is ..
#' @return Downloaded files.
#' @importFrom magrittr %>%
#' @export
download_chebi_compound <-
  function(url = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/",
           path = ".") {
    path <- file.path(path, "data")
    dir.create(path)
    message("Download chebiId_inchi.tsv...\n")
    download.file(
      url = paste0(url, "chebiId_inchi.tsv"),
      destfile = file.path(path, "chebiId_inchi.tsv")
    )
    message("Done.\n")

    message("Download chemical_data.tsv...\n")
    download.file(
      url = paste0(url, "chemical_data.tsv"),
      destfile = file.path(path, "chemical_data.tsv")
    )
    message("Done.\n")

    message("Download comments.tsv...\n")
    download.file(url = paste0(url, "comments.tsv"),
                  destfile = file.path(path, "comments.tsv"))
    message("Done.\n")

    message("Download compound_origins.tsv...\n")
    download.file(
      url = paste0(url, "compound_origins.tsv"),
      destfile = file.path(path, "compound_origins.tsv")
    )
    message("Done.\n")

    message("Download compounds.tsv.gz...\n")
    download.file(
      url = paste0(url, "compounds.tsv.gz"),
      destfile = file.path(path, "compounds.tsv.gz")
    )
    message("Done.\n")

    message("Download database_accession.tsv...\n")
    download.file(
      url = paste0(url, "database_accession.tsv"),
      destfile = file.path(path, "database_accession.tsv")
    )
    message("Done.\n")

    message("Download names.tsv.gz...\n")
    download.file(url = paste0(url, "names.tsv.gz"),
                  destfile = file.path(path, "names.tsv.gz"))
    message("Done.\n")

    message("Download names.tsv.gz...\n")
    download.file(url = paste0(url, "names.tsv.gz"),
                  destfile = file.path(path, "names.tsv.gz"))
    message("Done.\n")
  }



#' @title Read the CHEBI compound database from download_chebi_compound function
#' @description Read the CHEBI compound database from download_chebi_compound function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with download_chebi_compound function.
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @export
read_chebi_compound <-
  function(path = ".") {
    path <- file.path(path, "data")

    chebiid_inchi <-
      readr::read_delim(file.path(path, "chebiId_inchi.tsv")) %>%
      dplyr::mutate(CHEBI_ID = paste0("CHEBI:", CHEBI_ID))

    chemical_data <-
      readr::read_delim(file.path(path, "chemical_data.tsv")) %>%
      dplyr::mutate(COMPOUND_ID = paste0("CHEBI:", COMPOUND_ID))

    comments <-
      readr::read_delim(file.path(path, "comments.tsv")) %>%
      dplyr::mutate(COMPOUND_ID = paste0("CHEBI:", COMPOUND_ID))

    compound_origins <-
      readr::read_delim(file.path(path, "compound_origins.tsv")) %>%
      dplyr::mutate(COMPOUND_ID = paste0("CHEBI:", COMPOUND_ID))

    compound_origins <-
      compound_origins %>%
      plyr::dlply(.variables = .(COMPOUND_ID)) %>%
      purrr::map(function(x) {
        x$SPECIES <- paste(x$SPECIES_TEXT, collapse = "{}")
        x$SPECIES_ACCESSION <-
          paste(x$SPECIES_ACCESSION, collapse = "{}")
        x %>%
          dplyr::select(COMPOUND_ID, SPECIES, SPECIES_ACCESSION) %>%
          dplyr::distinct(COMPOUND_ID, .keep_all = TRUE)
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    compounds <-
      readr::read_delim(gzfile(file.path(path, "compounds.tsv.gz")))

    database_accession <-
      readr::read_delim(file.path(path, "database_accession.tsv")) %>%
      dplyr::mutate(COMPOUND_ID = paste0("CHEBI:", COMPOUND_ID))

    database_accession <-
      database_accession %>%
      dplyr::select(-c(ID, SOURCE)) %>%
      tidyr::pivot_wider(
        names_from = TYPE,
        values_from = ACCESSION_NUMBER,
        values_fn = function(x) {
          paste(unique(x), collapse = "{}")
        }
      ) %>%
      dplyr::rename(
        KEGG.ID = "KEGG COMPOUND accession",
        CAS.ID = "CAS Registry Number",
        FOODB.ID = "FooDB accession",
        METACYC.ID = "MetaCyc accession",
        HMDB.ID = "HMDB accession",
        LIPIDMAPS.ID = "LIPID MAPS instance accession",
        KEGG_DRUG.ID = "KEGG DRUG accession",
        KEGG_GLYCAN.ID = "KEGG GLYCAN accession",
        WIKIPEDIA.ID = "Wikipedia accession",
        DRUGBANK.ID = "DrugBank accession",
        PUBCHEM.ID = "Pubchem accession"
      ) %>%
      dplyr::select(
        c(
          COMPOUND_ID,
          KEGG.ID,
          CAS.ID,
          FOODB.ID,
          METACYC.ID,
          HMDB.ID,
          LIPIDMAPS.ID,
          KEGG_DRUG.ID,
          KEGG_GLYCAN.ID,
          WIKIPEDIA.ID,
          DRUGBANK.ID,
          PUBCHEM.ID
        )
      )

    names <-
      readr::read_delim(gzfile(file.path(path, "names.tsv.gz"))) %>%
      dplyr::mutate(COMPOUND_ID = paste0("CHEBI:", COMPOUND_ID))

    names <-
      names %>%
      dplyr::filter(LANGUAGE == "en") %>%
      dplyr::select(COMPOUND_ID, NAME) %>%
      plyr::dlply(.variables = .(COMPOUND_ID)) %>%
      purrr::map(function(x) {
        if (nrow(x) == 1) {
          return(x)
        } else{
          x$NAME <-
            paste(unique(x$NAME), collapse = "{}")
          x <-
            x %>%
            dplyr::distinct(COMPOUND_ID, .keep_all = TRUE)
          return(x)
        }
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    FORMULA <-
      chemical_data %>%
      dplyr::filter(TYPE == "FORMULA") %>%
      dplyr::select(-c(ID, TYPE, SOURCE)) %>%
      dplyr::rename(FORMULA = CHEMICAL_DATA) %>%
      dplyr::distinct(COMPOUND_ID, .keep_all = TRUE)

    MASS <-
      chemical_data %>%
      dplyr::filter(TYPE == "MASS") %>%
      dplyr::select(-c(ID, TYPE, SOURCE)) %>%
      dplyr::rename(MASS = CHEMICAL_DATA) %>%
      dplyr::distinct(COMPOUND_ID, .keep_all = TRUE)

    MONOISOTOPIC_MASS <-
      chemical_data %>%
      dplyr::filter(TYPE == "MONOISOTOPIC MASS") %>%
      dplyr::select(-c(ID, TYPE, SOURCE)) %>%
      dplyr::rename(MONOISOTOPIC_MASS = CHEMICAL_DATA) %>%
      dplyr::distinct(COMPOUND_ID, .keep_all = TRUE)

    CHARGE <-
      chemical_data %>%
      dplyr::filter(TYPE == "CHARGE") %>%
      dplyr::select(-c(ID, TYPE, SOURCE)) %>%
      dplyr::rename(CHARGE = CHEMICAL_DATA) %>%
      dplyr::distinct(COMPOUND_ID, .keep_all = TRUE)

    chebi_data <-
      FORMULA %>%
      dplyr::full_join(MASS, by = "COMPOUND_ID") %>%
      dplyr::full_join(MONOISOTOPIC_MASS, by = "COMPOUND_ID") %>%
      dplyr::full_join(CHARGE, by = "COMPOUND_ID") %>%
      dplyr::full_join(chebiid_inchi, by = c("COMPOUND_ID" = "CHEBI_ID")) %>%
      dplyr::full_join(compound_origins, by = c("COMPOUND_ID")) %>%
      dplyr::full_join(compounds, by = c("COMPOUND_ID" = "CHEBI_ACCESSION")) %>%
      dplyr::full_join(database_accession, by = c("COMPOUND_ID")) %>%
      dplyr::full_join(names %>% dplyr::rename(Synonyms = NAME),
                       by = c("COMPOUND_ID"))

    invisible(chebi_data)
  }




#' @title Request CHEBI compound
#' @description Request CHEBI compound
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url url https://www.ebi.ac.uk/chebi/searchId.do?chebiId=
#' @param compound_id compound_id, for example, CHEBI:18358
#' @return Compound information (list)
#' @importFrom stringr str_replace_all str_extract str_replace str_split
#' @importFrom stringr str_trim
#' @importFrom dplyr select filter bind_rows
#' @importFrom tibble as_tibble
#' @importFrom rvest html_elements html_table
#' @importFrom purrr map
#' @export
#' @examples
#' request_chebi_compound(compound_id = "CHEBI:18358")
#' request_chebi_compound(compound_id = "CHEBI:18377")
#' request_chebi_compound(compound_id = "CHEBI:18399")
request_chebi_compound <-
  function(url = "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=",
           compound_id = "CHEBI:18358") {

    result <-
      tryCatch(
        rvest::read_html(x = paste0(url, compound_id)),
        error = function(e)
          NULL
      )

    if (is.null(result)) {
      message("Check url: ", paste0(url, compound_id))
      return(NULL)
    }

    main <-
      tryCatch(
        result %>%
          rvest::html_elements(".gridLayoutCellStructure+ td") %>%
          rvest::html_table() %>%
          `[[`(1) %>%
          dplyr::select(X1, X2) %>%
          dplyr::filter(!X1 %in% c("", "Supplier Information", "Download")) %>%
          t() %>%
          tibble::as_tibble(),
        error = function(e)
          NULL
      )

    if (is.null(main)) {
      message("Can't find ", compound_id, ", check it.")
      return(NULL)
    }

    colnames(main) <- as.character(main[1,])
    main <- main[-1, , drop = FALSE]

    base_info <-
      tryCatch(
        result %>%
          rvest::html_elements(".small-boxed-section+ .small-boxed-section") %>%
          rvest::html_table() %>%
          `[[`(1) %>%
          dplyr::select(X1, X2),
        error = function(e)
          NULL
      )

    if (is.null(base_info)) {
      base_info <-
        tryCatch(
          result %>%
            rvest::html_elements(".chebiTableContent+ .small-boxed-section") %>%
            rvest::html_table() %>%
            `[[`(1) %>%
            dplyr::select(X1, X2),
          error = function(e)
            NULL
        )
    }

    base_info <-
      base_info$X1[1] %>%
      stringr::str_split("\n") %>%
      `[[`(1) %>%
      stringr::str_trim()

    base_info <- base_info[base_info != ""]

    base_info <-
      matrix(base_info, nrow = 2) %>%
      tibble::as_tibble()

    colnames(base_info) <- as.character(base_info[1,])

    base_info <- base_info[-1, , drop = FALSE]

    others <-
      result %>%
      rvest::html_elements(css = ".chebiTableContent")

    others <-
      seq_along(others) %>%
      purrr::map(function(i) {
        # cat(i, " ")
        x <- others[[i]]
        x <-
          tryCatch(
            rvest::html_table(x),
            error = function(e)
              NULL
          )

        if (!is.null(x)) {
          if (nrow(x) > 0) {
            x <-
              x %>%
              dplyr::filter(X1 != "")
          } else{
            return(NULL)
          }
        }
        x
      })

    remain_idx <-
      lapply(others, is.null) %>%
      unlist() %>%
      `!`() %>%
      which()

    others <- others[remain_idx]

    ###source
    source <-
      lapply(others, function(x) {
        if (any(x$X1 == "Metabolite of Species")) {
          return(x)
        } else{
          return(NULL)
        }
      }) %>%
      dplyr::bind_rows()

    if (nrow(source) == 0) {
      source <- NULL
    } else{
      source$X1[2] <-
        source$X1[2] %>%
        stringr::str_split("\n") %>%
        `[[`(1) %>%
        stringr::str_trim() %>%
        paste(collapse = "")

      source$X2[2] <-
        source$X2[2] %>%
        stringr::str_split("\n") %>%
        `[[`(1) %>%
        stringr::str_trim() %>%
        paste(collapse = "")

      colnames(source) <- as.character(source[1,])
      source <- source[-1, , drop = FALSE]

    }

    ###Roles Classification
    roles_classification <-
      lapply(others, function(x) {
        if (any(x$X1 == "Roles Classification")) {
          return(x)
        } else{
          return(NULL)
        }
      }) %>%
      dplyr::bind_rows()

    roles_classification <-
      roles_classification[-1, , drop = FALSE]

    roles_classification$X2 <-
      roles_classification$X2 %>%
      stringr::str_split("\n") %>%
      lapply(function(x) {
        x <-
          x %>%
          stringr::str_trim()
        x <- x[x != ""]
        x %>%
          paste(collapse = " ")
      }) %>%
      unlist()

    roles_classification <-
      roles_classification %>%
      t() %>%
      tibble::as_tibble()

    colnames(roles_classification) <-
      as.character(roles_classification[1,])

    roles_classification <- roles_classification[-1, , drop = FALSE]

    ###IUPAC Name
    iupac_name <-
      lapply(others, function(x) {
        if (any(x$X1 == "IUPAC Name")) {
          return(x)
        } else{
          return(NULL)
        }
      }) %>%
      dplyr::bind_rows()

    colnames(iupac_name) <- as.character(iupac_name[1,])
    iupac_name <- iupac_name[-1, , drop = FALSE]

    ###Synonyms
    synonyms <-
      lapply(others, function(x) {
        if (any(x$X1 == "Synonyms")) {
          return(x)
        } else{
          return(NULL)
        }
      }) %>%
      dplyr::bind_rows()

    colnames(synonyms) <- as.character(synonyms[1,])
    synonyms <- synonyms[-1, , drop = FALSE]

    ###Manual Xrefs
    manual_xrefs <-
      lapply(others, function(x) {
        if (any(x$X1 == "Manual Xrefs")) {
          return(x)
        } else{
          return(NULL)
        }
      }) %>%
      dplyr::bind_rows()

    if (nrow(manual_xrefs) == 0) {
      manual_xrefs <- NULL
    } else{
      manual_xrefs <-
        manual_xrefs %>%
        dplyr::filter(!is.na(X2))

      colnames(manual_xrefs) <- as.character(manual_xrefs[1,])
      manual_xrefs <- manual_xrefs[-1, , drop = FALSE]

    }

    ###Registry Numbers
    registry_numbers <-
      lapply(others, function(x) {
        if (any(x$X1 == "Registry Numbers")) {
          return(x)
        } else{
          return(NULL)
        }
      }) %>%
      dplyr::bind_rows()

    if (nrow(registry_numbers) == 0) {
      registry_numbers <- NULL
    } else{
      registry_numbers <-
        registry_numbers %>%
        dplyr::filter(!is.na(X2))

      colnames(registry_numbers) <-
        as.character(registry_numbers[1,])
      registry_numbers <- registry_numbers[-1, , drop = FALSE]
    }

    list(
      main = main,
      base_info = base_info,
      source = source,
      roles_classification = roles_classification,
      iupac_name = iupac_name,
      synonyms = synonyms,
      manual_xrefs = manual_xrefs,
      registry_numbers = registry_numbers
    )
  }






#' @title Convert BIGG universal metabolites (data.frame,
#' from read_chebi_metabolite)
#' to metID format database
#' @description Convert BIGG universal metabolites (data.frame,
#' from read_chebi_metabolite)
#' to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data data.frame, from read_chebi_metabolite function
#' @param path Default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @importFrom dplyr one_of
#' @export

convert_chebi2metid <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    data <-
      data %>%
      dplyr::rename(
        Lab.ID = COMPOUND_ID,
        Formula = FORMULA,
        mz = MONOISOTOPIC_MASS,
        Charge = CHARGE,
        INCHI.ID = InChI,
        Species = SPECIES,
        Species.ID = SPECIES_ACCESSION,
        Status = STATUS,
        Database_source = SOURCE,
        Compound.name = NAME,
        Definition = DEFINITION,
        Star = STAR,
        Updated_date = MODIFIED_ON
      ) %>%
      dplyr::mutate(
        CHEBI.ID = Lab.ID,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "CHEBI"
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
      ) %>%
      as.data.frame()

    data[which(data == "null", arr.ind = TRUE)] <-
      NA

    data$mz <- as.numeric(data$mz)

    ###remove mz is NA
    data <-
      data %>%
      dplyr::filter(!is.na(mz))

    species <-
      sort(unique(unlist(
        stringr::str_split(data$Species[!is.na(data$Species)], "\\{\\}")
      )))

    species.id <-
      sort(unique(unlist(
        stringr::str_split(data$Species.ID[!is.na(data$Species.ID)], "\\{\\}")
      )))

    # library("taxize")
    # gnr_resolve(species[1])

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = data,
                     file = file.path(temp_file, "data.csv"))

    chebi_ms1 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "data.csv",
        source = "CHEBI",
        link = "https://www.ebi.ac.uk/chebi/init.do",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "data.csv"))
    unlink(temp_file)

    save(chebi_ms1, file = file.path(path, "chebi_ms1"))
    invisible(chebi_ms1)
  }
