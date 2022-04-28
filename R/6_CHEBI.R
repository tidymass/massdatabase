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
    path <- file.path(path, "CHEBI_compound")
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
    path <- file.path(path, "CHEBI_compound")

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
