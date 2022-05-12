#' @title Download SMPDB pathway data
#' @description Download SMPDB pathway data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .
#' @return SMPDB pathway database, csv format.
#' @importFrom utils unzip
#' @export

download_smpdb_pathway <-
  function(path = ".") {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    message("Downloading pathway information...")
    utils::download.file(
      "https://www.smpdb.ca/downloads/smpdb_pathways.csv.zip",
      destfile = file.path(path, "smpdb_pathways.csv.zip")
    )

    message("Downloading metabolite information...")
    utils::download.file(
      "https://www.smpdb.ca/downloads/smpdb_metabolites.csv.zip",
      destfile = file.path(path, "smpdb_metabolites.csv.zip")
    )

    message("Downloading protein information...")
    utils::download.file(
      "https://www.smpdb.ca/downloads/smpdb_proteins.csv.zip",
      destfile = file.path(path, "smpdb_proteins.csv.zip")
    )

    utils::unzip(
      zipfile = file.path(path, "smpdb_pathways.csv.zip"),
      exdir = file.path(path, "smpdb_pathways"),
      overwrite = TRUE
    )
    utils::unzip(
      zipfile = file.path(path, "smpdb_metabolites.csv.zip"),
      exdir = file.path(path, "smpdb_metabolites"),
      overwrite = TRUE
    )
    utils::unzip(
      zipfile = file.path(path, "smpdb_proteins.csv.zip"),
      exdir = file.path(path, "smpdb_proteins"),
      overwrite = TRUE
    )

    unlink(file.path(path, "smpdb_pathways.csv.zip"))
    unlink(file.path(path, "smpdb_metabolites.csv.zip"))
    unlink(file.path(path, "smpdb_proteins.csv.zip"))
  }



#' @title Read the SMPDB pathway database from download_smpdb_pathway function
#' @description Read the SMPDB pathway database from download_smpdb_pathway function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with download_smpdb_pathway function.
#' @param only_primarity_pathway Only remain primary pathway?
#' @param remain_disease_pathway remain disease pathway?
#' @param remain_drug_pathway remain drug pathway?
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @importFrom progress progress_bar
#' @importFrom utils data
#' @export
read_smpdb_pathway <-
  function(path = ".",
           only_primarity_pathway = TRUE,
           remain_disease_pathway = TRUE,
           remain_drug_pathway = TRUE) {

    smpdb_pahtway <-
      readr::read_csv(file.path(path, "smpdb_pathways/smpdb_pathways.csv"),
                      show_col_types = FALSE)

    data("smpdb_primary_pathway_id", envir = environment())

    smpdb_pahtway$Subject[match(smpdb_primary_pathway_id, smpdb_pahtway$`SMPDB ID`)] <-
      paste(smpdb_pahtway$Subject[match(smpdb_primary_pathway_id, smpdb_pahtway$`SMPDB ID`)],
            "primary_pathway", sep = "{}")

    if (only_primarity_pathway) {
      smpdb_pahtway <-
        smpdb_pahtway %>%
        dplyr::filter(stringr::str_detect(Subject, "primary_pathway")) %>%
        dplyr::arrange(`SMPDB ID`)
    }

    if (!remain_disease_pathway) {
      smpdb_pahtway <-
        smpdb_pahtway %>%
        dplyr::filter(!stringr::str_detect(Subject, "Disease")) %>%
        dplyr::arrange(`SMPDB ID`)
    }

    if (!remain_drug_pathway) {
      smpdb_pahtway <-
        smpdb_pahtway %>%
        dplyr::filter(!stringr::str_detect(Subject, "Drug")) %>%
        dplyr::arrange(`SMPDB ID`)
    }

    message("Reading metabolites data...")
    file <- dir(file.path(path, "smpdb_metabolites"))

    file <-
      file[file %in% paste0(smpdb_pahtway$`SMPDB ID`, "_metabolites.csv")] %>%
      sort()

    pb <-
      progress::progress_bar$new(total = length(file))

    compound_list <-
      seq_len(length(file)) %>%
      purrr::map(function(i) {
        pb$tick()
        result <-
        readr::read_csv(file.path(path, "smpdb_metabolites", file[i]),
                        show_col_types = FALSE) %>%
          dplyr::select("HMDB ID", "KEGG ID", "Metabolite Name") %>%
          dplyr::rename(HMDB.ID = "HMDB ID",
                        KEGG.ID = "KEGG ID",
                        Compound.name = "Metabolite Name") %>%
          as.data.frame()
        result
      })

    names(compound_list) <-
      stringr::str_replace(file, "_metabolites.csv", "")

    message("Reading proteins data...")

    file <- dir(file.path(path, "smpdb_proteins"))

    file <-
      file[file %in% paste0(smpdb_pahtway$`SMPDB ID`, "_proteins.csv")] %>%
      sort()

    pb <-
      progress::progress_bar$new(total = length(file))

    protein_list <-
      seq_len(length(file)) %>%
      purrr::map(function(i) {
        pb$tick()
        result <-
          readr::read_csv(file.path(path, "smpdb_proteins", file[i]),
                          show_col_types = FALSE) %>%
          dplyr::select("Uniprot ID", "HMDBP ID", "GenBank ID", "Gene Name", "Protein Name") %>%
          dplyr::rename(Uniprot.ID = "Uniprot ID",
                        HMDBP.ID = "HMDBP ID",
                        GenBank.ID = "GenBank ID",
                        Gene.name = "Gene Name",
                        Protein.name = "Protein Name") %>%
          as.data.frame()
        result
      })

    names(protein_list) <-
      stringr::str_replace(file, "_proteins.csv", "")

    pathway_id <-
      intersect(names(compound_list),
                names(protein_list))

    pathway_name <-
      smpdb_pahtway$Name[match(pathway_id, smpdb_pahtway$`SMPDB ID`)]

    describtion <-
      smpdb_pahtway$Description[match(pathway_id, smpdb_pahtway$`SMPDB ID`)] %>%
      purrr::map(function(x){
        x
      })

    pathway_class <-
      smpdb_pahtway$Subject[match(pathway_id, smpdb_pahtway$`SMPDB ID`)] %>%
      purrr::map(function(x){
        x
      })

    compound_list <- compound_list[pathway_id]

    protein_list <- protein_list[pathway_id]

    return(list(
      pathway_id = pathway_id,
      pathway_name = pathway_name,
      describtion = describtion,
      pathway_class = pathway_class,
      compound_list = compound_list,
      protein_list = protein_list
    ))
  }





#' @title Convert smpdb pathway data (list,
#' from read_smpdb_pathway function)
#' to metpath format database
#' @description Convert smpdb pathway data (list,
#' from read_smpdb_pathway function)
#' to metpath format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_smpdb_pathway function
#' @param path Default is .
#' @param threads threads
#' @return metpath pathway_database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importClassesFrom metpath pathway_database
#' @export

convert_kegg2metpath <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    pathway <-
      new(
        Class = "pathway_database",
        database_info = list(source = "SMPDB",
                             version = as.character(Sys.Date())),
        pathway_id = data$pathway_id,
        pathway_name = data$pathway_name,
        describtion = data$describtion,
        pathway_class = data$pathway_class,
        gene_list = list(),
        compound_list = data$compound_list,
        protein_list = data$protein_list,
        reference_list = list(),
        related_disease = list(),
        related_module = list()
      )

    if (length(pathway@gene_list) == 0) {
      pathway@gene_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }

    if (length(pathway@compound_list) == 0) {
      pathway@compound_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }

    if (length(pathway@protein_list) == 0) {
      pathway@protein_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }
    smpdb_pathway <- pathway
    rm(list = c("pathway"))
    save(smpdb_pathway, file = file.path(path, "smpdb_pathway"))
    invisible(smpdb_pathway)
  }
