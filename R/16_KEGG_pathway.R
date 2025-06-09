#' @title Request all pathway information in KEGG
#' @description Request all pathway information in KEGG
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @return A data frame.
#' @param organism organism. see https://www.genome.jp/kegg/pathway.html
#' @importFrom KEGGREST keggGet keggList
#' @importFrom magrittr %>%
#' @export
#' @examples
#' head(request_kegg_pathway_info, 3)

request_kegg_pathway_info <-
  function(organism = "hsa") {
    pathway_id <-
      KEGGREST::keggList(database = "pathway", organism = organism)
    result <-
      data.frame(KEGG.ID = names(pathway_id),
                 Pathway.name = unname(pathway_id)) %>%
      dplyr::mutate(KEGG.ID = stringr::str_replace(KEGG.ID, "path\\:", ""))

    return(result)
  }




#' @title Request one specific pathway information in KEGG
#' @description Request one specific pathway information in KEGG
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param pathway_id pathway id. For example, hsa00010
#' @importFrom KEGGREST keggGet
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_kegg_pathway(pathway_id = "hsa00010")
#' x[1:2]
#'
request_kegg_pathway <-
  function(pathway_id = "hsa00010") {
    x <-
      KEGGREST::keggGet(dbentries = pathway_id)[[1]]

    return(x)
  }



#' @title Download KEGG pathway data
#' @description Download KEGG pathway data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .
#' @param sleep Default is 1 second.
#' @param organism organism. see https://www.genome.jp/kegg/pathway.html
#' @return KEGG pathway database, rda format.
#' @export

download_kegg_pathway <-
  function(path = ".",
           sleep = 1,
           organism = "hsa") {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    kegg_id <-
      request_kegg_pathway_info(organism = organism)

    pb <- progress::progress_bar$new(total = nrow(kegg_id))

    kegg_pathway_database <-
      seq_along(kegg_id$KEGG.ID) %>%
      purrr::map(function(i) {
        pb$tick()
        Sys.sleep(time = sleep)
        KEGGREST::keggGet(dbentries = kegg_id$KEGG.ID[i])[[1]]
      })

    save(kegg_pathway_database,
         file = file.path(path, "kegg_pathway_database"))
  }



#' @title Read the KEGG pathway database from download_kegg_pathway function
#' @description Read the KEGG pathway database from download_kegg_pathway function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with download_kegg_pathway function.
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @export
read_kegg_pathway <-
  function(path = ".") {
    load(file.path(path, "kegg_pathway_database"))

    pb <-
      progress::progress_bar$new(total = length(kegg_pathway_database))

    kegg_data <-
      seq_len(length(kegg_pathway_database)) %>%
      purrr::map(function(i) {
        pb$tick()
        pathway_id <-
          kegg_pathway_database[[i]]$ENTRY

        pathway_name <-
          kegg_pathway_database[[i]]$PATHWAY_MAP %>%
          unname()

        describtion <-
          kegg_pathway_database[[i]]$DESCRIPTION

        pathway_class <-
          kegg_pathway_database[[i]]$CLASS

        gene_list <-
          kegg_pathway_database[[i]]$GENE

        if (is.null(gene_list)) {
          gene_list <- data.frame()
        } else{
          gene_list <-
            data.frame(
              KEGG.ID = gene_list[seq(1, length(gene_list) - 1, by = 2)],
              Gene.name = gene_list[seq(2, length(gene_list), by = 2)],
              stringsAsFactors = FALSE
            )
        }

        compound_list <-
          kegg_pathway_database[[i]]$COMPOUND

        compound_list <-
          data.frame(
            KEGG.ID = names(kegg_pathway_database[[i]]$COMPOUND),
            Compound.name = kegg_pathway_database[[i]]$COMPOUND,
            stringsAsFactors = FALSE
          )

        related_disease <-
          data.frame(
            Disease.ID = names(kegg_pathway_database[[i]]$DISEASE),
            Disease.name =  kegg_pathway_database[[i]]$DISEASE,
            stringsAsFactors = FALSE
          )

        related_module <-
          data.frame(
            Module.ID = names(kegg_pathway_database[[i]]$MODULE),
            Module.name = kegg_pathway_database[[i]]$MODULE,
            stringsAsFactors = FALSE
          )

        list(
          pathway_id = pathway_id,
          pathway_name = pathway_name,
          describtion = describtion,
          pathway_class = pathway_class,
          gene_list = gene_list,
          compound_list = compound_list,
          related_disease = related_disease,
          related_module = related_module
        )
      })
    return(kegg_data)
  }









#' @title Convert KEGG pathway data (list,
#' from read_kegg_pathway function)
#' to metpath format database
#' @description Convert KEGG pathway data (list,
#' from read_kegg_pathway function)
#' to metpath format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_kegg_pathway function
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

    describtion <-
      lapply(data, function(x)
      paste0(x$describtion, collapse = "{}"))

    describtion[describtion == ""] <- NA

    pathway_class <-
      lapply(data, function(x)
        paste0(x$pathway_class, collapse = "{}"))

    pathway_class[pathway_class == ""] <- NA

    pathway <-
      new(
        Class = "pathway_database",
        database_info = list(source = "KEGG",
                             version = as.character(Sys.Date())),
        pathway_id = lapply(data, function(x)
          x$pathway_id) %>% unlist(),
        pathway_name = lapply(data, function(x)
          x$pathway_name) %>% unlist(),
        describtion = describtion,
        pathway_class = pathway_class,
        gene_list = lapply(data, function(x)
          x$gene_list),
        compound_list = lapply(data, function(x)
          x$compound_list),
        protein_list = list(),
        reference_list = list(),
        related_disease = lapply(data, function(x)
          x$related_disease),
        related_module = lapply(data, function(x)
          x$related_module)
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
    kegg_pathway <- pathway
    rm(list = c("pathway"))
    save(kegg_pathway, file = file.path(path, "kegg_pathway"))
    invisible(kegg_pathway)
  }
