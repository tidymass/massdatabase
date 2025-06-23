#' Retrieve Available Organisms from WikiPathways
#'
#' This function queries WikiPathways to obtain a list of organisms for which pathways are available.
#'
#' @return A character vector of organism names available in the WikiPathways database.
#' @export
#' @examples
#' request_wikipathway_organisms_info()

request_wikipathway_organisms_info <-
  function() {
    result <- rjson::fromJSON(file = "https://www.wikipathways.org/json/listOrganisms.json")
    return(result$organisms)
  }

#' Retrieve Pathway Information for a Specific Organism
#'
#' This function fetches all available pathways from WikiPathways for a specified organism.
#'
#' @param organism A character string specifying the organism name. Default is "Homo sapiens".
#' @return A data frame containing pathway information, including pathway names, IDs, and species.
#' @export
#' @examples
#' request_wikipathway_info("Homo sapiens")

request_wikipathway_info <-
  function(organism = "Homo sapiens") {
    organism <-
      match.arg(organism, choices = wikipathway_organisms_info)
    result <-
      rjson::fromJSON(file = "https://www.wikipathways.org/json/listPathways.json")

    result <- result$organisms %>%
      purrr::map_dfr(~ .x$pathways)

    if (organism != "") {
      result <- dplyr::filter(result, species == organism)
    }
    if (nrow(result) == 0) {
      message("No resultults")
      return(NULL)
    }
    return(result)
  }


##' Retrieve Detailed Pathway Data from WikiPathways
##'
##' This function fetches detailed information about a specific pathway from WikiPathways using its pathway ID.
##'
##' @param pathway_id A character string specifying the WikiPathways ID (e.g., "WP5293").
##' @return A list containing:
##'   \itemize{
##'     \item \code{pathway_info} - Metadata about the pathway.
##'     \item \code{pathway_description} - A description of the pathway.
##'     \item \code{node_info} - A data frame containing gene/metabolite nodes.
##'     \item \code{edge_info} - A data frame containing interaction edges.
##'   }
##' @export
##' @examples
##' data = request_wikipathway("WP5293")
##' data$pathway_info
#
# request_wikipathway <-
#   function(pathway_id = "WP5293") {
#     url <-
#       paste0(
#         "https://www.wikipathways.org/wikipathways-assets/pathways/",
#         pathway_id,
#         "/",
#         pathway_id,
#         ".gpml"
#       )
#
#     # Read the XML file
#     result <-
#       tryCatch(
#         expr = xml2::read_xml(url),
#         error = function(e) {
#           message("Error: ", e)
#           return(NULL)
#         }
#       )
#
#     # Extract pathway information
#     pathway_name <- xml2::xml_attr(result, "Name")
#     data_source <- xml2::xml_attr(result, "Data-Source")
#     version <- xml2::xml_attr(result, "Version")
#     author <- xml2::xml_attr(result, "Author")
#     last_modified <- xml2::xml_attr(result, "Last-Modified")
#     organism <- xml2::xml_attr(result, "Organism")
#
#     # Print pathway information
#     pathway_info <- list(
#       Name = pathway_name,
#       DataSource = data_source,
#       Version = version,
#       Author = author,
#       LastModified = last_modified,
#       Organism = organism
#     )
#
#     # Extract pathway description
#     pathway_description <-
#       xml2::xml_find_all(result, ".//d1:Comment", ns = xml2::xml_ns(result))
#
#     # Get the text content of each comment
#     pathway_description <-
#       xml2::xml_text(pathway_description)
#
#     if (length(pathway_description) > 1) {
#       pathway_description <-
#         pathway_description[1]
#     }
#
#     # Extract BiopaxRef elements
#     biopax_refs <-
#       xml2::xml_find_all(result, ".//d1:BiopaxRef", ns = xml2::xml_ns(result))
#
#     # Get the text content of each BiopaxRef
#     biopax_ref_texts <-
#       xml2::xml_text(biopax_refs)
#
#     # Extract all DataNode elements
#     datanodes <-
#       xml2::xml_find_all(result, ".//d1:DataNode", ns = xml2::xml_ns(result))
#
#     # Extract Xref elements within DataNodes
#     xref_nodes <- xml2::xml_find_all(datanodes, ".//d1:Xref", ns = xml2::xml_ns(result))
#
#     # Extract TextLabel, GraphId, and Type attributes
#     node_info <- data.frame(
#       TextLabel = xml2::xml_attr(datanodes, "TextLabel"),
#       GraphId = xml2::xml_attr(datanodes, "GraphId"),
#       Type = xml2::xml_attr(datanodes, "Type"),
#       Database = xml2::xml_attr(xref_nodes, "Database"),
#       ID = xml2::xml_attr(xref_nodes, "ID"),
#       stringsAsFactors = FALSE
#     )
#
#
#     # # Extract all Graphics elements inside DataNodes
#     # graphics_nodes <-
#     #   xml2::xml_find_all(result, ".//d1:DataNode/d1:Graphics", ns = xml2::xml_ns(result))
#     #
#     # # Extract positions (CenterX, CenterY)
#     # graphics_info <-
#     #   data.frame(
#     #   GraphId = xml2::xml_attr(
#     #     xml2::xml_find_all(result, ".//d1:DataNode", ns = xml2::xml_ns(result)),
#     #     "GraphId"
#     #   ),
#     #   CenterX = as.numeric(xml2::xml_attr(graphics_nodes, "CenterX")),
#     #   CenterY = as.numeric(xml2::xml_attr(graphics_nodes, "CenterY")),
#     #   stringsAsFactors = FALSE
#     # )
#
#
#     # Extract all Interaction elements
#     interactions <-
#       xml2::xml_find_all(result, ".//d1:Interaction", ns = xml2::xml_ns(result))
#
#     # Extract GraphId and references to connecting nodes
#     edge_info <-
#       data.frame(
#         GraphId = xml2::xml_attr(interactions, "GraphId"),
#         Points = sapply(interactions, function(inter) {
#           points <- xml2::xml_find_all(inter, ".//d1:Point", ns = xml2::xml_ns(result))
#           paste(xml2::xml_attr(points, "GraphRef"), collapse = " -> ")
#         }),
#         stringsAsFactors = FALSE
#       )
#
#
#     return_result <-
#       list(
#         pathway_info = pathway_info,
#         pathway_description = pathway_description,
#         node_info = node_info,
#         edge_info = edge_info
#       )
#     return(return_result)
#   }

#' Retrieve Detailed Pathway Data from WikiPathways
#'
#' This function fetches detailed information about a specific pathway from WikiPathways using its pathway ID.
#'
#' @param pathway_id A character string specifying the WikiPathways ID (e.g., "WP5293").
#' @return A list containing:
#'   \itemize{
#'     \item \code{pathway_info} - Metadata about the pathway.
#'     \item \code{pathway_description} - A description of the pathway.
#'     \item \code{class_info} - A list of ontology terms associated with the "Pathway Ontology".
#'     \item \code{node_info} - A data frame containing gene/metabolite nodes.
#'     \item \code{edge_info} - A data frame containing interaction edges.
#'   }
#' @export
#' @examples
#' data = request_wikipathway("WP5293")
#' data$pathway_info

request_wikipathway <-

  function(pathway_id = "WP5293") {

    url <-
      paste0(
        "https://www.wikipathways.org/wikipathways-assets/pathways/",
        pathway_id,
        "/",
        pathway_id,
        ".gpml"
      )

    # Read the XML file
    result <-
      tryCatch(
        expr = xml2::read_xml(url),
        error = function(e) {
          message("Error: ", e)
          return(NULL)
        }
      )

    # Extract ontology info and only keep "Pathway Ontology"
    ontology_terms <- xml2::xml_find_all(result, ".//bp:openControlledVocabulary", ns = xml2::xml_ns(result))
    # ontology_info_all <- purrr::map_df(ontology_terms, function(x) {
    #   tibble::tibble(
    #     TERM = xml2::xml_text(xml2::xml_find_first(x, ".//bp:TERM", xml2::xml_ns(result))),
    #     ID = xml2::xml_text(xml2::xml_find_first(x, ".//bp:ID", xml2::xml_ns(result))),
    #     Ontology = xml2::xml_text(xml2::xml_find_first(x, ".//bp:Ontology", xml2::xml_ns(result)))
    #   )
    # })

    # Check if ontology_terms exist
    if (length(ontology_terms) > 0) {
      ontology_info_all <- purrr::map_df(ontology_terms, function(x) {
        TERM <- xml2::xml_text(xml2::xml_find_first(x, ".//bp:TERM", xml2::xml_ns(result)))
        ID <- xml2::xml_text(xml2::xml_find_first(x, ".//bp:ID", xml2::xml_ns(result)))
        Ontology <- xml2::xml_text(xml2::xml_find_first(x, ".//bp:Ontology", xml2::xml_ns(result)))

        tibble::tibble(
          TERM = ifelse(length(TERM) == 0, NA_character_, TERM),
          ID = ifelse(length(ID) == 0, NA_character_, ID),
          Ontology = ifelse(length(Ontology) == 0, NA_character_, Ontology)
        )
      })
    } else {
      ontology_info_all <- tibble::tibble(TERM = character(0), ID = character(0), Ontology = character(0))
    }


    # Only keep the "Pathway Ontology" entry
    filtered_ontology <- dplyr::filter(ontology_info_all, Ontology == "Pathway Ontology")
    # term_column <- dplyr::pull(filtered_ontology, TERM)
    # ontology_info <- term_column[1]

    # check if filtered_ontology is not empty
    # if (nrow(filtered_ontology) > 0) {
    #   term_column <- dplyr::pull(filtered_ontology, TERM)
    #   ontology_info <- term_column[1]
    # } else {
    #   ontology_info <- NA_character_
    # }

    # Get all TERM values as a list (each element is a character string)
    if (nrow(filtered_ontology) > 0) {
      ontology_info <- as.list(dplyr::pull(filtered_ontology, TERM))
    } else {
      ontology_info <- list()
    }


    # Extract pathway information
    pathway_name <- xml2::xml_attr(result, "Name")
    data_source <- xml2::xml_attr(result, "Data-Source")
    version <- xml2::xml_attr(result, "Version")
    author <- xml2::xml_attr(result, "Author")
    last_modified <- xml2::xml_attr(result, "Last-Modified")
    class <- xml2::xml_attr(result, "Communiti")
    organism <- xml2::xml_attr(result, "Organism")

    # Print pathway information
    pathway_info <- list(
      Name = pathway_name,
      DataSource = data_source,
      Version = version,
      Author = author,
      LastModified = last_modified,
      Organism = organism
    )

    # Extract pathway description
    pathway_description <-
      xml2::xml_find_all(result, ".//d1:Comment", ns = xml2::xml_ns(result))

    # Get the text content of each comment
    pathway_description <-
      xml2::xml_text(pathway_description)

    if (length(pathway_description) > 1) {
      pathway_description <-
        pathway_description[1]
    }

    # Extract BiopaxRef elements
    biopax_refs <-
      xml2::xml_find_all(result, ".//d1:BiopaxRef", ns = xml2::xml_ns(result))

    # Get the text content of each BiopaxRef
    biopax_ref_texts <-
      xml2::xml_text(biopax_refs)

    # Extract all DataNode elements
    datanodes <-
      xml2::xml_find_all(result, ".//d1:DataNode", ns = xml2::xml_ns(result))

    # Extract Xref elements within DataNodes
    xref_nodes <- xml2::xml_find_all(datanodes, ".//d1:Xref", ns = xml2::xml_ns(result))

    # Extract TextLabel, GraphId, and Type attributes
    node_info <- data.frame(
      TextLabel = xml2::xml_attr(datanodes, "TextLabel"),
      GraphId = xml2::xml_attr(datanodes, "GraphId"),
      Type = xml2::xml_attr(datanodes, "Type"),
      Database = xml2::xml_attr(xref_nodes, "Database"),
      ID = xml2::xml_attr(xref_nodes, "ID"),
      stringsAsFactors = FALSE
    )

    is_chebi    <- node_info$Database == "ChEBI"

    idx_chebi <- which( is_chebi & node_info$ID != "" & !is.na(node_info$ID) )
    if (length(idx_chebi) > 0) {
      node_info$ID[idx_chebi] <- paste0(
        "CHEBI:",
        sub("^CHEBI:", "", node_info$ID[idx_chebi], ignore.case = TRUE)
      )
    }

    is_pubchem  <- node_info$Database == "PubChem-compound"

    idx_pubchem <- which( is_pubchem & node_info$ID != "" & !is.na(node_info$ID) )
    if (length(idx_pubchem) > 0) {
      node_info$ID[idx_pubchem] <- paste0(
        "pubchem.compound:",
        sub("^(CID:|pubchem\\.compound:)", "", node_info$ID[idx_pubchem], ignore.case = TRUE)
      )
    }


    # # Extract all Graphics elements inside DataNodes
    # graphics_nodes <-
    #   xml2::xml_find_all(result, ".//d1:DataNode/d1:Graphics", ns = xml2::xml_ns(result))
    #
    # # Extract positions (CenterX, CenterY)
    # graphics_info <-
    #   data.frame(
    #   GraphId = xml2::xml_attr(
    #     xml2::xml_find_all(result, ".//d1:DataNode", ns = xml2::xml_ns(result)),
    #     "GraphId"
    #   ),
    #   CenterX = as.numeric(xml2::xml_attr(graphics_nodes, "CenterX")),
    #   CenterY = as.numeric(xml2::xml_attr(graphics_nodes, "CenterY")),
    #   stringsAsFactors = FALSE
    # )


    # Extract all Interaction elements
    interactions <-
      xml2::xml_find_all(result, ".//d1:Interaction", ns = xml2::xml_ns(result))

    # Extract GraphId and references to connecting nodes
    edge_info <-
      data.frame(
        GraphId = xml2::xml_attr(interactions, "GraphId"),
        Points = sapply(interactions, function(inter) {
          points <- xml2::xml_find_all(inter, ".//d1:Point", ns = xml2::xml_ns(result))
          paste(xml2::xml_attr(points, "GraphRef"), collapse = " -> ")
        }),
        stringsAsFactors = FALSE
      )


    return_result <-
      list(
        pathway_info = pathway_info,
        pathway_description = pathway_description,
        class_info = ontology_info,
        node_info = node_info,
        edge_info = edge_info
      )
    return(return_result)
  }



wikipathway_organisms_info <-
  c(
    "Anopheles gambiae",
    "Arabidopsis thaliana",
    "Bacillus subtilis",
    "Bos taurus",
    "Caenorhabditis elegans",
    "Canis familiaris",
    "Caulobacter vibrioides",
    "Danio rerio",
    "Drosophila melanogaster",
    "Equus caballus",
    "Escherichia coli",
    "Gallus gallus",
    "Gibberella zeae",
    "Homo sapiens" ,
    "Hordeum vulgare",
    "Mus musculus",
    "Mycobacterium tuberculosis",
    "Oryza sativa",
    "Pan troglodytes" ,
    "Plasmodium falciparum",
    "Populus trichocarpa",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Solanum lycopersicum",
    "Sus scrofa",
    "Vitis vinifera",
    "Zea mays"
  )
