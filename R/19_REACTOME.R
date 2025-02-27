#' @title Retrieve Reactome Reactions for a Given Organism
#' @description This function retrieves reaction information for a specified organism from Reactome.
#' @param organism Character. The name of the organism. Default is "Homo sapiens".
#' @return A data frame containing:
#'   \item{reaction_id}{Character. The Reactome reaction ID.}
#'   \item{reaction_link}{Character. The URL link to the reaction in Reactome.}
#'   \item{reaction_name}{Character. The name of the reaction.}
#'   \item{organism}{Character. The organism associated with the reaction.}
#' @importFrom readr read_delim
#' @importFrom dplyr filter
#' @export
#' @author Xiaotao Shen (\email{xiaotao.shen@outlook.com})
#' @examples
#' # Retrieve reaction information for Homo sapiens
#' reactions <- request_reactome_reaction_info("Homo sapiens")
#' head(reactions)

request_reactome_reaction_info <-
  function(organism = "Homo sapiens") {
    organism_new <-
      match.arg(organism, choices = reactome_organisms_info)

    url <-
      "https://reactome.org/download/current/ChEBI2ReactomeReactions.txt"

    result <-
      readr::read_delim(url, delim = "\t", col_names = FALSE)

    result <-
      result[, c(2, 3, 4, 6)]

    colnames(result) <-
      c("reaction_id", "reaction_link", "reaction_name", "organism")

    result <-
      result %>%
      dplyr::filter(organism == organism_new)

    return(result)
  }


#' @title Request one specific reaction information in Reactome
#' @description Request one specific reaction information in Reactome
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param reaction_id reaction id. For example, reaction_id
#' @return Alist.
#' @importFrom curl curl_download
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_reactome_reaction(reaction_id = "R-HSA-8876188")
#' x$reactants
#' x$products

request_reactome_reaction <-
  function(reaction_id = "R-HSA-8876188") {
    url <-
      paste0("https://reactome.org/ContentService/exporter/event/",
             reaction_id,
             ".sbml")

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    curl::curl_download(url = url,
                        destfile = file.path(temp_file, "file.sbml"))

    result <-
      tryCatch(
        parse_reactome_reaction(file_name = file.path(temp_file, "file.sbml")),
        error = function(e) {
          return(NULL)
        }
      )

    invisible(result)
  }


#' @title Parse the sbml reaction data from Reactome
#' @description Parse the sbml reaction data from Reactome
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file_name file name of the data.
#' @return A data frame or list.
#' @importFrom dplyr filter
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @export
parse_reactome_reaction <-
  function(file_name) {
    result <-
      readLines(file_name)

    result <-
      XML::xmlToList(result)

    species <-
      result$model$listOfSpecies %>%
      lapply(function(x) {
        x$.attrs
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame() %>%
      dplyr::select(
        -c(
          boundaryCondition,
          compartment,
          constant,
          hasOnlySubstanceUnits,
          metaid,
          sboTerm
        )
      )

    reactants <-
      result$model$listOfReactions$reaction$listOfReactants %>%
      dplyr::bind_rows() %>%
      as.data.frame() %>%
      dplyr::select(-c(constant, id, sboTerm, stoichiometry))

    products <-
      result$model$listOfReactions$reaction$listOfProducts %>%
      dplyr::bind_rows() %>%
      as.data.frame() %>%
      dplyr::select(-c(constant, id, sboTerm, stoichiometry))

    reactants <-
      reactants %>%
      dplyr::left_join(species, by = c("species" = "id"))

    products <-
      products %>%
      dplyr::left_join(species, by = c("species" = "id"))

    result <-
      list(reactants = reactants, products = products)
    return(result)
  }


#' @title Retrieve Available Organisms in Reactome
#' @description This function retrieves the list of organisms available in the Reactome database.
#' @return A sorted vector of organism names.
#' @importFrom readr read_delim
#' @export
#' @author Xiaotao Shen (\email{xiaotao.shen@outlook.com})
#' @examples
#' organisms <- request_reactome_organisms_info()
#' print(organisms)

request_reactome_organisms_info <-
  function() {
    url <-
      "https://reactome.org/download/current/ReactomePathways.txt"

    result <-
      readr::read_delim(url, delim = "\t", col_names = FALSE)

    colnames(result) <-
      c("pathway_id", "pathway_name", "organism")
    return(sort(unique(result$organism)))
  }


#' @title Retrieve Reactome Pathways for a Given Organism
#' @description This function retrieves pathway information for a specified organism from Reactome.
#' @param organism Character. The name of the organism. Default is "Homo sapiens".
#' @return A data frame containing pathway ID, pathway name, and organism.
#' @importFrom readr read_delim
#' @importFrom dplyr filter
#' @export
#' @author Xiaotao Shen (\email{xiaotao.shen@outlook.com})
#' @examples
#' pathways <- request_reactome_pathway_info("Homo sapiens")
#' head(pathways)

request_reactome_pathway_info <-
  function(organism = "Homo sapiens") {
    organism_new <-
      match.arg(organism, choices = reactome_organisms_info)

    url <-
      "https://reactome.org/download/current/ReactomePathways.txt"

    result <-
      readr::read_delim(url, delim = "\t", col_names = FALSE)

    colnames(result) <-
      c("pathway_id", "pathway_name", "organism")

    result <-
      result %>%
      dplyr::filter(organism == organism_new)

    return(result)

  }


#' @title Retrieve Detailed Reactome Pathway Information
#' @description This function fetches and parses a given Reactome pathway in SBML format.
#' @param pathway_id Character. The Reactome pathway ID (e.g., "R-HSA-5652084").
#' @return A list containing pathway details, including name, description, and involved metabolites/proteins.
#' @importFrom curl curl_download
#' @export
#' @author Xiaotao Shen (\email{xiaotao.shen@outlook.com})
#' @examples
#' pathway_info <- request_reactome_pathway("R-HSA-5652084")
#' print(pathway_info$pathway_name)
#' print(pathway_info$pathway_description)

request_reactome_pathway <-
  function(pathway_id = "R-HSA-5652084") {
    url <-
      paste0("https://reactome.org/ContentService/exporter/event/",
             pathway_id,
             ".sbml")

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    curl::curl_download(url = url,
                        destfile = file.path(temp_file, "file.sbml"))

    result <-
      tryCatch(
        parse_reactome_pathway(file_name = file.path(temp_file, "file.sbml")),
        error = function(e) {
          return(NULL)
        }
      )

    return(result)

  }

#' @title Parse Reactome Pathway SBML File
#' @description This function extracts pathway details, including description, metabolites, and proteins from an SBML file.
#' @param file_name Character. The file path of the SBML file.
#' @return A list containing:
#'   \item{pathway_name}{Character. Name of the pathway.}
#'   \item{pathway_id}{Character. Reactome pathway ID.}
#'   \item{pathway_description}{Character. Text description of the pathway.}
#'   \item{component_info}{Data frame containing metabolites (CHEBI) and proteins (UniProt) involved in the pathway.}
#' @importFrom xml2 read_xml xml_ns xml_attr xml_find_first xml_find_all xml_text
#' @importFrom dplyr mutate case_when
#' @importFrom stringr str_extract str_detect
#' @export
#' @author Xiaotao Shen (\email{xiaotao.shen@outlook.com})
#' @examples
#' # Download a sample SBML file from Reactome and parse it
#' pathway_info <- request_reactome_pathway("R-HSA-5652084")
#' print(pathway_info$pathway_name)
#' print(pathway_info$component_info)
parse_reactome_pathway <-
  function(file_name) {
    sbml_content <- xml2::read_xml(file_name)
    # Extract namespaces
    ns <- xml2::xml_ns(sbml_content)

    # Extract pathway name and ID
    pathway_name <-
      xml2::xml_attr(xml2::xml_find_first(sbml_content, "//d1:model", ns),
                     "name")
    pathway_id <-
      xml2::xml_attr(xml2::xml_find_first(sbml_content, "//d1:model", ns), "id")

    ###pathway description
    notes_node <-
      xml2::xml_find_first(sbml_content,
                           "//d1:model/d1:notes/*[local-name()='p']",
                           ns)

    pathway_description <- xml2::xml_text(notes_node)

    # Extract all species (metabolites & proteins)
    species_nodes <-
      xml2::xml_find_all(sbml_content, "//d1:listOfSpecies/d1:species", ns)

    ####metabolites and proteins in the pathway
    get_external_ids <- function(species_node) {
      # Find all `rdf:li` nodes inside annotation
      xref_nodes <- xml_find_all(species_node, ".//rdf:Bag/rdf:li", ns)

      # Extract `rdf:resource` attributes
      xref_ids <- xml_attr(xref_nodes, "resource")

      # Return as a concatenated string (or list)
      return(paste(xref_ids, collapse = "; "))
    }

    node_id <-
      sapply(species_nodes, get_external_ids)

    node_id <-
      stringr::str_extract(node_id, "CHEBI:[0-9]{3,8}|uniprot:[a-zA-Z0-9]{3,20}")

    node_name = xml_attr(species_nodes, "name")

    component_info <-
      data.frame(node_id = node_id, node_name = node_name) %>%
      dplyr::mutate(node_type = case_when(
        stringr::str_detect(node_id, "CHEBI") ~ "metabolite",
        stringr::str_detect(node_id, "uniprot") ~ "protein"
      ))

    result <-
      list(
        pathway_name = pathway_name,
        pathway_id = pathway_id,
        pathway_description = pathway_description,
        component_info = component_info
      )

    return(result)

  }



reactome_organisms_info <-
  c(
    "Bos taurus",
    "Caenorhabditis elegans",
    "Canis familiaris",
    "Danio rerio",
    "Dictyostelium discoideum",
    "Drosophila melanogaster",
    "Gallus gallus",
    "Homo sapiens",
    "Mus musculus",
    "Mycobacterium tuberculosis",
    "Plasmodium falciparum",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Schizosaccharomyces pombe",
    "Sus scrofa",
    "Xenopus tropicalis"
  )
