#' @title Request BIGG version
#' @description Request BIGG version
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://foodb.ca/compounds".
#' @return version.
#' @importFrom curl curl
#' @importFrom stringr str_replace_all str_split str_trim
#' @importFrom magrittr %>%
#' @export
#' @examples
#' request_bigg_version()

request_bigg_version <-
  function(url = "http://bigg.ucsd.edu/api/v2/database_version") {
    options(warn = -1)
    x <-
      curl::curl(url = url)
    open(x)
    out <- readLines(x)
    close(x)
    # cat(out, sep = "\n")
    out <-
      out %>%
      stringr::str_replace_all("\"", "") %>%
      stringr::str_replace_all("\\{|\\}", "") %>%
      stringr::str_split(pattern = "\\,") %>%
      `[[`(1) %>%
      stringr::str_trim() %>%
      stringr::str_split("\\: ") %>%
      do.call(rbind, .) %>%
      as.data.frame()
    colnames(out) <- c("name", "value")
    return(out)
  }



#' @title Download BIGG model data
#' @description Download BIGG model data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param model_id model_id, for example, iND750
#' @return version.
#' @export

down_bigg_model <-
  function(model_id = "iND750") {
    url <-
      paste0('http://bigg.ucsd.edu/static/models/',
             model_id,
             '.json')
    system(command = paste('curl -O', url))
  }



#' @title Request BIGG model information
#' @description Request BIGG model information
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url url http://bigg.ucsd.edu/api/v2/models
#' @return model inforamtion
#' @importFrom curl curl
#' @importFrom stringr str_replace_all str_extract str_replace str_split 
#' @importFrom stringr str_trim
#' @export
#' @examples 
#' x <- request_bigg_model_info()
#' head(x)
request_bigg_model_info <-
  function(url = "http://bigg.ucsd.edu/api/v2/models") {
    x <-
      curl::curl(url = url)
    open(x)
    out <- readLines(x)
    close(x)
    # cat(out, sep = "\n")
    out <-
      out %>%
      stringr::str_replace_all("\"", "") %>%
      stringr::str_replace_all("\\{results\\: ", "")
    
    model_count <-
      stringr::str_extract(out, "results_count\\: [0-9]{2,3}") %>%
      stringr::str_replace("results_count\\: ", "") %>%
      as.numeric()
    
    out <-
      out %>%
      stringr::str_replace_all("\\, results_count\\: [0-9]{2,3}\\}", "") %>%
      stringr::str_replace(pattern = "^\\[", "") %>%
      stringr::str_replace(pattern = "\\]$", "")
    
    out <-
      out %>%
      stringr::str_split(pattern = "\\}, \\{") %>%
      `[[`(1) %>%
      stringr::str_replace_all("^\\{", "") %>%
      stringr::str_replace_all("\\}$", "") %>%
      lapply(function(x) {
        x <-
          stringr::str_split(x, "\\,")[[1]] %>%
          stringr::str_trim() %>%
          lapply(function(y) {
            stringr::str_split(y, "\\: ")[[1]]
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        
        x1 <- matrix(data = x[, 2], nrow = 1) %>%
          as.data.frame()
        
        colnames(x1) <- x$V1
        x1
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    out$gene_count <- as.numeric(out$gene_count)
    out$reaction_count <- as.numeric(out$reaction_count)
    out$metabolite_count <- as.numeric(out$metabolite_count)
    invisible(out)
  }



#' @title Request BIGG universal metabolite information
#' @description Request BIGG universal metabolite information
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url url http://bigg.ucsd.edu/api/v2/universal/metabolites
#' @return universal metabolite inforamtion
#' @importFrom curl curl
#' @importFrom stringr str_replace_all str_extract str_replace str_split 
#' @importFrom stringr str_trim
#' @export
#' @examples 
#' x <- request_bigg_universal_metabolite_info()
#' dim(x)
#' head(x)

request_bigg_universal_metabolite_info <-
  function(url = "http://bigg.ucsd.edu/api/v2/universal/metabolites") {
    x <-
      curl::curl(url = url)
    open(x)
    out <- readLines(x)
    close(x)
    
    out <-
      out %>%
      stringr::str_replace_all("\"", "") %>%
      stringr::str_replace_all("\\{results\\: ", "")
    
    metabolite_count <-
      stringr::str_extract(out, "results_count\\: [0-9]{2,5}") %>%
      stringr::str_replace("results_count\\: ", "") %>%
      as.numeric()
    
    out <-
      out %>%
      stringr::str_replace_all("\\, results_count\\: [0-9]{2,5}\\}", "") %>%
      stringr::str_replace(pattern = "^\\[", "") %>%
      stringr::str_replace(pattern = "\\]$", "")
    
    out <-
      out %>%
      stringr::str_split(pattern = "\\}, \\{") %>%
      `[[`(1) %>%
      stringr::str_replace_all("^\\{", "") %>%
      stringr::str_replace_all("\\}$", "") %>%
      lapply(function(x) {
        x <-
          x %>%
          stringr::str_split(", name\\: ") %>%
          `[[`(1) %>%
          stringr::str_split(", model_bigg_id\\: ") %>%
          unlist() %>%
          stringr::str_replace_all("bigg_id: ", "")
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(out) <- c("bigg_id", "name", "model_bigg_id")
    
    invisible(out)
  }


#' @title Request BIGG metabolite
#' @description Request BIGG metabolite
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url url http://bigg.ucsd.edu/api/v2/universal/metabolites
#' @param metabolite_id metabolite_id, for example, g3p
#' @param return_form data.frame or list.
#' @return universal metabolite inforamtion
#' @importFrom curl curl
#' @importFrom stringr str_replace_all str_extract str_replace str_split 
#' @importFrom stringr str_trim
#' @export
#' @examples 
#' x <- request_bigg_universal_metabolite(metabolite_id = "g3p", return_form = "list")
#' x
#' y <- request_bigg_universal_metabolite(metabolite_id = "g3p", return_form = "data.frame")
#' y
request_bigg_universal_metabolite <-
  function(url = "http://bigg.ucsd.edu/api/v2/universal/metabolites",
           metabolite_id = "g3p",
           return_form = c("list", "data.frame")) {
    return_form <- match.arg(return_form)
    url <- paste0(url, "/", metabolite_id)
    x <-
      curl::curl(url = url)
    open(x)
    out <- readLines(x)
    close(x)
    
    out <-
      out %>%
      stringr::str_replace_all("\"", "") %>%
      stringr::str_replace_all("\\{results\\: ", "")
    
    out <-
      out %>%
      stringr::str_replace(pattern = "^\\{", "") %>%
      stringr::str_replace(pattern = "\\}$", "")
    
    database_link <-
      stringr::str_extract(string = out, pattern = "database_links: .+\\}, bigg_id:") %>%
      stringr::str_replace("database_links: \\{", "") %>%
      stringr::str_replace("\\}, bigg_id:", "")
    
    database_link <-
      database_link %>%
      stringr::str_split("\\],") %>%
      `[[`(1) %>%
      stringr::str_trim() %>%
      stringr::str_replace("\\]", "") %>%
      lapply(function(x) {
        id <-
          x %>%
          stringr::str_extract_all("id: .{1,50}\\}") %>%
          `[[`(1) %>%
          stringr::str_replace("id: ", "") %>%
          stringr::str_replace("\\}", "") %>%
          paste(collapse = "{}") %>%
          stringr::str_replace_all("CHEBI:", "")
        
        database <-
          x %>%
          stringr::str_split("\\: \\[") %>%
          `[[`(1) %>%
          `[`(1)
        
        data.frame(database = database,
                   id = id)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    bigg_id <-
      stringr::str_extract(string = out, pattern = "\\}, bigg_id: [a-zA-Z0-9_-]{1,100}, formulae:") %>%
      stringr::str_replace(", formulae:", "") %>%
      stringr::str_replace("\\}, bigg_id: ", "")
    
    formula <-
      stringr::str_extract(string = out, pattern = "formulae: .{1,100}, old_identifiers") %>%
      stringr::str_replace("formulae: ", "") %>%
      stringr::str_replace(", old_identifiers", "") %>%
      stringr::str_replace("\\[", "") %>%
      stringr::str_replace("\\]", "") %>%
      stringr::str_replace_all(", ", "{}")
    
    old_identifiers <-
      stringr::str_extract(string = out, pattern = "old_identifiers: .{1,100}, charges:") %>%
      stringr::str_replace("old_identifiers: ", "") %>%
      stringr::str_replace(", charges:", "") %>%
      stringr::str_replace("\\[", "") %>%
      stringr::str_replace("\\]", "") %>%
      stringr::str_replace_all(", ", "{}")
    
    charges <-
      stringr::str_extract(string = out, pattern = "charges: .{1,100}, name:") %>%
      stringr::str_replace("charges: ", "") %>%
      stringr::str_replace(", name:", "") %>%
      stringr::str_replace("\\[", "") %>%
      stringr::str_replace("\\]", "") %>%
      stringr::str_replace_all(", ", "{}")
    
    name <-
      stringr::str_extract(string = out, pattern = "name: .{1,300}, compartments_in_models:") %>%
      stringr::str_replace("name: ", "") %>%
      stringr::str_replace(", compartments_in_models:", "")
    
    compartments_in_models <-
      stringr::str_extract(string = out, pattern = "compartments_in_models: .+\\}\\]") %>%
      stringr::str_replace("compartments_in_models: ", "") %>%
      stringr::str_replace("\\[", "") %>%
      stringr::str_replace("\\]", "")
    
    compartments_in_models <-
      compartments_in_models %>%
      stringr::str_split("\\}, \\{") %>%
      `[[`(1) %>%
      stringr::str_replace("^\\{", "") %>%
      stringr::str_replace("\\}\\]$", "") %>%
      stringr::str_replace("bigg_id: ", "") %>%
      stringr::str_split(", organism: |, model_bigg_id: ") %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(compartments_in_models) <-
      c("bigg_id", "organism", "model_bigg_id")
    
    if (return_form == "list") {
      return_result <-
        list(
          bigg_id = bigg_id,
          name = name,
          formula = formula,
          old_identifiers = old_identifiers,
          charges = charges,
          database_link = database_link,
          compartments_in_models = compartments_in_models
        )
    } else{
      organism <-
        compartments_in_models$organism %>%
        paste(collapse = "{}")
      model_bigg_id <-
        compartments_in_models$model_bigg_id %>%
        paste(collapse = "{}")
      
      return_result <-
        data.frame(
          bigg_id = bigg_id,
          name = name,
          formula = formula,
          old_identifiers = old_identifiers,
          charges = charges,
          organism = organism,
          model_bigg_id = model_bigg_id
        )
      
      database_link <-
        database_link %>%
        dplyr::mutate(
          database = case_when(
            database == "InChI Key" ~ "InChI_Key",
            database == "Reactome Compound" ~ "Reactome",
            database == "SEED Compound" ~ "SEED",
            database == "Human Metabolome Database" ~ "HMDB",
            database == "MetaNetX (MNX) Chemical" ~ "MetaNetX",
            database == "KEGG Compound" ~ "KEGG",
            database == "KEGG Drug" ~ "KEGG_Drug",
            database == "KEGG Glycan" ~ "KEGG_Glycan",
            TRUE ~ database
          )
        )
      
      missing_id <-
        setdiff(
          c(
            "InChI_Key",
            "Reactome",
            "SEED",
            "CHEBI",
            "HMDB",
            "MetaNetX",
            "KEGG",
            "KEGG_Drug",
            "KEGG_Glycan",
            "LipidMaps",
            "BioCyc"
          ),
          database_link$database
        )
      
      if (length(missing_id) > 0) {
        new_database_link <-
          data.frame(database = missing_id, id = NA)
        database_link <-
          rbind(database_link,
                new_database_link) %>%
          dplyr::filter(database != "")
      }
      
      database_link <-
        database_link %>%
        dplyr::arrange(database)
      
      id <-
        as.data.frame(matrix(database_link$id, nrow = 1))
      colnames(id) <- database_link$database
      
      return_result <-
        cbind(return_result,
              id)
    }
    invisible(return_result)
  }



#' @title Request BIGG universal reaction information
#' @description Request BIGG universal reaction information
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url url http://bigg.ucsd.edu/api/v2/universal/reactions
#' @return universal metabolite inforamtion
#' @importFrom curl curl
#' @importFrom stringr str_replace_all str_extract str_replace str_split 
#' @importFrom stringr str_trim
#' @export
#' @examples 
#' x <- request_bigg_universal_reaction_info()
#' dim(x)
#' head(x)
request_bigg_universal_reaction_info <-
  function(url = "http://bigg.ucsd.edu/api/v2/universal/reactions") {
    x <-
      curl::curl(url = url)
    open(x)
    out <- readLines(x)
    close(x)
    
    out <-
      out %>%
      stringr::str_replace_all("\"", "") %>%
      stringr::str_replace_all("\\{results\\: ", "")
    
    reaction_count <-
      stringr::str_extract(out, "results_count\\: [0-9]{2,5}") %>%
      stringr::str_replace("results_count\\: ", "") %>%
      as.numeric()
    
    out <-
      out %>%
      stringr::str_replace_all("\\, results_count\\: [0-9]{2,5}\\}", "") %>%
      stringr::str_replace(pattern = "^\\[", "") %>%
      stringr::str_replace(pattern = "\\]$", "")
    
    out <-
      out %>%
      stringr::str_split(pattern = "\\}, \\{") %>%
      `[[`(1) %>%
      stringr::str_replace_all("^\\{", "") %>%
      stringr::str_replace_all("\\}$", "") %>%
      lapply(function(x) {
        x <-
          x %>%
          stringr::str_split(", name\\: ") %>%
          `[[`(1) %>%
          stringr::str_split(", model_bigg_id\\: ") %>%
          unlist() %>%
          stringr::str_replace_all("bigg_id: ", "")
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(out) <- c("bigg_id", "name", "model_bigg_id")
    
    invisible(out)
  }





#' @title Request BIGG reaction
#' @description Request BIGG reaction
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url url http://bigg.ucsd.edu/api/v2/universal/reactions
#' @param reaction_id reaction_id, for example, ADA
#' @param return_form data.frame or list.
#' @return universal metabolite inforamtion
#' @importFrom curl curl
#' @importFrom stringr str_replace_all str_extract str_replace str_split 
#' @importFrom stringr str_trim
#' @export
#' @examples 
#' x <- request_bigg_universal_reaction(reaction_id = "ADA", return_form = "list")
#' x
#' y <- request_bigg_universal_reaction(reaction_id = "ADA", return_form = "data.frame")
#' y
request_bigg_universal_reaction <-
  function(url = "http://bigg.ucsd.edu/api/v2/universal/reactions",
           reaction_id = "ADA",
           return_form = c("list", "data.frame")) {
    return_form <- match.arg(return_form)
    url <- paste0(url, "/", reaction_id)
    x <-
      curl::curl(url = url)
    open(x)
    out <- readLines(x)
    close(x)
    
    out <-
      out %>%
      stringr::str_replace_all("\"", "") %>%
      stringr::str_replace_all("\\{results\\: ", "")
    
    out <-
      out %>%
      stringr::str_replace(pattern = "^\\{", "") %>%
      stringr::str_replace(pattern = "\\}$", "")
    
    models_containing_reaction <-
      stringr::str_extract(string = out, pattern = "models_containing_reaction: \\[.+\\}\\], reaction_string: ") %>%
      stringr::str_replace("models_containing_reaction: \\[", "") %>%
      stringr::str_replace("\\], reaction_string: ", "")
    
    models_containing_reaction <-
      models_containing_reaction %>%
      stringr::str_split("\\}, \\{") %>%
      `[[`(1) %>%
      stringr::str_trim() %>%
      stringr::str_replace("\\{|\\}", "") %>%
      stringr::str_split(", organism: ") %>%
      lapply(function(x) {
        x <-
          stringr::str_replace(x, "bigg_id: ", "")
        x
        
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(models_containing_reaction) <-
      c("bigg_id", "organism")
    
    bigg_id <-
      stringr::str_extract(string = out, pattern = "\\]\\}\\, bigg_id\\: .{1,20}, old_identifiers\\: ") %>%
      stringr::str_replace("\\]\\}, bigg_id\\: ", "") %>%
      stringr::str_replace("\\, old_identifiers\\: ", "")
    
    old_identifiers <-
      stringr::str_extract(string = out, pattern = "old_identifiers: .{1,50}, name: ") %>%
      stringr::str_replace("old_identifiers: ", "") %>%
      stringr::str_replace(", name: ", "") %>%
      stringr::str_replace("\\[", "") %>%
      stringr::str_replace("\\]", "") %>%
      stringr::str_replace_all(", ", "{}")
    
    pseudoreaction <-
      stringr::str_extract(string = out, pattern = "pseudoreaction: .{1,20}") %>%
      stringr::str_replace("pseudoreaction: ", "")
    
    name <-
      stringr::str_extract(string = out, pattern = "name: .{1,100}, pseudoreaction:") %>%
      stringr::str_replace("name: ", "") %>%
      stringr::str_replace(", pseudoreaction:", "")
    
    reaction_string <-
      stringr::str_extract(string = out, pattern = "reaction_string: .{1,100}, metabolites: ") %>%
      stringr::str_replace("reaction_string: ", "") %>%
      stringr::str_replace(", metabolites: ", "")
    
    metabolites <-
      stringr::str_extract(string = out, pattern = "metabolites: \\[.+, database_links: ") %>%
      stringr::str_replace("metabolites: \\[\\{", "") %>%
      stringr::str_replace("\\}\\], database_links: ", "")
    
    metabolites <-
      metabolites %>%
      stringr::str_split("\\}, \\{") %>%
      `[[`(1) %>%
      stringr::str_trim() %>%
      stringr::str_split(", name: |, compartment_bigg_id: |, stoichiometry: ") %>%
      lapply(function(x) {
        x <-
          stringr::str_replace_all(x, "bigg_id: ", "")
        x
        
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(metabolites) <-
      c("bigg_id", "name", "compartment_bigg_id", "stoichiometry")
    
    database_link <-
      stringr::str_extract(string = out, pattern = "database_links:.+\\}, bigg_id: ") %>%
      stringr::str_replace("database_links: \\{", "") %>%
      stringr::str_replace("\\]\\}, bigg_id:", "") %>%
      stringr::str_trim()
    
    database_link <-
      database_link %>%
      stringr::str_split("\\],") %>%
      `[[`(1) %>%
      stringr::str_trim() %>%
      stringr::str_replace("\\]", "") %>%
      lapply(function(x) {
        id <-
          x %>%
          stringr::str_extract_all("id: .{1,50}\\}") %>%
          `[[`(1) %>%
          stringr::str_replace("id: ", "") %>%
          stringr::str_replace("\\}", "") %>%
          paste(collapse = "{}") %>%
          stringr::str_replace_all("CHEBI:", "")
        
        database <-
          x %>%
          stringr::str_split("\\: \\[") %>%
          `[[`(1) %>%
          `[`(1)
        
        data.frame(database = database,
                   id = id)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    
    if (return_form == "list") {
      return_result <-
        list(
          models_containing_reaction = models_containing_reaction,
          bigg_id = bigg_id,
          old_identifiers = old_identifiers,
          pseudoreaction = pseudoreaction,
          name = name,
          reaction_string = reaction_string,
          metabolites = metabolites,
          database_link = database_link
        )
    } else{
      models_bigg_id <-
        models_containing_reaction$bigg_id %>%
        paste(collapse = "{}")
      models_bigg_organism <-
        models_containing_reaction$organism %>%
        paste(collapse = "{}")
      
      metabolite_bigg_id <-
        metabolites$bigg_id %>%
        paste(collapse = "{}")
      
      metabolite_name <-
        metabolites$name %>%
        paste(collapse = "{}")
      
      metabolites_compartment_bigg_id <-
        metabolites$compartment_bigg_id %>%
        paste(collapse = "{}")
      
      metabolites_stoichiometry <-
        metabolites$stoichiometry %>%
        paste(collapse = "{}")
      
      return_result <-
        data.frame(
          bigg_id = bigg_id,
          old_identifiers = old_identifiers,
          pseudoreaction = pseudoreaction,
          name = name,
          reaction_string = reaction_string,
          models_bigg_id = models_bigg_id,
          models_bigg_organism = models_bigg_organism,
          metabolite_bigg_id = metabolite_bigg_id,
          metabolite_name = metabolite_name,
          metabolites_compartment_bigg_id = metabolites_compartment_bigg_id,
          metabolites_stoichiometry = metabolites_stoichiometry
        )
      
      database_link <-
        database_link %>%
        dplyr::mutate(
          database = case_when(
            database == "EC Number" ~ "EC_Number",
            database == "Reactome Reaction" ~ "Reactome",
            database == "SEED Reaction" ~ "SEED",
            database == "Human Metabolome Database" ~ "HMDB",
            database == "MetaNetX (MNX) Equation" ~ "MetaNetX",
            database == "KEGG Reaction" ~ "KEGG",
            TRUE ~ database
          )
        )
      
      missing_id <-
        setdiff(
          c(
            "RHEA",
            "KEGG",
            "BioCyc",
            "Reactome",
            "Reactome",
            "SEED",
            "MetaNetX"
          ),
          database_link$database
        )
      
      if (length(missing_id) > 0) {
        new_database_link <-
          data.frame(database = missing_id, id = NA)
        database_link <-
          rbind(database_link,
                new_database_link)
      }
      
      database_link <-
        database_link %>%
        dplyr::arrange(database)
      
      id <-
        as.data.frame(matrix(database_link$id, nrow = 1))
      colnames(id) <- database_link$database
      
      return_result <-
        cbind(return_result,
              id)
    }
    invisible(return_result)
  }
