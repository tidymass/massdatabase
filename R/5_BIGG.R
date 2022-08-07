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
      tryCatch(
        curl::curl(url = url),
        error = function(e)
          NULL
      )
    open(x)
    out <- readLines(x, warn = FALSE)
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

download_bigg_model <-
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
#' @return model information
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
    out <- readLines(x, warn = FALSE)
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
#' @return universal metabolite information
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
    out <- readLines(x, warn = FALSE)
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
#' @return universal metabolite information
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
    out <- readLines(x, warn = FALSE)
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


#' @title Download BIGG universal metabolites
#' @description Download BIGG universal metabolites
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .
#' @param sleep Default is 1 second.
#' @return bigg_universal_metabolite database.
#' @export

download_bigg_universal_metabolite <-
  function(path = ".",
           sleep = 1) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    metabolite_info <-
      request_bigg_universal_metabolite_info()

    pb <- progress::progress_bar$new(total = nrow(metabolite_info))

    bigg_universal_metabolit_database <-
      seq_along(metabolite_info$bigg_id) %>%
      purrr::map(function(i) {
        pb$tick()
        Sys.sleep(time = sleep)
        request_bigg_universal_metabolite(metabolite_id = metabolite_info$bigg_id[i],
                                          return_form = "data.frame")
      })
    save(
      bigg_universal_metabolit_database,
      file = file.path(path, "bigg_universal_metabolit_database")
    )
  }



#' @title Read the BIGG universal metabolite database from
#' download_bigg_universal_metabolite function
#' @description Read the BIGG universal metabolite database from
#' download_bigg_universal_metabolite function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with
#' download_bigg_universal_metabolite function.
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @export
read_bigg_universal_metabolite <-
  function(path = ".") {
    load(file.path(path, "bigg_universal_metabolit_database"))

    pb <-
      progress::progress_bar$new(total =
                                   length(bigg_universal_metabolit_database))

    bigg_universal_metabolit_database <-
      bigg_universal_metabolit_database %>%
      dplyr::bind_rows() %>%
      as.data.frame()
    return(bigg_universal_metabolit_database)
  }



#' @title Convert BIGG universal metabolites (data.frame,
#' from read_bigg_universal_metabolite)
#' to metID format database
#' @description Convert BIGG universal metabolites (data.frame,
#' from read_bigg_universal_metabolite)
#' to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data data.frame, from read_bigg_universal_metabolite function
#' @param path Default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_bigg_universal2metid <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    ####Formula and mz
    data[which(data == "", arr.ind = TRUE)] <- NA

    data <-
      data %>%
      dplyr::filter(!is.na(formula))

    ####add formula to H or remove H based on charge
    message("Extract formula...")
    pb <-
      progress::progress_bar$new(total = nrow(data))

    formula <-
      seq_len(nrow(data)) %>%
      purrr::map(function(idx) {
        # cat(idx, " ")
        pb$tick()
        charge = data$charges[idx]
        formula <- data$formula[idx]

        if (is.na(charge)) {
          return(data$formula[idx])
        }

        if (charge == "0") {
          return(data$formula[idx])
        }

        if (length(grep("\\{\\}", charge)) == 0) {
          charge <- as.numeric(charge)
          adduct = paste0(ifelse(charge > 0, "M-", "M+"), abs(charge), "H")
          formula <-
            masstools::sum_formula(formula = formula, adduct = adduct)
          return(formula)
        }

        if (length(grep("\\{\\}", charge)) > 0) {
          charge <- as.numeric(stringr::str_split(charge, "\\{\\}")[[1]])
          formula <- stringr::str_split(formula, "\\{\\}")[[1]]
          if (length(charge) != length(formula)) {
            return(NA)
          }

          if (any(charge) == 0) {
            return(formula[which(charge == 0)])
          }

          adduct = paste0(ifelse(charge[1] > 0, "M-", "M+"), abs(charge[1]), "H")
          formula <-
            masstools::sum_formula(formula = formula[1], adduct = adduct)
          return(formula)
        }
      }) %>%
      unlist()

    data$formula <- formula

    data$HMDB <-
      data$HMDB %>%
      purrr::map(function(x) {
        if (is.na(x)) {
          return(x)
        }

        x <- stringr::str_split(x, "\\{\\}")[[1]]
        x <- stringr::str_replace(x, "HMDB", "HMDB00")
        x <- paste(x, collapse = "{}")
        x
      }) %>%
      unlist()

    data <-
      data %>%
      dplyr::filter(!is.na(formula))

    ###add mz
    message("Calculating mz...")
    pb <-
      progress::progress_bar$new(total = nrow(data))
    mz <-
      seq_len(nrow(data)) %>%
      purrr::map(function(i) {
        pb$tick()
        temp <-
          tryCatch(
            Rdisop::getMass(Rdisop::getMolecule(data$formula[i])),
            error = function(e)
              NA
          )
      }) %>%
      unlist()

    data$mz <- mz

    data <-
      data %>%
      dplyr::filter(!is.na(mz))

    data <-
      data %>%
      dplyr::rename(
        BIGG.ID = bigg_id,
        Compound.name = name,
        Formula = formula,
        BIOCYC.ID = BioCyc,
        CHEBI.ID = CHEBI,
        HMDB.ID = HMDB,
        INCHIKEY.ID = InChI_Key,
        KEGG.ID = KEGG,
        METANETX.ID = MetaNetX,
        REACTOME.ID = Reactome,
        SEED.ID = SEED,
        KEGG_DRUG.ID = KEGG_Drug,
        KEGG_GLYCAN.ID = KEGG_Glycan,
        LIPIDMAPS.ID = LipidMaps
      ) %>%
      dplyr::mutate(
        Lab.ID = BIGG.ID,
        CAS.ID = NA,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "BIGG"
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

    data$Synonyms <-
      data$Compound.name %>%
      stringr::str_replace_all("; ", "\\{\\}")

    data$Compound.name <-
      data$Synonyms %>%
      stringr::str_split("\\{\\}") %>%
      lapply(function(x) {
        x[1]
      }) %>%
      unlist()

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = data,
                     file = file.path(temp_file, "data.csv"))


    bigg_ms1 =
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "data.csv",
        source = "BIGG",
        link = "http://bigg.ucsd.edu/",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "data.csv"))
    unlink(temp_file)

    save(bigg_ms1, file = file.path(path, "bigg_ms1"))
    invisible(bigg_ms1)
  }

#' @title Request BIGG universal reaction information
#' @description Request BIGG universal reaction information
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url url http://bigg.ucsd.edu/api/v2/universal/reactions
#' @return universal metabolite information
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
    out <- readLines(x, warn = FALSE)
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









#' @title Request BIGG reaction information for specific species (model)
#' @description Request BIGG reaction information for specific species
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param model_bigg_id please use the model bigg id using the
#' request_bigg_model_info function. Default is iAB_RBC_283 (human)
#' @return reaction information
#' @importFrom curl curl
#' @importFrom stringr str_replace_all str_extract str_replace str_split
#' @importFrom stringr str_trim
#' @export
#' @examples
#' x <- request_bigg_universal_reaction_info()
#' dim(x)
#' head(x)
request_bigg_reaction_info <-
  function(model_bigg_id = "iAB_RBC_283") {
    url <-
      paste0("http://bigg.ucsd.edu/api/v2/models/",
             model_bigg_id,
             "/reactions")
    x <-
      curl::curl(url = url)
    open(x)
    out <- readLines(x, warn = FALSE)
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
#' @return universal metabolite information
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
    out <- readLines(x, warn = FALSE)
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










#' @title Read the BIGG model from download_bigg_model function
#' @description Read the BIGG model from download_bigg_model function
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Default is .. Should be same with
#' download_bigg_model function.
#' @param model model
#' @return A data frame
#' @importFrom magrittr %>%
#' @importFrom plyr dlply .
#' @importFrom readr read_delim
#' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map
#' @importFrom rjson fromJSON
#' @export
read_bigg_model <-
  function(path = ".",
           model = "iND750") {
    result <-
      rjson::fromJSON(file = file.path(path, paste0(model, ".json")))
    return(result)
  }
