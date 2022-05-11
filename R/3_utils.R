#' @title Convert species to source
#' @description Convert species to source
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param x source
#' @param match_table match_table
#' @return A data.frame
#' @importFrom magrittr %>%
#' @importFrom dplyr case_when everything select filter
#' @importFrom purrr map map2 walk
#' @importFrom crayon green
#' @export

convert_species2source <-
  function(x,
           match_table) {
    if (is.na(x)) {
      return(
        data.frame(
          From_human = "No",
          From_animal = "No",
          From_microbiota = "No",
          From_archaea = "No",
          From_bacteria = "No",
          From_fungi = "No",
          From_food = "No",
          From_plant = "No",
          From_drug = "No",
          From_environment = "No",
          From_eukaryota = "No",
          From_other = "Yes"
        )
      )
    }

    x <-
      stringr::str_split(x, "\\{\\}")[[1]]

    match_table <-
      as.data.frame(match_table)
    x <-
      as.character(match_table[, 2])[match(x, as.character(match_table[, 1]))]
    x <- x[!is.na(x)]

    if (length(x) == 0) {
      return(
        data.frame(
          From_human = "No",
          From_animal = "No",
          From_microbiota = "No",
          From_archaea = "No",
          From_bacteria = "No",
          From_fungi = "No",
          From_food = "No",
          From_plant = "No",
          From_drug = "No",
          From_environment = "No",
          From_eukaryota = "No",
          From_virus = "No",
          From_other = "Yes"
        )
      )
    }

    ###
    From_human = "No"
    From_animal = "No"
    From_microbiota = "No"
    From_archaea = "No"
    From_bacteria = "No"
    From_fungi = "No"
    From_food = "No"
    From_plant = "No"
    From_drug = "No"
    From_environment = "No"
    From_eukaryota = "No"
    From_virus = "No"
    From_other = "Yes"

    if (any(x == "Human")) {
      From_human <- "Yes"
      From_other = "No"
    }

    if (any(x == "Animalia")) {
      From_animal <- "Yes"
      From_other = "No"
    }

    if (any(x == "Plantae") |
        any(x == "Archaeplastida") | any(x == "Viridiplantae")) {
      From_plant <- "Yes"
      From_other = "No"
    }

    if (any(x == "Bacteria")) {
      From_microbiota <- "Yes"
      From_bacteria <- "Yes"
      From_other = "No"
    }

    if (any(x == "Fungi")) {
      From_microbiota <- "Yes"
      From_fungi <- "Yes"
      From_other = "No"
    }

    if (any(x == "Archaea")) {
      From_microbiota <- "Yes"
      From_archaea <- "Yes"
      From_other = "No"
    }


    if (any(x == "Eukaryota")) {
      From_eukaryota <- "Yes"
      From_other = "No"
    }

    if (any(x == "Food")) {
      From_food <- "Yes"
      From_other = "No"
    }

    if (any(x == "Environment")) {
      From_environment <- "Yes"
      From_other = "No"
    }

    if (any(x == "Virus")) {
      From_virus <- "Yes"
      From_other = "No"
    }

    if (any(x == "Drug")) {
      From_drug <- "Yes"
      From_other = "No"
    }


    if (any(x == "Food_plant")) {
      From_food <- "Yes"
      From_plant <- "Yes"
      From_other = "No"
    }

    data.frame(
      From_human = From_human,
      From_animal = From_animal,
      From_microbiota = From_microbiota,
      From_archaea = From_archaea,
      From_bacteria = From_bacteria,
      From_fungi = From_fungi,
      From_food = From_food,
      From_plant = From_plant,
      From_drug = From_drug,
      From_environment = From_environment,
      From_eukaryota = From_eukaryota,
      From_virus = From_virus,
      From_other = From_other
    )
  }



lipid_class_table <-
  data.frame(
    lipid_class = c(
      "All data",
      "Acylglycerol",
      "Bile Acid",
      "Fatty acid",
      "Long chain alcohol",
      "Long chain aldehyde",
      "Long chain base and Ceramide",
      "Eicosanoid",
      "Ether type lipid",
      "Carotenoid",
      "Coenzyme Q",
      "Vitamin A",
      "Vitamin D",
      "Vitamin E",
      "Vitamin F",
      "Vitamin K",
      "Glycosphingolipid",
      "Glycoglycerolipid and others",
      "Isoprenoid",
      "Lipid peroxide",
      "Lipoamino acid",
      "Lipopolysaccharide",
      "Lipoprotein",
      "Mycolic acid",
      "Glycerophospholipid",
      "Sphingophospholipid",
      "Steroid",
      "Wax"
    ),
    url = c(
      "ALL",
      "NAG",
      "BBA",
      "DFA",
      "DLL",
      "DLD",
      "DLB",
      "XPR",
      "EEL",
      "VCA",
      "VCQ",
      "VVA",
      "VVD",
      "VVE",
      "VVF",
      "VVK",
      "GSG",
      "GCG",
      "IIP",
      "OPO",
      "ALA",
      "CLS",
      "TLP",
      "MMA",
      "PGP",
      "PSP",
      "SST",
      "WWA"
    )
  )


#' @title show_progresser
#' @description show_progresser
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param index index for loop
#' @param progresser progresser
#' @return A data.frame
#' @importFrom magrittr %>%
#' @importFrom dplyr case_when everything select filter
#' @importFrom purrr map map2 walk
#' @importFrom crayon green
#' @export

show_progresser <-
  function(index = 1:1000,
           progresser = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) {
    idx <-
      seq(
        from = 1,
        to = max(index),
        length.out = length(progresser)
      ) %>%
      round()

    data.frame(idx = idx,
               progresser = paste0(progresser, "%"))

  }





keep_one_from_multiple <-
  function(df) {
    df %>%
      apply(1, function(x) {
        x <- as.character(x)
        x <- x[!is.na(x)]
        if (length(x) == 0) {
          return(NA)
        } else{
          x[1]
        }
      })
  }


standard_hmdb_id <-
  function(id) {
    id %>%
      purrr::map(function(x) {
        # cat(x, " ")
        if (is.na(x)) {
          return(NA)
        }
        if (nchar(x) == 9) {
          x %>%
            stringr::str_replace("HMDB", "HMDB00") %>%
            return()
        } else{
          return(x)
        }
      }) %>%
      unlist() %>%
      unname()
  }


update_metid_database_info <-
  function(database,
           ref_database,
           by = c("CAS.ID",
                  "HMDB.ID",
                  "KEGG.ID"),
           combine_columns = c("CAS.ID", "HMDB.ID", "KEGG.ID", "PUBCHEM.ID"),
           new_columns = c("Kingdom", "Super_class", "Class", "Sub_class")) {
    all_names <- c(by, combine_columns) %>%
      unique()
    if (any(!all_names %in% colnames(database@spectra.info))) {
      stop(paste(all_names[which(!all_names %in% colnames(database@spectra.info))], collapse = ", "),
           " are not in database.")
    }

    if (any(!all_names %in% colnames(ref_database@spectra.info))) {
      stop(paste(all_names[which(!all_names %in% colnames(ref_database@spectra.info))], collapse = ", "),
           " are not in ref_database.")
    }

    if (any(!new_columns %in% colnames(ref_database@spectra.info))) {
      stop(paste(new_columns[which(!new_columns %in% colnames(ref_database@spectra.info))], collapse = ", "),
           " are not in ref_database.")
    }

    database@spectra.info <-
      database@spectra.info %>%
        as.data.frame()

    ref_database@spectra.info <-
      ref_database@spectra.info %>%
      as.data.frame()

    idx <-
      by %>%
      purrr::map(function(x) {
        match(database@spectra.info[, x,],
              ref_database@spectra.info[, x],
              incomparables = NA)
      }) %>%
      dplyr::bind_cols()

    idx <-
      idx %>%
      keep_one_from_multiple %>%
      as.numeric()

    ###combine columns
    message("Combining columns...")
    for (x in combine_columns) {
      value <-
        data.frame(database@spectra.info[, x],
                   ref_database@spectra.info[idx, x]) %>%
        keep_one_from_multiple()
      database@spectra.info[, x] <- value
    }
    message("Done.")


    ###new columns
    if (!is.null(new_columns)) {
      if (length(new_columns) > 0) {
        message("Adding new columns...")

        database@spectra.info <-
          database@spectra.info %>%
          dplyr::select(!dplyr::one_of(new_columns))

        for (x in new_columns) {
          value <- ref_database@spectra.info[idx, x]
          database@spectra.info <-
            database@spectra.info %>%
            dplyr::mutate(x = value)
          colnames(database@spectra.info)[ncol(database@spectra.info)] <-
            x
        }

        message("Done.")
      }
    }

    return(database)
  }




merge_same_source <-
  function(source_system) {
    id <-
      source_system %>%
      dplyr::select(tidyselect::ends_with("ID"))

    source <-
      source_system %>%
      dplyr::select(tidyselect::starts_with("From"))

    from <-
      grep("From_", colnames(source), value = TRUE)

    duplicated_from <-
      from %>%
      stringr::str_extract("^From_[a-zA-Z]{1,20}")

    unique_from <-
      duplicated_from[duplicated(duplicated_from)] %>%
      unique()

    if (length(unique_from) == 0) {
      return(source_system)
    }

    new_source <-
      seq_along(unique_from) %>%
      purrr::map(function(i) {
        idx <-
          which(duplicated_from == unique_from[i])
        x <-
        source[, idx] %>%
          apply(1, function(y) {
            y <- as.character(y)
            y <- y[!is.na(y)]
            if (length(y) == 0) {
              return(NA)
            }

            if (any(y == "Yes")) {
              return("Yes")
            }
            return("No")
          })
        unname(x)
      }) %>%
      dplyr::bind_cols() %>%
      as.data.frame()

    colnames(new_source) <- unique_from

    remove_idx <-
      which(duplicated_from %in% unique_from)

    source <-
      cbind(source[, -remove_idx],
            new_source)
    cbind(id, source)
  }





update_metid_database_source_system <-
  function(database,
           source_system,
           by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
           prefer = c("database", "source_system")) {
    prefer <- match.arg(prefer)

    source <-
    source_system %>%
      dplyr::select(tidyselect::starts_with("From_"))

    database@spectra.info <-
      database@spectra.info %>%
      as.data.frame()

    idx <-
      by %>%
      purrr::map(function(x) {
        match(database@spectra.info[, x ,drop = TRUE],
              source_system[, x,drop = TRUE],
              incomparables = NA)
      }) %>%
      dplyr::bind_cols()

    idx <-
      idx %>%
      keep_one_from_multiple %>%
      as.numeric()

    database_source <-
      database@spectra.info %>%
      dplyr::select(tidyselect::starts_with("From"))

    if(ncol(database_source) == 0){
      database_source <- NULL
    }

    if(prefer == "database"){
      final_source <-
        data.frame(database_source, source[idx,])
    }else{
      final_source <-
        data.frame(source[idx,], database_source)
    }


    final_source <-
    merge_same_source(final_source)

    rownames(final_source) <- NULL

    database@spectra.info <-
      database@spectra.info %>%
      dplyr::select(-tidyselect::starts_with("From"))

    database@spectra.info <-
      cbind(database@spectra.info, final_source)

    return(database)
  }
