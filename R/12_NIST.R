#' @title Convert GNPS data (list) to metID format database
#' @description Convert GNPS data (list) to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_msp_data function
#' @param path Default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_nist2metid <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    message("Extracting MS1 inforamtion...")
    all_names <-
      data %>%
      purrr::map(function(x) {
        x$info$key
      }) %>%
      unlist() %>%
      unique() %>%
      sort()

    progresser <-
      show_progresser(index = seq_along(data),
                      progresser = c(1, seq(10, 100, 10)))

    ms1_info <-
      seq_along(data) %>%
      purrr::map(function(i) {
        # cat(i, " ")
        if (i %in% progresser$idx) {
          message(progresser$progresser[which(i == progresser$idx)], " ",
                  appendLF = FALSE)
        }
        x <- data[[i]]
        x <-
          x$info %>%
          dplyr::arrange(key)

        if (sum(duplicated(x$key)) == 0) {
          x <-
            t(x) %>%
            as.data.frame()
          colnames(x) <- as.character(x[1, ])
          x <- x[-1, , drop = FALSE]
          new_name <-
            setdiff(all_names, colnames(x))
          if (length(new_name) > 0) {
            new_x <-
              matrix(NA, nrow = 1, ncol = length(new_name)) %>%
              as.data.frame()
            colnames(new_x) <- new_name
            x <-
              cbind(x, new_x) %>%
              as.data.frame()
            x <- x[, all_names]
          }
          return(x)
        }

        x <-
          x %>%
          plyr::dlply(.variables = .(key)) %>%
          lapply(function(y) {
            if (nrow(y) == 1) {
              return(y)
            }
            y$value <-
              paste(y$value, collapse = "{}")
            return(y[1, , drop = FALSE])
          }) %>%
          dplyr::bind_rows() %>%
          as.data.frame() %>%
          dplyr::arrange(key)

        x <-
          t(x) %>%
          as.data.frame()
        colnames(x) <- as.character(x[1, ])
        x <- x[-1, , drop = FALSE]
        new_name <-
          setdiff(all_names, colnames(x))
        if (length(new_name) > 0) {
          new_x <-
            matrix(NA, nrow = 1, ncol = length(new_name)) %>%
            as.data.frame()
          colnames(new_x) <- new_name
          x <-
            cbind(x, new_x) %>%
            as.data.frame()
          x <- x[, all_names]

        }
        return(x)
      })

    ms1_info <-
      ms1_info %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    rownames(ms1_info) <- NULL

    message("Done.")

    message("Extracting MS2 inforamtion...")

    spectra_data <-
      seq_along(data) %>%
      purrr::map(function(i) {
        # cat(i, " ")
        if (i %in% progresser$idx) {
          message(progresser$progresser[which(i == progresser$idx)], " ",
                  appendLF = FALSE)
        }
        data[[i]]$spec
      })

    message("Done.")

    message("Organizing...")

    key_names <-
      stringr::str_to_lower(colnames(ms1_info)) %>%
      stringr::str_replace_all("#", "") %>%
      stringr::str_replace_all("-", "") %>%
      stringr::str_replace_all("_", "")

    colnames(ms1_info) <- key_names

    if (any(key_names == "cas")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(CAS.ID = "cas")
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(CAS.ID = NA)
    }

    if (any(colnames(ms1_info) == "inchi")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(INCHI.ID = inchi)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(INCHI.ID = NA)
    }

    if (any(colnames(ms1_info) == "inchikey")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(INCHIKEY.ID = inchikey)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(INCHIKEY.ID = NA)
    }

    if (any(colnames(ms1_info) == "precursortype")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Adduct = precursortype)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Adduct = NA)
    }

    if (any(colnames(ms1_info) == "synon")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Synonyms = Synon)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Synonyms = NA)
    }

    if (any(colnames(ms1_info) == "collisionenergy")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(CE = collisionenergy)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(CE = NA)
    }

    if (any(colnames(ms1_info) == "ionmode")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Polarity = ionmode)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Polarity = NA)
    }

    if (any(colnames(ms1_info) == "precursormz")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Precursor_mz = precursormz)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Precursor_mz = NA)
    }

    if (any(colnames(ms1_info) == "smiles")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(SMILES.ID = smiles)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(SMILES.ID = NA)
    }

    if (any(colnames(ms1_info) == "formula")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Formula = formula)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Formula = NA)
    }

    if (any(colnames(ms1_info) == "instrument")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Instrument = instrument)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Instrument = NA)
    }

    if (any(colnames(ms1_info) == "instrumenttype")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Instrument_type = instrumenttype)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Instrument_type = NA)
    }

    if (any(colnames(ms1_info) == "exactmass")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(mz = exactmass)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(mz = NA)
    }

    if (any(colnames(ms1_info) == "db")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Lab.ID = db)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Lab.ID = NA)
    }

    if (any(colnames(ms1_info) == "nist")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(NIST.ID = nist)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(NIST.ID = NA)
    }

    if (any(colnames(ms1_info) == "name")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Compound.name = name)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Compound.name = NA)
    }

    ms1_info <-
      ms1_info %>%
      dplyr::mutate(
        HMDB.ID = NA,
        KEGG.ID = NA,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "NIST"
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

    ms1_info <-
      ms1_info %>%
      dplyr::mutate(
        Polarity =
          case_when(
            Polarity == "POSITIVE" ~ "Positive",
            Polarity == "NEGATIVE" ~ "Negative",
            Polarity == "P" ~ "Positive",
            Polarity == "N" ~ "Negative",
            TRUE ~ Polarity
          )
      )

    ms1_info <-
      ms1_info %>%
      dplyr::mutate(mz = as.numeric(mz),
                    Precursor_mz = as.numeric(Precursor_mz))

    ###calculate mz
    if (any(is.na(ms1_info$mz))) {
      message("Calculating accurate mass...")

      ms1_info$mz[is.na(ms1_info$mz)] <-
        seq_len(sum(is.na(ms1_info$mz))) %>%
        purrr::map(function(i) {
          formula <-
            ms1_info$Formula[is.na(ms1_info$mz)][i]
          if (!is.na(formula)) {
            temp_mz <-
              tryCatch(
                expr = Rdisop::getMass(Rdisop::getMolecule(formula)),
                error = function(e)
                  NA
              )
            return(temp_mz)
          } else{
            adduct <-
              ms1_info$Adduct[is.na(ms1_info$mz)][i] %>%
              stringr::str_replace("\\]\\+", "") %>%
              stringr::str_replace("\\]\\-", "") %>%
              stringr::str_replace("\\[", "") %>%
              stringr::str_replace("\\]", "")
            Precursor_mz <-
              ms1_info$Precursor_mz[is.na(ms1_info$mz)][i]
            temp_mz <-
              masstools::convert_precursor_mz2accurate_mass(precursor_mz = Precursor_mz,
                                                            adduct = adduct)
            return(temp_mz)
          }
        }) %>%
        unlist() %>%
        as.numeric()
      message("Done.")
    }

    ###new lab.ID
    if (any(is.na(ms1_info$Lab.ID))) {
      ms1_info$Lab.ID[is.na(ms1_info$Lab.ID)] <-
        paste("NIST", 1:sum(is.na(ms1_info$Lab.ID)), sep = "_")
      ms1_info$Lab.ID <-
        masstools::name_duplicated(ms1_info$Lab.ID)
    }

    rownames(ms1_info) <- NULL

    remove_idx <-
      which(is.na(ms1_info$mz))

    if (length(remove_idx) > 0) {
      ms1_info <-
        ms1_info[-remove_idx, ]

      spectra_data <-
        spectra_data[-remove_idx]
    }

    ms1_info[which(ms1_info == "None", arr.ind = TRUE)] <- NA
    ms1_info[which(ms1_info == "", arr.ind = TRUE)] <- NA

    ms1_info2 <-
      ms1_info %>%
      plyr::dlply(.variables = .(Lab.ID)) %>%
      purrr::map(function(y) {
        if (sum(is.na(y$CE)) > 0) {
          y$CE[is.na(y$CE)] <-
            paste("Unknown", 1:length(y$CE[is.na(y$CE)]), sep = "_")
        }
        y
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    ms1_info2 <-
      ms1_info2[match(ms1_info$Lab.ID, ms1_info2$Lab.ID), ]

    progresser <-
      show_progresser(index = seq_along(spectra_data),
                      progresser = c(1, seq(10, 100, 10)))

    spectra_data2 <-
      seq_along(spectra_data) %>%
      purrr::map(function(i) {
        if (i %in% progresser$idx) {
          message(progresser$progresser[which(i == progresser$idx)], " ",
                  appendLF = FALSE)
        }
        x <- spectra_data[[i]]
        x <- list(x)
        names(x) <-
          ms1_info2$CE[i]
        x
      })

    names(spectra_data2) <- ms1_info2$Lab.ID

    ######positive mode
    # ms1_info2$Lab.ID == names(spectra_data2)

    index_pos <- which(ms1_info2$Polarity == "Positive")
    index_neg <- which(ms1_info2$Polarity == "Negative")

    spectra_data_pos <- spectra_data2[index_pos]
    spectra_data_neg <- spectra_data2[index_neg]

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = ms1_info2,
                     file = file.path(temp_file, "ms1_info2.csv"))

    nist_ms2 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "ms1_info2.csv",
        source = "NIST",
        link = "https://www.nist.gov/blogs/taking-measure/shaping-future-nist-mass-spectral-library",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "ms1_info2.csv"))
    unlink(temp_file)

    nist_ms2@spectra.data$Spectra.positive <-
      spectra_data_pos

    nist_ms2@spectra.data$Spectra.negative <-
      spectra_data_neg
    save(nist_ms2, file = file.path(path, "nist_ms2"))
    invisible(nist_ms2)
  }
