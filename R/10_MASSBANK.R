#' @title Convert MassBank data (list) to metID format database
#' @description Convert MassBank data (list) to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_msp_data function
#' @param source riken or nist
#' @param path Default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_massbank2metid <-
  function(data,
           source = c("riken", "nist"),
           path = ".",
           threads = 5) {
    browser()
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    source <-
      match.arg(source)
    if (source == "nist") {
      convert_massbank2metid_nist(data = data,
                                  path = path,
                                  threads = threads)
    }

    if (source == "riken") {
      convert_massbank2metid_riken(data = data,
                                   path = path,
                                   threads = threads)
    }
  }





#' @title Convert MassBank data (list, from NIST) to metID format database
#' @description Convert MassBank data (list) to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_msp_data function
#' @param path default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_massbank2metid_nist <-
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

    ms1_info <-
      ms1_info %>%
      dplyr::rename(
        Lab.ID = `DB#`,
        mz = ExactMass,
        Compound.name = Name,
        INCHI.ID = InChI,
        INCHIKEY.ID = InChIKey,
        Polarity = Ion_mode,
        Adduct = Precursor_type,
        Precursor_mz = PrecursorMZ,
        SMILES.ID = SMILES,
        Splash = Splash,
        CE = Collision_energy,
        Synonyms = Synon
      ) %>%
      dplyr::mutate(
        MASSBANK.ID = Lab.ID,
        CAS.ID = NA,
        HMDB.ID = NA,
        KEGG.ID = NA,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "MASSBANK"
      ) %>%
      dplyr::select(-c(MW, "Num Peaks", Spectrum_type)) %>%
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
      dplyr::mutate(Polarity =
                      case_when(
                        Polarity == "POSITIVE" ~ "Positive",
                        Polarity == "NEGATIVE" ~ "Negative"
                      ))

    ms1_info <-
      ms1_info %>%
      dplyr::mutate(mz = as.numeric(mz),
                    Precursor_mz = as.numeric(Precursor_mz))

    remove_idx <-
      which(is.na(ms1_info$mz))

    if (length(remove_idx) > 0) {
      ms1_info <-
        ms1_info[-remove_idx, ]

      spectra_data <-
        spectra_data[-remove_idx]
    }

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
    ms1_info2$Lab.ID == names(spectra_data2)

    index_pos <- which(ms1_info2$Polarity == "Positive")
    index_neg <- which(ms1_info2$Polarity == "Negative")

    spectra_data_pos <- spectra_data2[index_pos]
    spectra_data_neg <- spectra_data2[index_neg]

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = ms1_info2,
                     file = file.path(temp_file, "ms1_info2.csv"))

    massbank_ms2 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "ms1_info2.csv",
        source = "MassBank",
        link = "https://massbank.eu/MassBank/",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "ms1_info2.csv"))
    unlink(temp_file)

    massbank_ms2@spectra.data$Spectra.positive <-
      spectra_data_pos

    massbank_ms2@spectra.data$Spectra.negative <-
      spectra_data_neg

    save(massbank_ms2, file = file.path(path, "massbank_ms2"))
    invisible(massbank_ms2)
  }



#' @title Convert MassBank data (list, from RIKEN) to metID format database
#' @description Convert MassBank data (list) to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_msp_data function
#' @param path default is .
#' @param threads threads
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_massbank2metid_riken <-
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

    rownames(ms1_info) <- NULL

    ms1_info <-
      ms1_info %>%
      dplyr::rename(
        Compound.name = NAME,
        INCHI.ID = INCHI,
        INCHIKEY.ID = INCHIKEY,
        Polarity = IONMODE,
        Adduct = ADDUCTIONNAME,
        Precursor_mz = PRECURSORMZ,
        SMILES.ID = SMILES,
        Formula = FORMULA,
        Instrumnet = INSTRUMENT,
        Instrumnet_type = INSTRUMENTTYPE,
        Links = LINKS
      ) %>%
      dplyr::mutate(Lab.ID = paste("MassBank_RIKEN", 1:nrow(ms1_info), sep = "_")) %>%
      dplyr::mutate(
        MASSBANK.ID = Lab.ID,
        CAS.ID = NA,
        HMDB.ID = NA,
        KEGG.ID = NA,
        mz = NA,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "MASSBANK_RIKEN"
      ) %>%
      dplyr::select(-c("Num Peaks")) %>%
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
      dplyr::mutate(Polarity =
                      case_when(
                        Polarity == "POSITIVE" ~ "Positive",
                        Polarity == "NEGATIVE" ~ "Negative"
                      ))


    message("Calculating m/z...")

    ms1_info$mz <-
      seq_along(ms1_info$Formula) %>%
      purrr::map(function(i) {
        if (i %in% progresser$idx) {
          message(progresser$progresser[which(i == progresser$idx)], " ",
                  appendLF = FALSE)
        }
        x <-  ms1_info$Formula[i]
        x <-
          tryCatch(
            Rdisop::getMass(Rdisop::getMolecule(x)),
            error = function(e)
              NA
          )
        x
      }) %>%
      unlist() %>%
      as.numeric()

    message("Done.")

    message("Organizing...")

    ms1_info <-
      ms1_info %>%
      dplyr::mutate(mz = as.numeric(mz),
                    Precursor_mz = as.numeric(Precursor_mz))

    remove_idx <-
      which(is.na(ms1_info$mz))

    if (length(remove_idx) > 0) {
      ms1_info <-
        ms1_info[-remove_idx, ]

      spectra_data <-
        spectra_data[-remove_idx]
    }

    ms1_info[which(ms1_info == "", arr.ind = TRUE)] <- NA

    ms1_info$CE <- NA

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
    ms1_info2$Lab.ID == names(spectra_data2)

    index_pos <- which(ms1_info2$Polarity == "Positive")
    index_neg <- which(ms1_info2$Polarity == "Negative")

    spectra_data_pos <- spectra_data2[index_pos]
    spectra_data_neg <- spectra_data2[index_neg]

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = ms1_info2,
                     file = file.path(temp_file, "ms1_info2.csv"))

    massbank_ms2 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "ms1_info2.csv",
        source = "MassBank",
        link = "https://massbank.eu/MassBank/",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "ms1_info2.csv"))
    unlink(temp_file)

    massbank_ms2@spectra.data$Spectra.positive <-
      spectra_data_pos

    massbank_ms2@spectra.data$Spectra.negative <-
      spectra_data_neg
    save(massbank_ms2, file = file.path(path, "massbank_ms2"))
    invisible(massbank_ms2)
  }
