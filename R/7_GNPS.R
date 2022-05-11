#' @title Download GNPS MS2 database (msp)
#' @description Download GNPS MS2 database (msp)
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://gnps-external.ucsd.edu/gnpslibrary".
#' @param gnps_library See here: https://gnps-external.ucsd.edu/gnpslibrary.
#' @param path Default is .
#' @return Downloaded files.
#' @importFrom magrittr %>%
#' @export
download_gnps_spectral_library <-
  function(url = "https://gnps-external.ucsd.edu/gnpslibrary",
           gnps_library = "HMDB",
           path = ".") {
    message("Downloading...\n")
    download.file(
      url = paste0(url, "/", gnps_library, ".msp"),
      destfile = file.path(path, paste0(gnps_library, ".msp"))
    )
    message("Done.\n")
  }

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

convert_gnps2metid <-
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
          colnames(x) <- as.character(x[1,])
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
        colnames(x) <- as.character(x[1,])
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

    if (any(colnames(ms1_info) == "INCHI")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(INCHI.ID = INCHI)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(INCHI.ID = NA)
    }

    if (any(colnames(ms1_info) == "INCHIKEY")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(INCHIKEY.ID = INCHIKEY)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(INCHIKEY.ID = NA)
    }

    if (any(colnames(ms1_info) == "PRECURSORTYPE")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Adduct = PRECURSORTYPE)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Adduct = NA)
    }

    if (any(colnames(ms1_info) == "Synon")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Synonyms = Synon)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Synonyms = NA)
    }

    if (any(colnames(ms1_info) == "COLLISIONENERGY")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(CE = COLLISIONENERGY)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(CE = NA)
    }

    if (any(colnames(ms1_info) == "IONMODE")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Polarity = IONMODE)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Polarity = NA)
    }

    if (any(colnames(ms1_info) == "PRECURSORMZ")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Precursor_mz = PRECURSORMZ)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Precursor_mz = NA)
    }

    if (any(colnames(ms1_info) == "SMILES")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(SMILES.ID = SMILES)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(SMILES.ID = NA)
    }

    if (any(colnames(ms1_info) == "FORMULA")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Formula = FORMULA)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Formula = NA)
    }

    if (any(colnames(ms1_info) == "INSTRUMENT")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Instrument = INSTRUMENT)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Instrument = NA)
    }

    if (any(colnames(ms1_info) == "INSTRUMENTTYPE")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Instrument_type = INSTRUMENTTYPE)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Instrument_type = NA)
    }

    if (any(colnames(ms1_info) == "ExactMass")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(mz = ExactMass)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(mz = NA)
    }

    ms1_info$Lab.ID <-
      ms1_info$Comment %>%
      purrr::map(function(x) {
        stringr::str_split(x, ";")[[1]][1] %>%
          stringr::str_replace("^DB#=", "")
      }) %>%
      unlist()

    ms1_info <-
      ms1_info %>%
      dplyr::rename(Compound.name = NAME) %>%
      dplyr::mutate(
        GNPS.ID = Lab.ID,
        CAS.ID = NA,
        HMDB.ID = NA,
        KEGG.ID = NA,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "GNPS"
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

    if (all(is.na(ms1_info$mz))) {
      message("Calculating accurate mass...")
      ms1_info$mz <-
        purrr::map2(
          .x = ms1_info$Precursor_mz,
          .y = ms1_info$Adduct,
          .f = function(x, y) {
            y <-
              y %>%
              stringr::str_replace("\\]\\+", "") %>%
              stringr::str_replace("\\]\\-", "") %>%
              stringr::str_replace("\\[", "") %>%
              stringr::str_replace("\\]", "")
            masstools::convert_precursor_mz2accurate_mass(precursor_mz = x,
                                                          adduct = y)
          }
        ) %>%
        unlist() %>%
        as.numeric()
      message("Done.")
    }

    rownames(ms1_info) <- NULL

    remove_idx <-
      which(is.na(ms1_info$mz))

    if (length(remove_idx) > 0) {
      ms1_info <-
        ms1_info[-remove_idx,]

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
      ms1_info2[match(ms1_info$Lab.ID, ms1_info2$Lab.ID),]

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

    gnps_ms2 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "ms1_info2.csv",
        source = "GNPS",
        link = "https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "ms1_info2.csv"))
    unlink(temp_file)

    gnps_ms2@spectra.data$Spectra.positive <-
      spectra_data_pos

    gnps_ms2@spectra.data$Spectra.negative <-
      spectra_data_neg
    save(gnps_ms2, file = file.path(path, "gnps_ms2"))
    invisible(gnps_ms2)
  }
