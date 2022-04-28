#' @title Convert MoNA data (list) to metID format database
#' @description Convert MoNA data (list) to metID format database
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

convert_mona2metid <-
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

    if (any(colnames(ms1_info) == "InChI")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(INCHI.ID = InChI)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(INCHI.ID = NA)
    }

    if (any(colnames(ms1_info) == "InChIKey")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(INCHIKEY.ID = InChIKey)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(INCHIKEY.ID = NA)
    }

    if (any(colnames(ms1_info) == "Precursor_type")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Adduct = Precursor_type)
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

    if (any(colnames(ms1_info) == "Collision_energy")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(CE = Collision_energy)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(CE = NA)
    }

    if (any(colnames(ms1_info) == "Ion_mode")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Polarity = Ion_mode)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Polarity = NA)
    }

    if (any(colnames(ms1_info) == "PrecursorMZ")) {
      ms1_info <-
        ms1_info %>%
        dplyr::rename(Precursor_mz = PrecursorMZ)
    } else{
      ms1_info <-
        ms1_info %>%
        dplyr::mutate(Precursor_mz = NA)
    }

    ms1_info <-
      ms1_info %>%
      dplyr::rename(
        Lab.ID = `DB#`,
        mz = ExactMass,
        Compound.name = Name
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
                        Polarity == "NEGATIVE" ~ "Negative",
                        Polarity == "P" ~ "Positive",
                        Polarity == "N" ~ "Negative",
                        TRUE ~ Polarity
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

    mona_ms2 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "ms1_info2.csv",
        source = "MoNA",
        link = "https://mona.fiehnlab.ucdavis.edu/",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "ms1_info2.csv"))
    unlink(temp_file)

    mona_ms2@spectra.data$Spectra.positive <-
      spectra_data_pos

    mona_ms2@spectra.data$Spectra.negative <-
      spectra_data_neg
    save(mona_ms2, file = file.path(path, "mona_ms2"))
    invisible(mona_ms2)
  }
