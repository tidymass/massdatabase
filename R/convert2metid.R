#' @title Convert data (list, from read_msp_data function) to metID format database
#' @description Convert data (list, from read_msp_data function) to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data list, from read_msp_data function
#' @param path Default is .
#' @param threads threads
#' @param Submitter Default is "Xiaotao Shen"
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert2metid <-
  function(data,
           path = ".",
           threads = 5,
           Submitter = "Xiaotao Shen") {
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

    final_names <-
      c(
        "INCHI.ID",
        "INCHIKEY.ID",
        "Adduct",
        "Synonyms",
        "CE",
        "Polarity",
        "Precursor_mz",
        "SMILES.ID",
        "Formula",
        "Instrument",
        "Instrument_type",
        "mz",
        "Compound.name",
        "RT",
        "CAS.ID",
        "HMDB.ID",
        "KEGG.ID"
      )

    raw_names <-
      list(
        c("INCHI"),
        c("INCHIKEY"),
        c("PRECURSORTYPE"),
        c("Synon"),
        c("COLLISIONENERGY"),
        c("IONMODE"),
        c("PRECURSORMZ"),
        c("SMILES"),
        c("FORMULA"),
        c("INSTRUMENT"),
        c("INSTRUMENTTYPE"),
        c("ExactMass"),
        c("NAME"),
        c("RETENTIONTIME"),
        c('CAS'),
        c("HMDB"),
        c("KEGG")
      )

    for (i in seq_along(raw_names)) {
      idx <-
        which(colnames(ms1_info) %in% raw_names[[i]])

      if (length(idx) > 0) {
        colnames(ms1_info)[match(colnames(ms1_info)[idx], colnames(ms1_info))] <-
          final_names[i]
      } else{
        ms1_info <-
          data.frame(ms1_info, new = NA)
        colnames(ms1_info)[ncol(ms1_info)] <-
          final_names[i]
      }
    }

    ms1_info$Lab.ID <-
      masstools::name_duplicated(ms1_info$Compound.name)

    ms1_info <-
      ms1_info %>%
      dplyr::mutate(mz.pos = NA,
                    mz.neg = NA,
                    Submitter = Submitter) %>%
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
        purrr::map(
          seq_len(nrow(ms1_info)),
          .f = function(i) {
            x <- ms1_info$Precursor_mz[i]
            y <- ms1_info$Adduct[i]
            z <- ms1_info$Formula[i]
            mz <-
              tryCatch(
                Rdisop::getMass(Rdisop::getMolecule(z)),
                error = function(e) {
                  NA
                }
              )
            if (is.na(mz)) {
              y <-
                y %>%
                stringr::str_replace("\\]\\+", "") %>%
                stringr::str_replace("\\]\\-", "") %>%
                stringr::str_replace("\\[", "") %>%
                stringr::str_replace("\\]", "")
              mz <-
                masstools::convert_precursor_mz2accurate_mass(precursor_mz = x,
                                                              adduct = y)
            }
            mz
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

    in_house_ms2 <-
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "ms1_info2.csv",
        source = Submitter,
        link = "https://metid.tidymass.org/",
        creater = Submitter,
        email = "shenxt1990@outlook.com",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "ms1_info2.csv"))
    unlink(temp_file)

    in_house_ms2@spectra.data$Spectra.positive <-
      spectra_data_pos

    in_house_ms2@spectra.data$Spectra.negative <-
      spectra_data_neg
    save(in_house_ms2, file = file.path(path, "in_house_ms2"))
    invisible(in_house_ms2)
  }
