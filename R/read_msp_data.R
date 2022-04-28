#---------------------------------------------------------------------------
#' @title read_msp_data
#' @description Read MSP data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param source massbank, mona.
#' @param threads threads
#' @return Return ms2 data. This is a list.
#' @export

read_msp_data <-
  function(file,
           source = c("massbank", "mona", "gnps", "nist"),
           threads = 5) {
    source <- match.arg(source)
    msp_data <- readr::read_lines(file)
    if (source == "massbank") {
      data <-
        read_msp_data_massbank(file = file,
                               threads = threads)
    }

    if (source == "mona") {
      data <-
        read_msp_data_mona(file = file,
                           threads = threads)
    }

    if (source == "gnps") {
      data <-
        read_msp_data_gnps(file = file,
                           threads = threads)
    }

    if (source == "nist") {
      data <-
        read_msp_data_nist(file = file,
                           threads = threads)
    }

    return(data)

  }

#---------------------------------------------------------------------------
#' @title read_msp_data_massbank
#' @description Read MSP data from massbank
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param threads threads
#' @return Return ms2 data. This is a list.
#' @export

read_msp_data_massbank <-
  function(file,
           threads = 5) {
    msp_data <- readr::read_lines(file)
    if (length(grep("BEGIN IONS", msp_data)) > 0) {
      msp_data <- msp_data[msp_data != ""]
      temp_idx1 <- grep("BEGIN IONS", msp_data)
      temp_idx2 <- grep("END IONS", msp_data)
      if (length(temp_idx2) < length(temp_idx1)) {
        temp_idx2 <- c(temp_idx2, length(msp_data))
      }

      temp_idx <-
        purrr::map2(.x = temp_idx1, temp_idx2, function(x, y) {
          c(x + 1, y - 1)
        })
      # future::plan(strategy = future::multisession, workers = threads)
      ms2_spec <- purrr::map(
        .x = temp_idx,
        .f = function(x) {
          temp_spec <- msp_data[x[1]:x[2]]
          temp_spec <- temp_spec
          spec_info <-
            temp_spec[stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec <-
            temp_spec[!stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec_info <- stringr::str_split(spec_info, "\\=") %>%
            do.call(rbind, .)
          mz <-
            as.numeric(spec_info[grep("MASS|MZ", spec_info[, 1]), 2])
          rt <-
            as.numeric(spec_info[grep("RT|RETETION", spec_info[, 1]), 2])
          spec_info <- c(mz = mz, rt = rt)

          spec <- purrr::map(
            .x = spec,
            .f = function(x) {
              stringr::str_split(x, " ")[[1]] %>% as.numeric()
            }
          ) %>%
            do.call(rbind, .) %>%
            as.data.frame()

          spec <-
            spec %>% as.matrix()

          colnames(spec) <- c("mz", "intensity")

          spec <- list(info = spec_info, spec = spec)
          spec
        }
      )
      return(ms2_spec)
    } else{
      # n.tot <- length(msp_data)
      n_null <- which(msp_data == '')
      temp_idx1 <- c(1, n_null[-length(n_null)])
      temp_idx2 <- n_null - 1

      temp_idx <-
        1:length(temp_idx1) %>%
        purrr::map(function(i) {
          c(temp_idx1[i], temp_idx2[i])
        })

      progresser_index <-
        seq(from = 1,
            to = length(temp_idx),
            length.out = 10) %>%
        round()

      progresser_index <-
        temp_idx[progresser_index] %>%
        lapply(function(y) {
          y[1]
        }) %>%
        unlist()

      progresser_index <-
        data.frame(
          idx = seq(10, 100, by = 10),
          progresser_index = progresser_index,
          stringsAsFactors = FALSE
        )

      info_spec <-
        purrr::map(
          temp_idx,
          .f =
            function(idx) {
              if (idx[1] %in% progresser_index$progresser_index) {
                message(paste0(progresser_index$idx[idx[1] == progresser_index$progresser_index], "% "),
                        appendLF = FALSE)
              }
              temp_msp_data <- msp_data[idx[1]:idx[2]]
              temp_msp_data <-
                temp_msp_data[temp_msp_data != ""]
              info_idx <-
                grep("[A-Za-z]", temp_msp_data)

              temp_info <-
                temp_msp_data[info_idx] %>%
                stringr::str_split(":", n = 2) %>%
                purrr::map(function(x) {
                  stringr::str_trim(x)
                }) %>%
                do.call(rbind, .) %>%
                as.data.frame()

              rownames(temp_info) <- NULL
              colnames(temp_info) <- c("key", "value")

              temp_spec <- temp_msp_data[-info_idx]

              if (length(temp_spec) != 0) {
                if (length(grep(" ", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(temp_spec, split = ' ') %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                if (length(grep("\t", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(x = temp_spec, split = "\t") %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                temp_spec <-
                  temp_spec %>%
                  apply(1, as.numeric) %>%
                  t() %>%
                  as.data.frame()

                colnames(temp_spec) <-
                  c('mz', 'intensity')

                rownames(temp_spec) <- NULL

                temp_spec <-
                  temp_spec %>%
                  dplyr::filter(intensity != 0)
              } else{
                temp_spec <- NULL
              }

              list('info' = temp_info,
                   'spec' = temp_spec)
            }
        )
    }

    remove_idx <-
      info_spec %>%
      purrr::map(function(x) {
        is.null(x$spec)
      }) %>%
      unlist() %>%
      which()

    if (length(remove_idx) > 0) {
      info_spec <- info_spec[-remove_idx]
    }

    info_spec <- info_spec
  }









#---------------------------------------------------------------------------
#' @title read_msp_data_mona
#' @description Read MSP data from MoNA
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param threads threads
#' @return Return ms2 data. This is a list.
#' @export

read_msp_data_mona <-
  function(file,
           threads = 5) {
    msp_data <- readr::read_lines(file)
    if (length(grep("BEGIN IONS", msp_data)) > 0) {
      msp_data <- msp_data[msp_data != ""]
      temp_idx1 <- grep("BEGIN IONS", msp_data)
      temp_idx2 <- grep("END IONS", msp_data)
      if (length(temp_idx2) < length(temp_idx1)) {
        temp_idx2 <- c(temp_idx2, length(msp_data))
      }

      temp_idx <-
        purrr::map2(.x = temp_idx1, temp_idx2, function(x, y) {
          c(x + 1, y - 1)
        })
      # future::plan(strategy = future::multisession, workers = threads)
      ms2_spec <- purrr::map(
        .x = temp_idx,
        .f = function(x) {
          temp_spec <- msp_data[x[1]:x[2]]
          temp_spec <- temp_spec
          spec_info <-
            temp_spec[stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec <-
            temp_spec[!stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec_info <- stringr::str_split(spec_info, "\\=") %>%
            do.call(rbind, .)
          mz <-
            as.numeric(spec_info[grep("MASS|MZ", spec_info[, 1]), 2])
          rt <-
            as.numeric(spec_info[grep("RT|RETETION", spec_info[, 1]), 2])
          spec_info <- c(mz = mz, rt = rt)

          spec <- purrr::map(
            .x = spec,
            .f = function(x) {
              stringr::str_split(x, " ")[[1]] %>% as.numeric()
            }
          ) %>%
            do.call(rbind, .) %>%
            as.data.frame()

          spec <-
            spec %>% as.matrix()

          colnames(spec) <- c("mz", "intensity")

          spec <- list(info = spec_info, spec = spec)
          spec
        }
      )
      return(ms2_spec)
    } else{
      # n.tot <- length(msp_data)
      n_null <- which(msp_data == '')
      temp_idx1 <- c(1, n_null[-length(n_null)])
      temp_idx2 <- n_null - 1

      temp_idx <-
        1:length(temp_idx1) %>%
        purrr::map(function(i) {
          x <-
            c(temp_idx1[i], temp_idx2[i])
        })

      remove_idx <-
        temp_idx %>%
        purrr::map(function(x) {
          (x[1] == x[2])
        }) %>%
        unlist() %>%
        which()

      if (length(remove_idx) > 0) {
        temp_idx <-
          temp_idx[-remove_idx]
      }

      progresser_index <-
        seq(from = 1,
            to = length(temp_idx),
            length.out = 10) %>%
        round()

      progresser_index <-
        temp_idx[progresser_index] %>%
        lapply(function(y) {
          y[1]
        }) %>%
        unlist()

      progresser_index <-
        data.frame(
          idx = seq(10, 100, by = 10),
          progresser_index = progresser_index,
          stringsAsFactors = FALSE
        )

      info_spec <-
        purrr::map(
          temp_idx,
          .f =
            function(idx) {
              # cat(idx[1], " ")
              if (idx[1] %in% progresser_index$progresser_index) {
                message(paste0(progresser_index$idx[idx[1] == progresser_index$progresser_index], "% "),
                        appendLF = FALSE)
              }
              temp_msp_data <- msp_data[idx[1]:idx[2]]
              temp_msp_data <-
                temp_msp_data[temp_msp_data != ""]
              info_idx <-
                grep("[A-Za-z]", temp_msp_data)

              temp_info <-
                temp_msp_data[info_idx] %>%
                stringr::str_split(":", n = 2) %>%
                purrr::map(function(x) {
                  stringr::str_trim(x)
                }) %>%
                do.call(rbind, .) %>%
                as.data.frame()

              rownames(temp_info) <- NULL
              colnames(temp_info) <- c("key", "value")

              temp_spec <- temp_msp_data[-info_idx]

              if (length(temp_spec) != 0) {
                if (length(grep(" ", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(temp_spec, split = ' ') %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                if (length(grep("\t", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(x = temp_spec, split = "\t") %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                temp_spec <-
                  temp_spec %>%
                  apply(1, as.numeric) %>%
                  t() %>%
                  as.data.frame()

                colnames(temp_spec) <-
                  c('mz', 'intensity')

                rownames(temp_spec) <- NULL

                temp_spec <-
                  temp_spec %>%
                  dplyr::filter(intensity != 0)
              } else{
                temp_spec <- NULL
              }

              list('info' = temp_info,
                   'spec' = temp_spec)
            }
        )
    }

    remove_idx <-
      info_spec %>%
      purrr::map(function(x) {
        is.null(x$spec)
      }) %>%
      unlist() %>%
      which()

    if (length(remove_idx) > 0) {
      info_spec <- info_spec[-remove_idx]
    }

    info_spec <- info_spec
  }







#---------------------------------------------------------------------------
#' @title read_msp_data_gnps
#' @description Read MSP data from GNPS
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param threads threads
#' @return Return ms2 data. This is a list.
#' @export

read_msp_data_gnps <-
  function(file,
           threads = 5) {
    msp_data <- readr::read_lines(file, progress = TRUE)
    if (length(grep("BEGIN IONS", msp_data)) > 0) {
      msp_data <- msp_data[msp_data != ""]
      temp_idx1 <- grep("BEGIN IONS", msp_data)
      temp_idx2 <- grep("END IONS", msp_data)
      if (length(temp_idx2) < length(temp_idx1)) {
        temp_idx2 <- c(temp_idx2, length(msp_data))
      }

      temp_idx <-
        purrr::map2(.x = temp_idx1, temp_idx2, function(x, y) {
          c(x + 1, y - 1)
        })
      # future::plan(strategy = future::multisession, workers = threads)
      ms2_spec <- purrr::map(
        .x = temp_idx,
        .f = function(x) {
          temp_spec <- msp_data[x[1]:x[2]]
          temp_spec <- temp_spec
          spec_info <-
            temp_spec[stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec <-
            temp_spec[!stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec_info <- stringr::str_split(spec_info, "\\=") %>%
            do.call(rbind, .)
          mz <-
            as.numeric(spec_info[grep("MASS|MZ", spec_info[, 1]), 2])
          rt <-
            as.numeric(spec_info[grep("RT|RETETION", spec_info[, 1]), 2])
          spec_info <- c(mz = mz, rt = rt)

          spec <- purrr::map(
            .x = spec,
            .f = function(x) {
              stringr::str_split(x, " ")[[1]] %>% as.numeric()
            }
          ) %>%
            do.call(rbind, .) %>%
            as.data.frame()

          spec <-
            spec %>% as.matrix()

          colnames(spec) <- c("mz", "intensity")

          spec <- list(info = spec_info, spec = spec)
          spec
        }
      )
      return(ms2_spec)
    } else{
      # n.tot <- length(msp_data)
      n_null <- which(msp_data == '')
      temp_idx1 <- c(1, n_null[-length(n_null)])
      temp_idx2 <- n_null - 1

      temp_idx <-
        1:length(temp_idx1) %>%
        purrr::map(function(i) {
          x <-
            c(temp_idx1[i], temp_idx2[i])
        })

      remove_idx <-
        temp_idx %>%
        purrr::map(function(x) {
          (x[1] == x[2])
        }) %>%
        unlist() %>%
        which()

      if (length(remove_idx) > 0) {
        temp_idx <-
          temp_idx[-remove_idx]
      }

      progresser_index <-
        seq(from = 1,
            to = length(temp_idx),
            length.out = 10) %>%
        round()

      progresser_index <-
        temp_idx[progresser_index] %>%
        lapply(function(y) {
          y[1]
        }) %>%
        unlist()

      progresser_index <-
        data.frame(
          idx = seq(10, 100, by = 10),
          progresser_index = progresser_index,
          stringsAsFactors = FALSE
        )

      info_spec <-
        purrr::map(
          temp_idx,
          .f =
            function(idx) {
              # cat(idx[1], " ")
              if (idx[1] %in% progresser_index$progresser_index) {
                message(paste0(progresser_index$idx[idx[1] == progresser_index$progresser_index], "% "),
                        appendLF = FALSE)
              }
              temp_msp_data <- msp_data[idx[1]:idx[2]]
              temp_msp_data <-
                temp_msp_data[temp_msp_data != ""]
              info_idx <-
                grep("[A-Za-z]", temp_msp_data)

              temp_info <-
                temp_msp_data[info_idx] %>%
                stringr::str_split(":", n = 2) %>%
                purrr::map(function(x) {
                  stringr::str_trim(x)
                }) %>%
                do.call(rbind, .) %>%
                as.data.frame()

              rownames(temp_info) <- NULL
              colnames(temp_info) <- c("key", "value")

              temp_spec <- temp_msp_data[-info_idx]

              if (length(temp_spec) != 0) {
                if (length(grep(" ", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(temp_spec, split = ' ') %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                if (length(grep("\t", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(x = temp_spec, split = "\t") %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                temp_spec <-
                  temp_spec %>%
                  apply(1, as.numeric) %>%
                  t() %>%
                  as.data.frame()

                colnames(temp_spec) <-
                  c('mz', 'intensity')

                rownames(temp_spec) <- NULL

                temp_spec <-
                  temp_spec %>%
                  dplyr::filter(intensity != 0)
              } else{
                temp_spec <- NULL
              }

              list('info' = temp_info,
                   'spec' = temp_spec)
            }
        )
    }

    remove_idx <-
      info_spec %>%
      purrr::map(function(x) {
        is.null(x$spec)
      }) %>%
      unlist() %>%
      which()

    if (length(remove_idx) > 0) {
      info_spec <- info_spec[-remove_idx]
    }

    info_spec <- info_spec
  }







#---------------------------------------------------------------------------
#' @title read_msp_data_nist
#' @description Read MSP data from NIST
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of ms2 files. MS2 file must be msp format.
#' @param threads threads
#' @return Return ms2 data. This is a list.
#' @export

read_msp_data_nist <-
  function(file,
           threads = 5) {
    msp_data <- readr::read_lines(file, progress = TRUE)

    if (sum(stringr::str_detect(string = msp_data[1:1000], pattern = "BEGIN IONS")) > 0) {
      msp_data <- msp_data[msp_data != ""]
      temp_idx1 <- grep("BEGIN IONS", msp_data)
      temp_idx2 <- grep("END IONS", msp_data)
      if (length(temp_idx2) < length(temp_idx1)) {
        temp_idx2 <- c(temp_idx2, length(msp_data))
      }

      temp_idx <-
        purrr::map2(.x = temp_idx1, temp_idx2, function(x, y) {
          c(x + 1, y - 1)
        })
      # future::plan(strategy = future::multisession, workers = threads)
      ms2_spec <- purrr::map(
        .x = temp_idx,
        .f = function(x) {
          temp_spec <- msp_data[x[1]:x[2]]
          temp_spec <- temp_spec
          spec_info <-
            temp_spec[stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec <-
            temp_spec[!stringr::str_detect(temp_spec, "[A-Za-z]")]
          spec_info <- stringr::str_split(spec_info, "\\=") %>%
            do.call(rbind, .)
          mz <-
            as.numeric(spec_info[grep("MASS|MZ", spec_info[, 1]), 2])
          rt <-
            as.numeric(spec_info[grep("RT|RETETION", spec_info[, 1]), 2])
          spec_info <- c(mz = mz, rt = rt)

          spec <- purrr::map(
            .x = spec,
            .f = function(x) {
              stringr::str_split(x, " ")[[1]] %>% as.numeric()
            }
          ) %>%
            do.call(rbind, .) %>%
            as.data.frame()

          spec <-
            spec %>% as.matrix()

          colnames(spec) <- c("mz", "intensity")

          spec <- list(info = spec_info, spec = spec)
          spec
        }
      )
      return(ms2_spec)
    } else{
      # n.tot <- length(msp_data)
      n_null <- which(msp_data == '')
      temp_idx1 <- c(1, n_null[-length(n_null)])
      temp_idx2 <- n_null - 1

      temp_idx <-
        1:length(temp_idx1) %>%
        purrr::map(function(i) {
          x <-
            c(temp_idx1[i], temp_idx2[i])
        })

      remove_idx <-
        temp_idx %>%
        purrr::map(function(x) {
          (x[1] == x[2])
        }) %>%
        unlist() %>%
        which()

      if (length(remove_idx) > 0) {
        temp_idx <-
          temp_idx[-remove_idx]
      }

      progresser_index <-
        seq(from = 1,
            to = length(temp_idx),
            length.out = 10) %>%
        round()

      progresser_index <-
        temp_idx[progresser_index] %>%
        lapply(function(y) {
          y[1]
        }) %>%
        unlist()

      progresser_index <-
        data.frame(
          idx = seq(10, 100, by = 10),
          progresser_index = progresser_index,
          stringsAsFactors = FALSE
        )

      info_spec <-
        purrr::map(
          temp_idx,
          .f =
            function(idx) {
              # cat(idx[1], " ")
              if (idx[1] %in% progresser_index$progresser_index) {
                message(paste0(progresser_index$idx[idx[1] == progresser_index$progresser_index], "% "),
                        appendLF = FALSE)
              }
              temp_msp_data <- msp_data[idx[1]:idx[2]]
              temp_msp_data <-
                temp_msp_data[temp_msp_data != ""]

              if (stringr::str_detect(tail(temp_msp_data, 1), "\t")) {
                temp_spec <-
                  stringr::str_extract(
                    temp_msp_data,
                    "[0-9]{1,4}[\\.]{0,1}[0-9]{0,10}\t[0-9]{1,4}[\\.]{0,1}[0-9]{0,10}"
                  )
              } else{
                temp_spec <-
                  stringr::str_extract(
                    temp_msp_data,
                    "[0-9]{1,4}[\\.]{0,1}[0-9]{0,10} [0-9]{1,4}[\\.]{0,1}[0-9]{0,10}"
                  )
              }

              temp_info <-
                temp_msp_data[is.na(temp_spec)] %>%
                stringr::str_split(":", n = 2) %>%
                purrr::map(function(x) {
                  stringr::str_trim(x)
                }) %>%
                do.call(rbind, .) %>%
                as.data.frame()

              rownames(temp_info) <- NULL
              colnames(temp_info) <- c("key", "value")

              temp_spec <- temp_spec[!is.na(temp_spec)]

              if (length(temp_spec) != 0) {
                if (length(grep(" ", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(temp_spec, split = ' ') %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                if (length(grep("\t", temp_spec[1])) == 1) {
                  temp_spec <-
                    strsplit(x = temp_spec, split = "\t") %>%
                    do.call(rbind, .) %>%
                    as.data.frame()
                }

                temp_spec <-
                  temp_spec %>%
                  apply(1, as.numeric) %>%
                  t() %>%
                  as.data.frame()

                colnames(temp_spec) <-
                  c('mz', 'intensity')

                rownames(temp_spec) <- NULL

                temp_spec <-
                  temp_spec %>%
                  dplyr::filter(intensity != 0)
              } else{
                temp_spec <- NULL
              }

              list('info' = temp_info,
                   'spec' = temp_spec)
            }
        )
    }

    remove_idx <-
      info_spec %>%
      purrr::map(function(x) {
        is.null(x$spec)
      }) %>%
      unlist() %>%
      which()

    if (length(remove_idx) > 0) {
      info_spec <- info_spec[-remove_idx]
    }

    info_spec <- info_spec
  }
