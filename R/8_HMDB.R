#' @title Request one specific the metabolite information in HMDB
#' @description Request one specific the metabolite information in HMDB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://hmdb.ca/metabolites".
#' @param metabolite_id metabolite id. For example, HMDB0000001
#' @param return_form data.frame or list.
#' @return A data frame or list.
#' @importFrom XML xmlTreeParse xmlToList
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_hmdb_metabolite(metabolite_id = "HMDB0000001", return_form = "list")
#' x[1:2]
#' y = request_hmdb_metabolite(metabolite_id = "HMDB0000001", return_form = "data.frame")
#' head(y)

request_hmdb_metabolite <-
  function(url = "https://hmdb.ca/metabolites",
           metabolite_id = "HMDB0000001",
           return_form = c("list" , "data.frame")) {
    return_form <- match.arg(return_form)
    result <-
      readLines(paste0(url, "/", metabolite_id, ".xml"),
                warn = FALSE)
    result <-
      XML::xmlTreeParse(file = result, asText = TRUE)
    result <-
      XML::xmlToList(result)

    # names(result)

    if (return_form == "list") {
      result <- result
    } else{
      version <- result$version
      monisotopic_moleculate_weight <-
        ifelse(
          is.null(result$monisotopic_moleculate_weight),
          NA,
          as.numeric(result$monisotopic_moleculate_weight)
        )
      average_molecular_weight <-
        ifelse(
          is.null(result$average_molecular_weight),
          NA,
          as.numeric(result$average_molecular_weight)
        )
      chemical_formula <-
        ifelse(is.null(result$chemical_formula),
               NA,
               result$average_molecular_weight)
      synonyms <-
        ifelse(is.null(result$synonyms),
               NA,
               paste(unlist(result$synonyms), collapse = "{}"))
      description <-
        ifelse(is.null(result$description),
               NA,
               result$description)
      result <-
        data.frame(
          version = version,
          creation_date = result$creation_date,
          update_date = result$update_date,
          accession = result$accession,
          name = result$name,
          description = description,
          synonyms = synonyms,
          chemical_formula = chemical_formula,
          average_molecular_weight = average_molecular_weight,
          monisotopic_moleculate_weight = monisotopic_moleculate_weight,
          iupac_name = ifelse(is.null(result$iupac_name), NA, result$iupac_name),
          traditional_iupac = ifelse(
            is.null(result$traditional_iupac),
            NA,
            result$traditional_iupac
          ),
          cas_registry_number = ifelse(
            is.null(result$cas_registry_number),
            NA,
            result$cas_registry_number
          ),
          smiles = ifelse(is.null(result$smiles), NA, result$smiles),
          inchi = ifelse(is.null(result$inchi), NA, result$inchi),
          inchikey = ifelse(is.null(result$inchikey), NA, result$inchikey),
          state = ifelse(is.null(result$state), NA, result$state),
          pathways = ifelse(is.null(result$pathways), NA, result$pathways),
          hmdb_id = ifelse(is.null(result$hmdb_id), NA, result$hmdb_id),
          pubchem_compound_id = ifelse(
            is.null(result$pubchem_compound_id),
            NA,
            result$pubchem_compound_id
          ),
          chemspider_id = ifelse(is.null(result$chemspider_id),
                                 NA,
                                 result$chemspider_id),
          kegg_id = ifelse(is.null(result$kegg_id), NA, result$kegg_id),
          chebi_id = ifelse(is.null(result$chebi_id), NA, result$chebi_id),
          biocyc_id = ifelse(is.null(result$biocyc_id), NA, result$biocyc_id),
          het_id = ifelse(is.null(result$het_id), NA, result$het_id),
          wikipidia = ifelse(is.null(result$wikipidia), NA, result$wikipidia),
          vmh_id = ifelse(is.null(result$vmh_id), NA, result$vmh_id),
          fbonto_id = ifelse(is.null(result$fbonto_id), NA, result$fbonto_id),
          hmdb_id = ifelse(is.null(result$hmdb_id), NA, result$hmdb_id),
          # general_references = ifelse(
          #   is.null(result$general_references),
          #   NA,
          #   result$general_references
          # ),
          flavors = ifelse(is.null(result$flavors), NA, result$flavors),
          enzymes = ifelse(is.null(result$enzymes), NA, result$enzymes)
        )
    }
    invisible(result)
  }




#' @title Request one specific reaction information in HMDB
#' @description Request one specific reaction information in HMDB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://hmdb.ca/reactions".
#' @param reaction_id reaction id. For example, 241
#' @return A data frame .
#' @importFrom XML xmlTreeParse xmlToList
#' @importFrom magrittr %>%
#' @importClassesFrom metid databaseClass
#' @export
#' @examples
#' request_hmdb_reaction(reaction_id = "1")

request_hmdb_reaction <-
  function(url = "https://hmdb.ca/reactions",
           reaction_id = "241") {
    final_url <-
      paste0(url, "/", reaction_id)

    result <-
      tryCatch(
        rvest::read_html(x = final_url),
        error = function(e)
          NULL
      )

    if (is.null(result)) {
      message("Check url: ", final_url)
      return(NULL)
    }

    equation <-
      tryCatch(
        result %>%
          rvest::html_elements(".panel-heading") %>%
          rvest::html_text(),
        error = function(e)
          NULL
      )

    main <-
      tryCatch(
        result %>%
          rvest::html_elements(".panel-body") %>%
          rvest::html_text2(),
        error = function(e)
          NULL
      )

    if (!is.null(main)) {
      main <-
        main %>%
        stringr::str_replace_all(pattern = "\\\t", "") %>%
        stringr::str_split("\\\n") %>%
        `[[`(1)
      main <-
        main[!stringr::str_detect(main, "\\=")]

      enzymes_idx <- which(main == "Enzymes")
      External_Links_idx <- which(main == "External Links")
      Status_idx <- which(main == "Status")
      Comments_idx <- which(main == "Comments")

      if (length(enzymes_idx) > 0) {
        Enzymes <- main[enzymes_idx + 1]
      } else{
        Enzymes <- NA
      }

      if (length(External_Links_idx) > 0) {
        External_Links <- main[External_Links_idx + 1]
      } else{
        External_Links <- NA
      }

      if (length(Status_idx) > 0) {
        Status <- main[Status_idx + 1]
      } else{
        Status <- NA
      }

      if (length(Comments_idx) > 0) {
        Comments <- main[Comments_idx + 1]
      } else{
        Comments <- NA
      }

      info <-
        data.frame(Enzymes,
                   External_Links,
                   Status,
                   Comments)
    } else{
      info <-
        data.frame(
          Enzymes = NA,
          External_Links = NA,
          Status = NA,
          Comments = NA
        )
    }

    info <-
      cbind(equation, info)

    info
  }


#' @title Search metabolite in HMDB
#' @description Search metabolite in HMDB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param name metabolite name
#' @param mz mz
#' @param data_database data_database, metID format
#' @param mz_error_cutoff mz_error_cutoff
#' @param similarity_score_cutoff similarity_score_cutoff
#' @return a data frmae
#' @importFrom XML xmlTreeParse xmlToList
#' @importFrom magrittr %>%
#' @export

search_hmdb_database <-
  function(name = "5-methoxysalicylic acid",
           mz = 168.042625,
           data_database,
           mz_error_cutoff = 100,
           similarity_score_cutoff = 0.5) {
    metabolite_database <-
      data_database@spectra.info

    mz <- as.numeric(mz)

    if (is.na(mz)) {
      mz <- 0
    }

    mz_error <-
      abs(metabolite_database$mz - mz) * 10 ^ 6 / ifelse(mz < 40, 40, mz)

    idx <-
      which(mz_error < mz_error_cutoff)

    if (length(idx) == 0) {
      result <-
        data.frame(
          query = "",
          Compound.name = "",
          Synonyms = "",
          similarity_score = "",
          Lab.ID = "",
          HMDB.ID = "",
          idx = "",
          mz_error = 0
        ) %>%
        dplyr::filter(mz_error > 1)
      return(result)
    }

    # metabolite_database$Compound.name[idx]
    similarity_score1 <-
      get_words_similarity(word1 = name,
                           word2 = metabolite_database$Compound.name[idx])

    similarity_score2 <-
      metabolite_database$Synonyms[idx] %>%
      purrr::map(function(x) {
        if (is.na(x)) {
          return(data.frame(Synonyms = NA, score2 = 0))
        }
        x <- stringr::str_split(x, pattern = "\\{\\}")[[1]]
        score <- get_words_similarity(word1 = name, word2 = x)
        temp_idx <- which.max(score)
        data.frame(Synonyms = x[temp_idx], score2 = score[temp_idx])
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()

    similarity_score <-
      data.frame(Compound.name =  metabolite_database$Compound.name[idx],
                 score1 = similarity_score1,
                 similarity_score2)

    similarity_score <-
      similarity_score %>%
      dplyr::rowwise() %>%
      dplyr::mutate(similarity_score = max(c(score1, score2))) %>%
      dplyr::select(-c(score1, score2))

    metabolite_database[idx, , drop = FALSE] %>%
      dplyr::select(Lab.ID, HMDB.ID) %>%
      cbind(similarity_score, .) %>%
      dplyr::mutate(idx, mz_error = mz_error[idx]) %>%
      dplyr::filter(similarity_score > similarity_score_cutoff) %>%
      dplyr::arrange(dplyr::desc(similarity_score)) %>%
      dplyr::mutate(query = name) %>%
      dplyr::select(query, Compound.name, Synonyms, dplyr::everything())
  }


#' @title Convert HMDB data (MS1 or MS2) to metID format database
#' @description Convert HMDB data (MS1 or MS2) to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data for MS1, it is a data.from from read_xml_data_hmdb.
#' @param path Default is .
#' @param threads threads
#' @param ms1_or_ms2 MS1 or MS2 data.
#' @return metid database class
#' @importFrom magrittr %>%
#' @importFrom plyr . dlply
#' @importFrom metid construct_database
#' @export

convert_hmdb2metid <-
  function(data,
           path = ".",
           threads = 5,
           ms1_or_ms2 = c("ms1", "ms2")) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    ms1_or_ms2 <- match.arg(ms1_or_ms2)

    if (ms1_or_ms2 == "ms1") {
      data <-
        data %>%
        dplyr::filter(!is.na(monisotopic_molecular_weight))

      data[which(data == "", arr.ind = TRUE)] <- NA

      data <-
        data %>%
        dplyr::rename(
          mz = average_molecular_weight,
          Create_date = creation_date,
          Updated_date = update_date,
          Lab.ID = accession,
          Compound.name = name,
          Description = description,
          Synonyms = synonyms,
          Formula = chemical_formula,
          IUPAC_name = iupac_name,
          Traditional_IUPAC_name = traditional_iupac,
          CAS.ID = cas_registry_number,
          SMILES.ID = smiles,
          INCHI.ID = inchi,
          INCHIKEY.ID = inchikey,
          Kingdom = kingdom,
          Super_class = super_class,
          Class = class,
          Sub_class = sub_class,
          State = state,
          Biospecimen_locations = biospecimen_locations,
          Cellular_locations = cellular_locations,
          Tissue_locations = tissue_locations,
          CHEMSPIDER.ID = chemspider_id,
          DRUGBANK.ID = drugbank_id,
          FOODB.ID = foodb_id,
          PUBCHEM.ID = pubchem_compound_id,
          CHEBI.ID = chebi_id,
          KEGG.ID = kegg_id,
          BIOCYC.ID = biocyc_id,
          BIGG.ID = bigg_id,
          WIKIPEDIA.ID = wikipedia_id,
          METLIN.ID = metlin_id
        ) %>%
        dplyr::mutate(
          HMDB.ID = Lab.ID,
          RT = NA,
          mz.pos = NA,
          mz.neg = NA,
          Submitter = "HMDB",
          From_human = "Yes"
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

      temp_file <- tempfile()
      dir.create(temp_file, showWarnings = FALSE)

      readr::write_csv(x = data,
                       file = file.path(temp_file, "data.csv"))

      hmdb_ms1 <-
        metid::construct_database(
          path = temp_file,
          version =  as.character(Sys.Date()),
          metabolite.info.name = "data.csv",
          source = "HMDB",
          link = "https://hmdb.ca/",
          creater = "Xiaotao Shen",
          email = "shenxt@stanford.edu",
          rt = FALSE,
          threads = threads
        )

      save(hmdb_ms1, file = file.path(path, "hmdb_ms1"))
      invisible(hmdb_ms1)
    } else{
      remove_idx <-
        data %>%
        lapply(function(x) {
          nrow(x$ms2)
        }) %>%
        unlist() %>%
        `==`(0) %>%
        which()

      if (length(remove_idx) > 0) {
        data <-
          data[-remove_idx]
      }

      spectra_info <-
        data %>%
        purrr::map(function(x) {
          x$ms1_info
        }) %>%
        dplyr::bind_rows() %>%
        as.data.frame()

      spectra_data <-
        data %>%
        purrr::map(function(x) {
          x$ms2
        })

      spectra_info[which(spectra_info == "NA", arr.ind = TRUE)] <-
        NA
      spectra_info[which(spectra_info == "n/a", arr.ind = TRUE)] <-
        NA
      spectra_info[which(spectra_info == "N/A", arr.ind = TRUE)] <-
        NA

      spectra_info <-
        spectra_info %>%
        dplyr::select(HMDB.ID,
                      Instrument_type,
                      Polarity,
                      collision_energy_voltage,
                      adduct)

      remove_idx <-
        which(is.na(spectra_info$Polarity))

      if (length(remove_idx) > 0) {
        spectra_info <-
          spectra_info[-remove_idx, ]

        spectra_data <-
          spectra_data[-remove_idx]
      }

      spectra_info <-
        spectra_info %>%
        dplyr::mutate(
          Polarity = case_when(
            Polarity == "positive" ~ "Positive",
            Polarity == "negative" ~ "Negative",
            TRUE ~ Polarity
          )
        )

      spectra_info$Lab.ID <-
        masstools::name_duplicated(spectra_info$HMDB.ID) %>%
        paste("shen", sep = "_")

      spectra_info2 <-
        spectra_info %>%
        plyr::dlply(.variables = .(HMDB.ID)) %>%
        purrr::map(function(y) {
          if (sum(is.na(y$collision_energy_voltage)) > 0) {
            y$collision_energy_voltage[is.na(y$collision_energy_voltage)] <-
              paste("Unknown",
                    1:length(y$collision_energy_voltage[is.na(y$collision_energy_voltage)]),
                    sep = "_")
          }
          y
        }) %>%
        dplyr::bind_rows() %>%
        as.data.frame()

      spectra_info2 <-
        spectra_info2[match(spectra_info$Lab.ID, spectra_info2$Lab.ID), ]

      spectra_data2 <-
        1:length(spectra_data) %>%
        purrr::map(function(i) {
          x <- spectra_data[[i]]
          x <- list(x)
          names(x) <-
            spectra_info2$collision_energy_voltage[i]
          x
        })

      names(spectra_data2) <- spectra_info2$Lab.ID

      ######positive mode
      spectra_info2$Lab.ID == names(spectra_data2)

      index_pos <- which(spectra_info2$Polarity == "Positive")
      index_neg <- which(spectra_info2$Polarity == "Negative")

      spectra_info_pos <- spectra_info2[index_pos, ]
      spectra_data_pos <- spectra_data2[index_pos]

      spectra_info_neg <- spectra_info2[index_neg, ]
      spectra_data_neg <- spectra_data2[index_neg]

      colnames(spectra_info2)
      colnames(hmdb_ms1@spectra.info)

      spectra_info2 <-
        spectra_info2 %>%
        dplyr::rename(CE = "collision_energy_voltage")

      spectra_info2 <-
        spectra_info2 %>%
        dplyr::left_join(hmdb_ms1@spectra.info %>% dplyr::select(-Lab.ID),
                         by = c("HMDB.ID"))

      temp_file <- tempfile()
      dir.create(temp_file, showWarnings = FALSE)

      readr::write_csv(x = spectra_info2,
                       file = file.path(temp_file, "spectra_info2.csv"))

      hmdb_ms2 <-
        metid::construct_database(
          path = temp_file,
          version =  as.character(Sys.Date()),
          metabolite.info.name = "spectra_info2.csv",
          source = "HMDB",
          link = "https://hmdb.ca/",
          creater = "Xiaotao Shen",
          email = "shenxt@stanford.edu",
          rt = FALSE,
          threads = threads
        )

      hmdb_ms2@spectra.data$Spectra.positive <- spectra_data_pos
      hmdb_ms2@spectra.data$Spectra.negative <- spectra_data_neg

      save(hmdb_ms2, file = file.path(path, "hmdb_ms2"))
      message("Done.")
      invisible(hmdb_ms2)
    }

  }
