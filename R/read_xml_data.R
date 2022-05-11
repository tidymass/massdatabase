#---------------------------------------------------------------------------
#' @title read_xml_data
#' @description Read xml data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The vector of names of files.
#' @param source hmdb
#' @param ms1_or_ms2 ms1 or ms2 data.
#' @param path where is the file(s).
#' @param threads threads
#' @return Return ms2 data. This is a list.
#' @export

read_xml_data <-
  function(file,
           source = c("hmdb"),
           ms1_or_ms2 = c("ms1", "ms2"),
           path = ".",
           threads = 5) {
    source <- match.arg(source)
    xml_data <- readr::read_lines(file)
    if (source == "hmdb") {
      data <-
        read_xml_data_hmdb(
          file = file,
          ms1_or_ms2 = ms1_or_ms2,
          path = path,
          threads = threads
        )
    }

    return(data)

  }

#---------------------------------------------------------------------------
#' @title read_xml_data_hmdb
#' @description Read xml data from HMDB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param file The file name(s).
#' @param threads threads
#' @param path where is/are the file(s)
#' @param ms1_or_ms2 MS1 or MS2 data.
#' @return Return ms2 data. This is a list.
#' @importFrom xml2 as_list read_xml
#' @export

read_xml_data_hmdb <-
  function(file,
           threads = 5,
           path = ".",
           ms1_or_ms2 = c("ms1", "ms2")) {
    ms1_or_ms2 <- match.arg(ms1_or_ms2)
    message("Reading xml data...")
    if (ms1_or_ms2 == "ms1") {
      data <- xml2::read_xml(x = file.path(path, file))
      data <- xml2::as_list(data)
      message("Done.")

      hmdb_metabolite <-
        vector(mode = "list", length = length(data$hmdb))

      progresser <-
        masstools::show_progresser(index = seq_along(data$hmdb))
      message("Extracting information...")
      for (i in seq_along(data$hmdb)) {
        # cat(i, " ")
        if (i %in% progresser$idx) {
          message(progresser$progresser[which(progresser$idx == i)], " ",
                  appendLF = FALSE)
        }
        x <- data$hmdb[[i]]
        biospecimen_locations <-
          paste(stringr::str_trim(unname(
            unlist(x$biological_properties$biospecimen_locations)
          )),
          collapse = "{}")
        cellular_locations <-
          paste(stringr::str_trim(unlist(
            x$biological_properties$cellular_locations
          )),
          collapse = "{}")
        tissue_locations <-
          paste(stringr::str_trim(unlist(
            x$biological_properties$tissue_locations
          )),
          collapse = "{}")

        average_molecular_weight <-
          ifelse(is.null(unlist(x$average_molecular_weight)),
                 NA,
                 as.numeric(unlist(x$average_molecular_weight)))

        monisotopic_molecular_weight <-
          ifelse(is.null(unlist(x$monisotopic_molecular_weight)),
                 NA,
                 as.numeric(unlist(x$monisotopic_molecular_weight)))

        secondary_accessions <-
          paste(stringr::str_trim(unlist(x$secondary_accessions)),
                collapse = "{}")

        synonyms <-
          paste(stringr::str_trim(unname(unlist(x$synonyms))), collapse = "{}")

        hmdb_metabolite[[i]] <-
          data.frame(
            version = ifelse(is.null(unlist(x$version)), NA, unlist(x$version)),
            creation_date = ifelse(is.null(unlist(
              x$creation_date
            )), NA, unlist(x$creation_date)),
            update_date = ifelse(is.null(unlist(
              x$update_date
            )), NA, unlist(x$update_date)),
            accession = ifelse(is.null(unlist(x$accession)), NA, unlist(x$accession)),
            status = ifelse(is.null(unlist(x$status)), NA, unlist(x$status)),
            secondary_accessions = secondary_accessions,
            name = unlist(x$name),
            description = ifelse(is.null(unlist(
              x$description
            )), NA, unlist(x$description)),
            synonyms = synonyms,
            chemical_formula = ifelse(is.null(unlist(
              x$chemical_formula
            )), NA, unlist(x$chemical_formula)),
            average_molecular_weight = average_molecular_weight,
            monisotopic_molecular_weight = monisotopic_molecular_weight,
            iupac_name = ifelse(is.null(unlist(x$iupac_name)), NA, unlist(x$iupac_name)),
            traditional_iupac = ifelse(is.null(unlist(
              x$traditional_iupac
            )), NA, unlist(x$traditional_iupac)),
            cas_registry_number = ifelse(is.null(unlist(
              x$cas_registry_number
            )), NA, unlist(x$cas_registry_number)),
            smiles = ifelse(is.null(unlist(x$smiles)), NA, unlist(x$smiles)),
            inchi = ifelse(is.null(unlist(x$inchi)), NA, unlist(x$inchi)),
            inchikey = ifelse(is.null(unlist(x$inchikey)), NA, unlist(x$inchikey)),
            kingdom = ifelse(is.null(unlist(
              x$taxonomy$kingdom
            )),
            NA, unlist(x$taxonomy$kingdom)),
            super_class = ifelse(is.null(unlist(
              x$taxonomy$super_class
            )),
            NA, unlist(x$taxonomy$super_class)),
            class = ifelse(is.null(unlist(
              x$taxonomy$class
            )),
            NA, unlist(x$taxonomy$class)),
            sub_class = ifelse(is.null(unlist(
              x$taxonomy$sub_class
            )),
            NA, unlist(x$taxonomy$sub_class)),
            state = ifelse(is.null(unlist(x$state)), NA, unlist(x$state)),
            biospecimen_locations = biospecimen_locations,
            cellular_locations = cellular_locations,
            tissue_locations = tissue_locations,
            chemspider_id = ifelse(is.null(unlist(
              x$chemspider_id
            )), NA, unlist(x$chemspider_id)),
            drugbank_id = ifelse(is.null(unlist(
              x$drugbank_id
            )), NA, unlist(x$drugbank_id)),
            foodb_id = ifelse(is.null(unlist(x$foodb_id)), NA, unlist(x$foodb_id)),
            pubchem_compound_id = ifelse(is.null(unlist(
              x$pubchem_compound_id
            )), NA, unlist(x$pubchem_compound_id)),
            chebi_id = ifelse(is.null(unlist(x$chebi_id)), NA, unlist(x$chebi_id)),
            kegg_id = ifelse(is.null(unlist(x$kegg_id)), NA, unlist(x$kegg_id)),
            biocyc_id = ifelse(is.null(unlist(x$biocyc_id)), NA, unlist(x$biocyc_id)),
            bigg_id = ifelse(is.null(unlist(x$bigg_id)), NA, unlist(x$bigg_id)),
            wikipedia_id = ifelse(is.null(unlist(
              x$wikipedia_id
            )), NA, unlist(x$wikipedia_id)),
            metlin_id = ifelse(is.null(unlist(x$metlin_id)), NA, unlist(x$metlin_id))
          )
      }
      message("Done.")

      hmdb_metabolite <-
        hmdb_metabolite %>%
        dplyr::bind_rows() %>%
        as.data.frame()

      hmdb_metabolite
    } else{
      progresser <-
        masstools::show_progresser(index = seq_along(file)) %>%
        dplyr::distinct(idx, .keep_all = TRUE)
      message("Extracting information...")
      hmdb_ms2 <-
        seq_along(file) %>%
        purrr::map(function(i) {
          if (i %in% progresser$idx) {
            message(progresser$progresser[which(progresser$idx == i)], " ",
                    appendLF = FALSE)
          }
          x <- file[i]
          data <-
            read_xml(file.path(path, x)) %>%
            xml2::as_list()
          Instrument_type <-
            unlist(data$`ms-ms`$`instrument-type`)
          Instrument_type <-
            ifelse(is.null(Instrument_type), NA, Instrument_type)
          Polarity <-
            unlist(data$`ms-ms`$`ionization-mode`)
          Polarity <-
            ifelse(is.null(Polarity), NA, Polarity)
          collision_energy_level <-
            unlist(data$`ms-ms`$`collision-energy-level`)
          collision_energy_level <-
            ifelse(is.null(collision_energy_level),
                   NA,
                   collision_energy_level)
          collision_energy_voltage <-
            unlist(data$`ms-ms`$`collision-energy-voltage`)
          collision_energy_voltage <-
            ifelse(is.null(collision_energy_voltage),
                   NA,
                   collision_energy_voltage)

          chromatography_type <-
            unlist(data$`ms-ms`$`chromatography-type`)
          chromatography_type <-
            ifelse(is.null(chromatography_type), NA, chromatography_type)
          analyzer_type <-
            unlist(data$`ms-ms`$`analyzer-type`)
          analyzer_type <-
            ifelse(is.null(analyzer_type), NA, analyzer_type)
          ionization_type <-
            unlist(data$`ms-ms`$`ionization-type`)
          ionization_type <-
            ifelse(is.null(ionization_type), NA, ionization_type)
          charge_type <-
            unlist(data$`ms-ms`$`charge-type`)
          charge_type <-
            ifelse(is.null(charge_type), NA, charge_type)
          adduct <-
            unlist(data$`ms-ms`$adduct)
          adduct <-
            ifelse(is.null(adduct), NA, adduct)
          adduct_type <-
            unlist(data$`ms-ms`$`adduct-type`)
          adduct_type <-
            ifelse(is.null(adduct_type), NA, adduct_type)
          adduct_mass <-
            unlist(data$`ms-ms`$`adduct-mass`)
          adduct_mass <-
            ifelse(is.null(adduct_mass), NA, adduct_mass)
          ms1_info <-
            data.frame(
              HMDB.ID = stringr::str_extract(x, "HMDB[0-9]{7,9}"),
              Instrument_type = Instrument_type,
              Polarity = Polarity,
              collision_energy_level = collision_energy_level,
              collision_energy_voltage = collision_energy_voltage,
              chromatography_type = chromatography_type,
              analyzer_type = analyzer_type,
              ionization_type = ionization_type,
              charge_type = charge_type,
              adduct = adduct,
              adduct_type = adduct_type,
              adduct_mass = adduct_mass
            )

          ms2 <-
            data$`ms-ms`$`ms-ms-peaks`
          if (is.null(ms2)) {
            ms2 <- data.frame()
          } else{
            ms2 <-
              lapply(ms2, function(y) {
                unlist(y)
              }) %>%
              dplyr::bind_rows() %>%
              as.data.frame() %>%
              dplyr::select(`mass-charge`, intensity) %>%
              dplyr::rename(mz = `mass-charge`)
          }
          list(ms1_info = ms1_info,
               ms2 = ms2)
        })
      message("Done.")
      hmdb_ms2
    }
  }
