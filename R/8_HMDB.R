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
      readLines(paste0(url, "/", metabolite_id, ".xml"))
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
        ifelse(is.null(result$monisotopic_moleculate_weight),
               NA,
               as.numeric(result$monisotopic_moleculate_weight))
      average_molecular_weight <-
        ifelse(is.null(result$average_molecular_weight),
               NA,
               as.numeric(result$average_molecular_weight))
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


#' @title Search metabolite in HMDB
#' @description Search metabolite in HMDB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param name metabolite name
#' @param mz mz
#' @param hmdb_metabolite_database hmdb_metabolite_database, metID format
#' @param mz_error_cutoff mz_error_cutoff
#' @param similarity_score_cutoff similarity_score_cutoff
#' @return a data frmae
#' @importFrom XML xmlTreeParse xmlToList
#' @importFrom magrittr %>%
#' @export

search_hmdb_database <-
  function(name = "5-methoxysalicylic acid",
           mz = 168.042625,
           hmdb_metabolite_database,
           mz_error_cutoff = 100,
           similarity_score_cutoff = 0.5) {
    metabolite_database <-
      hmdb_metabolite_database@spectra.info

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
