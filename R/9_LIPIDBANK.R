#' @title Download lIPIDBANK database
#' @description Download lIPIDBANK database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://lipidbank.jp/download/".
#' @param lipid_class See here: https://lipidbank.jp/index.html.
#' @param path Default is .
#' @return Downloaded files.
#' @importFrom magrittr %>%
#' @export
download_lipidbank_lipid_class <-
  function(url = "https://lipidbank.jp/download/",
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
           path = ".") {
    lipid_class <- match.arg(lipid_class)

    bre <-
      lipid_class_table$url[match(lipid_class, lipid_class_table$lipid_class)]

    url <- paste0("https://lipidbank.jp/download/",
                  bre, ".xlsx")

    message("Downloading...\n")
    download.file(url = url,
                  destfile = file.path(path, paste0(lipid_class, ".xlsx")))
    message("Done.\n")

  }


#' @title Request Lipidbank database
#' @description Request one specific the metabolite information in HMDB
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param url Default is "https://hmdb.ca/metabolites".
#' @param lipid_class lipid_class
#' @return A data frame or list.
#' @importFrom XML xmlTreeParse xmlToList
#' @importFrom magrittr %>%
#' @export
#' @examples
#' x = request_lipidbank_lipid_class(lipid_class = "Wax")
#' head(x)

request_lipidbank_lipid_class <-
  function(url = "https://lipidbank.jp/",
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
           )) {
    lipid_class <- match.arg(lipid_class)
    bre <-
      lipid_class_table$url[match(lipid_class, lipid_class_table$lipid_class)]

    url <- paste0(url, bre, ".html")

    result <-
      tryCatch(
        rvest::read_html(url) %>%
          rvest::html_table(),
        error = function(e) {
          NULL
        }
      )

    if (is.null(result)) {
      message("Check your internet.")
      return(NULL)
    }

    result <- result[[2]] %>%
      as.data.frame()

    colnames(result) <-
      result[1, ] %>%
      as.character()

    result <-
      result[-1, ]

    result
  }


#' @title Convert lipid bank data (data.frame) to metID format database
#' @description Convert lipid bank data (data.frame) to metID format database
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data data.frame
#' @return Downloaded files.
#' @importFrom magrittr %>%
#' @export

convert_lipidbank2metid <-
  function(data) {
    data <-
      data %>%
      dplyr::select(
        -c(
          Image,
          "MASS SPECTRA",
          "OTHER SPECTRA",
          "NMR SPECTRA",
          "IR SPECTRA",
          "UV SPECTRA",
          "Download cdx file / Mol format file",
          "BOILING POINT",
          "DENSITY",
          "REFRACTIVE INDEX",
          "CHEMICAL SYNTHESIS",
          "OPTICAL ROTATION",
          "SOLUBILITY",
          "CHROMATOGRAM DATA",
          "CHEMICAL SYNTHESIS",
          "METABOLISM",
          "GENETIC INFORMATION",
          "NOTE",
          "REFERENCES",
          "Id",
          "MELTING POINT"
        )
      )

    data <-
      data %>%
      dplyr::rename(
        Lab.ID = `DATA No`,
        Compound.name = NAME,
        Informant = INFORMANT,
        Synonyms = `COMMON NAME`,
        Symbol = SYMBOL,
        Main_class_lipidbank = `Lipid class`,
        Average.mass = `MOL.WT(average)`,
        Formula = FORMULA,
        Biological_activity = `BIOOGICAL ACTIVITY`,
        Source = SOURCE
      ) %>%
      dplyr::mutate(
        LIPIDBANK.ID = Lab.ID,
        CAS.ID = NA,
        HMDB.ID = NA,
        KEGG.ID = NA,
        RT = NA,
        mz = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "LIPIDBANK"
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

    data <-
      data %>%
      dplyr::distinct(Lab.ID, .keep_all = TRUE) %>%
      dplyr::filter(!is.na(Compound.name) & !is.na(Formula))


    data$Compound.name <-
      data$Compound.name %>%
      stringr::str_replace_all(' \\<\\<[ a-zA-Z\\.\\/0-9]{1,40}\\>\\> ', "") %>%
      stringr::str_replace_all('\\<\\<[ a-zA-Z\\.\\/0-9]{1,40}\\>\\>', "") %>%
      stringr::str_replace_all('\\\"', "") %>%
      stringr::str_replace_all('\\.$', "") %>%
      stringr::str_replace_all(' \\/ ', "{}") %>%
      stringr::str_replace_all(' \\/', "{}") %>%
      stringr::str_replace_all('\\/ ', "{}") %>%
      stringr::str_replace_all('\\/', "{}") %>%
      stringr::str_replace_all('\\/ \\/', "{}") %>%
      stringr::str_replace_all('\\{\\}$', "") %>%
      stringr::str_replace_all('^\\{\\}', "") %>%
      stringr::str_trim()

    data$Synonyms <-
      data$Synonyms %>%
      stringr::str_replace_all(' \\<\\<[ a-zA-Z\\.\\/0-9]{1,40}\\>\\> ', "") %>%
      stringr::str_replace_all('\\<\\<[ a-zA-Z\\.\\/0-9]{1,40}\\>\\>', "") %>%
      stringr::str_replace_all('\\\"', "") %>%
      stringr::str_replace_all('\\.$', "") %>%
      stringr::str_replace_all(' \\/ ', "{}") %>%
      stringr::str_replace_all(' \\/', "{}") %>%
      stringr::str_replace_all('\\/ ', "{}") %>%
      stringr::str_replace_all('\\/', "{}") %>%
      stringr::str_replace_all('\\/ \\/', "{}") %>%
      stringr::str_replace_all('\\{\\}$', "") %>%
      stringr::str_replace_all('^\\{\\}', "") %>%
      stringr::str_trim()

    data[which(data == "", arr.ind = TRUE)] <- NA

    data <-
      data %>%
      dplyr::filter(!is.na(Compound.name) | !is.na(Synonyms))

    new_name <-
      1:nrow(data) %>%
      purrr::map(function(i) {
        # cat(i, " ")
        name1 <- data$Compound.name[i]
        name2 <- data$Synonyms[i]

        name <-
          c(name1, name2)

        name <-
          name[!is.na(name)]

        name <-
          name %>%
          stringr::str_split(pattern = "\\{\\}") %>%
          unlist() %>%
          unique()

        name <-
          name[order(nchar(name), name)]

        Synonyms = paste(name, collapse = "{}")
        name3 <-
          name[!stringr::str_detect(name, "from|Lipid")]

        if (length(name3) == 0) {
          Compound.names <- name[1]
        } else{
          Compound.names <- name3[1]
        }

        data.frame(Compound.name = Compound.names,
                   Synonyms = Synonyms)
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    data$Compound.name <- new_name$Compound.name
    data$Synonyms <- new_name$Synonyms

    data <-
      data %>%
      dplyr::distinct(Compound.name, Synonyms, .keep_all = TRUE)

    data <-
      data %>%
      dplyr::filter(!is.na(Formula)) %>%
      dplyr::filter(nchar(Formula) > 1)

    data$mz <-
      data$Formula %>%
      purrr::map(function(x) {
        value <-
          tryCatch(
            Rdisop::getMass(Rdisop::getMolecule(x)),
            error = function(e) {
              NA
            }
          )

      }) %>%
      unlist() %>%
      as.numeric()

    data <-
      data %>%
      dplyr::filter(!is.na(mz))
    data$Average.mass <- as.numeric(data$Average.mass)

    ###species
    Species <-
      data$Source %>%
      purrr::map(function(x) {
        if (is.na(x)) {
          return(NA)
        }

        x <-
          x %>%
          stringr::str_replace_all("from African", "") %>%
          stringr::str_replace_all("from Berber", "") %>%
          stringr::str_replace_all("from Chinese", "") %>%
          stringr::str_replace_all("from India", "") %>%
          stringr::str_replace_all("from Korean", "") %>%
          stringr::str_replace_all("\\[Table [0-9]{1,5}\\]", "") %>%
          stringr::str_replace_all("other references\\:", "") %>%
          stringr::str_replace_all("\\(1\\% of total\\)", "")

        ######if it is separated by ;
        if (stringr::str_detect(x, ";")) {
          x2 <-
            x %>%
            # stringr::str_replace_all(x, "\\(Bacteroides\\)", "") %>%
            stringr::str_replace_all(" and ", " ") %>%
            stringr::str_replace_all("  ", " ") %>%
            stringr::str_replace("\\.$", "") %>%
            stringr::str_replace("\\;$", "") %>%
            stringr::str_replace("^\\;", "") %>%
            stringr::str_replace_all("\\/", "") %>%
            stringr::str_replace_all(",", " ") %>%
            stringr::str_replace_all("  ", " ") %>%
            stringr::str_replace_all(pattern = "\\<\\<[ a-zA-Z0-9\\.\\/]{5,50}\\>\\>", "") %>%
            stringr::str_trim()

          x2 <-
            x2 %>%
            stringr::str_split(";") %>%
            `[[`(1) %>%
            stringr::str_trim()

          x2 <- x2[x2 != ""]
          x2 <-
            x2 %>%
            stringr::str_replace_all("^\\.", "") %>%
            stringr::str_replace_all("\\.$", "") %>%
            stringr::str_trim()
          if (all(unlist(lapply(stringr::str_split(x2, " "), length)) <= 5)) {
            return(paste(x2, collapse = "{}"))
          }

        }

        x <-
          x %>%
          # stringr::str_replace_all(x, "\\(Bacteroides\\)", "") %>%
          stringr::str_replace_all(" and ", " ") %>%
          stringr::str_replace_all("  ", " ") %>%
          stringr::str_replace("\\.$", "") %>%
          stringr::str_replace_all("\\/", "") %>%
          stringr::str_replace_all(",", " ") %>%
          stringr::str_replace_all("  ", " ") %>%
          stringr::str_replace_all(";", " ") %>%
          stringr::str_replace_all("  ", " ") %>%
          stringr::str_trim()

        x2 <-
          stringr::str_split(x, pattern = "\\<\\<[ a-zA-Z0-9\\.\\/]{5,50}\\>\\>")[[1]] %>%
          stringr::str_replace("\\.$", "") %>%
          stringr::str_trim()

        x2 <- x2[x2 != ""]

        x2 <-
          x2 %>%
          lapply(function(y) {
            if (stringr::str_detect(y, "( from ){1}|( of ){1}|( in ){1}")) {
              y <-
                y %>%
                stringr::str_split(pattern = "( from )|( of )|( in )") %>%
                `[[`(1)

              y <- y[-1] %>%
                paste(collapse = " ")

              y <-
                y %>%
                stringr::str_split(" ") %>%
                `[[`(1) %>%
                stringr::str_trim()
              y <-
                y[!y %in% c("in",
                            "of",
                            "the",
                            "many",
                            "from",
                            "and",
                            "during",
                            "with",
                            "a",
                            "high",
                            "or")]
              y <-
                y %>%
                stringr::str_replace_all("^\\.", "") %>%
                stringr::str_replace_all("\\.$", "") %>%
                stringr::str_trim()
              return(paste(y, collapse = " "))

            } else{
              y <-
                y %>%
                stringr::str_replace_all("^\\.", "") %>%
                stringr::str_replace_all("\\.$", "") %>%
                stringr::str_trim()
              return(y)
            }
          }) %>%
          unlist()

        x2 <-
          x2 %>%
          stringr::str_replace_all("^\\.", "") %>%
          stringr::str_replace_all("\\.$", "") %>%
          stringr::str_trim()

        if (all(unlist(lapply(stringr::str_split(x2, " "), length)) <= 5)) {
          return(paste(x2, collapse = "{}"))
        }

        # species1 <-
        #   stringr::str_extract_all(x2, "[a-zA-Z]{3,20} [a-z]{3,20}") %>%
        #   unlist() %>%
        #   unique()
        #
        # species1 <-
        #   species1[!stringr::str_detect(species1, "from|in|of")]

        species2 <-
          stringr::str_extract_all(x2, "[A-Z]{1}\\. [a-z]{3,20}") %>%
          unlist() %>%
          unique()

        species3 <-
          stringr::str_extract_all(x2, "[A-Za-z]{3,20} [a-zA-Z]{1,5}\\.") %>%
          unlist() %>%
          unique()

        # species4 <-
        #   stringr::str_extract_all(
        #     x2,
        #     "Human|human|monkey|fish"
        #   ) %>%
        #   unlist() %>%
        #   unique()

        species5 <-
          x2

        species <-
          c(species2,
            species3,
            species5) %>%
          unique()

        species <-
          species %>%
          stringr::str_replace_all("^\\.", "") %>%
          stringr::str_replace_all("\\.$", "") %>%
          stringr::str_trim()

        species <- species[species != ""]

        if (length(species) == 0) {
          return(NA)
        }
        species <- paste(species, collapse = "{}")
        species
      }) %>%
      unlist()

    # temp <-
    #   data.frame(Source = data$Source,
    #              Species) %>%
    #   dplyr::filter(!is.na(Source) | !is.na(Species))
    #
    # openxlsx::write.xlsx(temp, file = "temp.xlsx", asTable = TRUE)

    data$Species <- Species
    #
    #   all_word <-
    #   Species[!is.na(Species)] %>%
    #     stringr::str_split("\\{\\}") %>%
    #     unlist() %>%
    #     unique()
    #
    #   ##get match table
    #
    #   match_table <-
    #     data.frame(species = stringr::str_to_lower(all_word)) %>%
    #     dplyr::mutate(source = dplyr::case_when(
    #       stringr::str_detect(species, "bacteria|sulfolobus|methanosphaera|pyrococcus") ~ "Bacteria",
    #       stringr::str_detect(species, "mucor|sponge|bacterium|brevundimonas|proteus") ~ "Bacteria",
    #       stringr::str_detect(species, " sp") ~ "Bacteria",
    #       stringr::str_detect(species, "pseudomonas|rhizobium|rhodovulum|rhodospirillum|rhodomicrobium") ~ "Bacteria",
    #       stringr::str_detect(species, "rhodocyclus|salmonella|schizophylum|acinetobacter|actinobacillus") ~ "Bacteria",
    #       stringr::str_detect(species, "aeromonas|aeropyrum|bordetella|campylobacter|thermococcus") ~ "Bacteria",
    #       stringr::str_detect(species, "yersinia|vibrio|vibrio|xanthomonas|ustilago|uredovora|achlya|acetobacter") ~ "Bacteria",
    #       stringr::str_detect(species, "aspergillus|providencia|porphyromonas|streptococcus|strain|bacilli") ~ "Bacteria",
    #       stringr::str_detect(species, "tubercle|anacystis|bacteroides|azospirillum|rhodococcus|rhodobacter|pyramimonas") ~ "Bacteria",
    #       stringr::str_detect(species, "penicillium|pectinatus|pasteurianus|nocardia|neurospora|neisseria|coli") ~ "Bacteria",
    #       stringr::str_detect(species, "fungal|fungus|fungu|schizosaccharomyces|yeast") ~ "Fungi",
    #       stringr::str_detect(species, "human|infant|patient|neonatal|newborn|serum|urine") ~ "Human",
    #       stringr::str_detect(species, "feces|barin|lung|spleen|retina|cortex|female|pregnant|plasma") ~ "Human",
    #       stringr::str_detect(species, "blood|brain|feces|liver|tissue|urinary|adrenal|gestation|testis|skin") ~ "Human",
    #       stringr::str_detect(species, "dog|rat|vertebrates|rabbit|snake|chiken|cattle|chicken") ~ "Animalia",
    #       stringr::str_detect(species, "kangaroos|opossum|pig|koala|mammalian|python") ~ "Animalia",
    #       stringr::str_detect(species, "pelicans|owls|baboon|kite|rana|bullfrog") ~ "Animalia",
    #       stringr::str_detect(species, "bufo|varanus|amyda|fish|moth|drasche|bollworm") ~ "Animalia",
    #       stringr::str_detect(species, "sheep|shark|mouse|animal|paca|moschatus|monkey") ~ "Animalia",
    #       stringr::str_detect(species, "menhaden|boar|bovine|porcine|rodent|rice|mytilus") ~ "Animalia",
    #       stringr::str_detect(species, "prasinophyceae|scallop|quinqueradiata|shigella|worm") ~ "Animalia",
    #       stringr::str_detect(species, "primates|alligator|amphiuma|animal|bird|caiman|whale") ~ "Animalia",
    #       stringr::str_detect(species, "xenopus|bee|grease|wool|latirostris|visceral|viscera|ant") ~ "Animalia",
    #       stringr::str_detect(species, "ascidiacea|marine|asteroidea|clam|sea|porifera|toad|tilapia|thunnus") ~ "Animalia",
    #       stringr::str_detect(species, "atelesto|sturgeon|coral|pernyi|alestes|amoebae|arapaima|potamon|crab|oyster") ~ "Animalia",
    #       stringr::str_detect(species, "ox ") ~ "Animalia",
    #       stringr::str_detect(species, "oncoryhnchus") ~ "Animalia",
    #       stringr::str_detect(species, "sertifer|cyperus|asteraceae|alga|caldariella|plant|codium") ~ "Plantae",
    #       stringr::str_detect(species, "rhamnus|cinchona|cotton|trillium|dioscorea|seed|root|sarsaparilla") ~ "Plantae",
    #       stringr::str_detect(species, "scilla|corn|ruscus|petals|ricinus|sansevieria") ~ "Plantae",
    #       stringr::str_detect(species, "hydnocarpus|palm|acnistus|amaroucium|tunicate|capsicum") ~ "Plantae",
    #       stringr::str_detect(species, "syringae|creeper|pterosperma|gonyaulax|thevetia") ~ "Plantae",
    #       stringr::str_detect(species, "solanaceae|liliaceae|trichosanthes|tragopogon|asclepias|umbellatum") ~ "Plantae",
    #       stringr::str_detect(species, "solanaceae|liliaceae|trichosanthes|tragopogon|asclepias|umbellatum") ~ "Plantae",
    #       stringr::str_detect(species, "thalictrum|leave|tree|leaves|flower|rhizomes|rambutan|racemosa|peucedanum") ~ "Plantae",
    #       stringr::str_detect(species, "neriifolia") ~ "Plantae",
    #       stringr::str_detect(species, "soil|digitonin|diginatin|water|tabacco") ~ "Environment",
    #       stringr::str_detect(species, "pine|ginseng|bean|peach|avocado|carrot|wheat|sultana") ~ "Food_plant",
    #       stringr::str_detect(species, "lemon|cabbage|cocoa|fruit|violet|grape|apple|berries|tanghinin|sunflower") ~ "Food_plant",
    #       stringr::str_detect(species, "apium|banana|peanut|potato|olive") ~ "Food_plant",
    #       stringr::str_detect(species, "milk|lobster|menhaden|marronnier|vinegar|tanghinin|sardine|butter") ~ "Food"
    #     ))

    data("match_table", envir = environment())

    source <-
      data$Species %>%
      stringr::str_to_lower() %>%
      purrr::map(function(x) {
        # cat(x, " ")
        convert_species2source(x = x,
                               match_table = match_table)
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()

    data <-
      cbind(data, source) %>%
      as.data.frame()

    invisible(data)

  }
