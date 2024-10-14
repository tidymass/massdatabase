#' #' @title Download the RHEA compound database
#' #' @description Download the RHEA compound database
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param url Default is "https://ftp.expasy.org/databases/rhea/biopax/rhea-biopax.owl.gz".
#' #' @param path Default is ..
#' #' @return Downloaded files.
#' #' @importFrom magrittr %>%
#' #' @export
#' download_rhea_reaction <-
#'   function(url = "https://ftp.expasy.org/databases/rhea/biopax/rhea-biopax.owl.gz",
#'            path = ".") {
#'     path <- file.path(path, "data")
#'     dir.create(path)
#'     message("Download rhea-biopax.owl.gz...\n")
#'     download.file(url = url,
#'                   destfile = file.path(path, "rhea-biopax.owl"))
#'     message("Done.\n")
#'   }
#'
#'
#' #' @title Read the RHEA reaction database from download_rhea_reaction function
#' #' @description Read the RHEA reaction database from download_rhea_reaction function
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param path Default is .. Should be same with download_rhea_reaction function.
#' #' @return A data frame
#' #' @importFrom magrittr %>%
#' #' @importFrom plyr dlply .
#' #' @importFrom readr read_delim
#' #' @importFrom dplyr mutate bind_rows select distinct rename full_join filter
#' #' @importFrom tidyr pivot_wider
#' #' @importFrom purrr map
#' #' @export
#' read_rhea_reaction <-
#'   function(path = ".") {
#'     path <- file.path(path, "data")
#'     owl <-
#'       readLines(file.path(path, "rhea-biopax.owl"))
#'   }
#'
#'
#' ####the code is from
#' ####https://github.com/cran/RbioRXN
#' ####credit should go there
#' #' @title Parse owl format data from RHEA
#' #' @description Parse owl format data from RHEA
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param owl From download_rhea_reaction function
#' #' @return A data.frame
#' #' @export
#' parse_rhea_owl_old <-
#'   function(owl) {
#'     ###reaction index
#'     idx <-
#'       grep('<bp:BiochemicalReaction rdf:about=\"http://rdf.rhea-db.org/',
#'            owl)
#'     pb <-
#'       progress::progress_bar$new(total = length(entry))
#'
#'
#'     result <- cbind(rhea_reaction,
#'                     I(equationWithChebi),
#'                     I(equationParticipant))
#'     result[is.na(result)] <- ''
#'     return(result)
#'   }
#'
#'
#'
#' # parse_rhea_owl_old <-
#' #   function(owl) {
#' #     entry = c(
#' #       'biochemicalReaction',
#' #       'transportReaction',
#' #       'equationWithCommonName',
#' #       'ecNumber',
#' #       'metacyc',
#' #       'kegg',
#' #       'sameParticipant',
#' #       'mapped',
#' #       'formuled',
#' #       'polymerization',
#' #       'chemicallyBalanced',
#' #       'iubmb',
#' #       'status',
#' #       'transport',
#' #       'direction',
#' #       'classOfReactions',
#' #       'closeBiochemicalReaction',
#' #       'closeTransport'
#' #     )
#' #
#' #     # define regular expression
#' #     regexp <- list()
#' #     # regexp[[entry[1]]] <-
#' #     #   '(<bp:biochemicalReaction rdf:about="http://identifiers\\.org/rhea/)(.*)(">)' # biochemical reaction
#' #     regexp[[entry[1]]] <-
#' #       '<bp:BiochemicalReaction rdf:about=\"http://rdf.rhea-db.org/' # biochemical reaction
#' #     # regexp[[entry[2]]] <-
#' #     #   '(<bp:transport.* rdf:about=")(.*)(">)' # transport reaction
#' #     regexp[[entry[2]]] <-
#' #       '<bp:Transport rdf:about=\"http://rdf.rhea-db.org/' # transport reaction
#' #     # regexp[[entry[3]]] <-
#' #     #   '(<bp:NAME .*>)(.*)(</bp:NAME>)' # reaction equation expressed by chemical name
#' #     regexp[[entry[3]]] <-
#' #       '<bp:displayName rdf:datatype \\=' # reaction equation expressed by chemical name
#' #
#' #     regexp[[entry[4]]] <-
#' #       '(<bp:EC-NUMBER .*>)(.*)(</bp:EC-NUMBER>)' # EC number
#' #     regexp[[entry[5]]] <-
#' #       '(<bp:XREF rdf:resource="#METACYC:)(.*)(" />)' # cross-reference to MetaCyc
#' #     regexp[[entry[6]]] <-
#' #       '(<bp:XREF rdf:resource="#KEGG_REACTION:)(.*)(" />)' # cross-reference to KEGG
#' #     regexp[[entry[7]]] <-
#' #       '(<bp:XREF rdf:resource="#rel/../RHEA:)(.*)(" />)' # same participants, different direction
#' #     regexp[[entry[8]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Mapped=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[9]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Formuled=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[10]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Polymerization=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[11]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Chemically balanced=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[12]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:IUBMB=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[13]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Status=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[14]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Transport=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[15]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Direction=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[16]]] <-
#' #       '(<bp:COMMENT rdf:datatype = .+>RHEA:Class of reactions=)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[17]]] <-
#' #       '</bp:biochemicalReaction>' # close tag
#' #     regexp[[entry[18]]] <-
#' #       '</bp:transport' # close tag
#' #
#' #     # pre-calculate condtions (TRUE/FALSE) for fast 'for' loop operation
#' #     condition <- list()
#' #
#' #     pb <-
#' #       progress::progress_bar$new(total = length(entry))
#' #
#' #     message("Matching entry...")
#' #     for (i in entry) {
#' #       pb$tick()
#' #       condition[[i]] <-
#' #         grepl(regexp[[i]], owl)
#' #     }
#' #
#' #     # parsing using regular expressions
#' #     message("parsing ", length(owl), ' lines')
#' #
#' #     rhea_reaction <- data.frame()
#' #     rhea_reaction_row <- list()
#' #
#' #     pb <-
#' #       progress::progress_bar$new(total = length(owl))
#' #
#' #     for (i in 1:length(owl)) {
#' #       pb$tick()
#' #       for (j in 1:length(entry)) {
#' #         # tag start
#' #         if (j %in% 1:2 && condition[[entry[j]]][i]) {
#' #           rhea_reaction_row <- list()
#' #           value <- unlist(strsplit(owl[i], split = '"'))[2]
#' #           value <- tail(unlist(strsplit(value, split = '/')), 1)
#' #           rhea_reaction_row[['rheaId']] <-
#' #             c(rhea_reaction_row[[entry[1]]], value)
#' #           rhea_reaction_row[['reactionType']] <- entry[j]
#' #         } else if (j %in% 3:16 && condition[[entry[j]]][i]) {
#' #           value = sub(regexp[[entry[j]]], '\\2', owl[i])
#' #           rhea_reaction_row[[entry[j]]] <-
#' #             paste(rhea_reaction_row[[entry[j]]], value, sep =
#' #                     ',')
#' #         }
#' #         # tag close
#' #         else if (j %in% 17:18 && condition[[entry[j]]][i]) {
#' #           rhea_reaction_row <-
#' #             as.data.frame(rhea_reaction_row, stringsAsFactors = FALSE)
#' #           rhea_reaction <-
#' #             rbind.fill(rhea_reaction, rhea_reaction_row)
#' #           rhea_reaction_row <- list()
#' #         }
#' #       }
#' #     }
#' #
#' #     # post process
#' #     for (i in names(rhea_reaction)) {
#' #       rhea_reaction[[i]] <- sub('^,', '', rhea_reaction[[i]])
#' #       rhea_reaction[[i]] <- gsub('%2b', '+', rhea_reaction[[i]])
#' #     }
#' #     rhea_reaction <- trim(rhea_reaction)
#' #     rhea_reaction$direction <-
#' #       gsub(' ', '-', rhea_reaction$direction)
#' #
#' #     numberOfBiochecmialReaction <-
#' #       length(unique(rhea_reaction[rhea_reaction$reactionType == 'biochemicalReaction', 'rheaId']))
#' #     numberOfTransportReaction <-
#' #       length(unique(rhea_reaction[rhea_reaction$reactionType == 'transportReaction', 'rheaId']))
#' #
#' #     message('# of biochemical reaction: ',
#' #             numberOfBiochecmialReaction)
#' #     message('# of transport reaction: ',
#' #             numberOfTransportReaction)
#' #
#' #     # replace html code
#' #     rhea_reaction[, 'equationWithCommonName'] <-
#' #       gsub(pattern = "&gt;",
#' #            replacement = ">",
#' #            x = rhea_reaction[, 'equationWithCommonName'])
#' #     rhea_reaction[, 'equationWithCommonName'] <-
#' #       gsub(pattern = "&lt;",
#' #            replacement = "<",
#' #            x = rhea_reaction[, 'equationWithCommonName'])
#' #     rhea_reaction[, 'equationWithCommonName'] <-
#' #       gsub(pattern = "&apos;",
#' #            replacement = "'",
#' #            x = rhea_reaction[, 'equationWithCommonName'])
#' #
#' #     # Sorting
#' #     rhea_reaction <-
#' #       rhea_reaction[order(rhea_reaction[, 'rheaId']), ]
#' #     rownames(rhea_reaction) <-
#' #       1:nrow(rhea_reaction)
#' #
#' #     ##### Parsing equation expressed with ChEBI ID #####
#' #
#' #     message("Parsing equation expressed with ChEBI ID")
#' #
#' #     entry <- c('compoundId', 'name', 'chebiId', 'close')
#' #
#' #     # define regular expression
#' #     regexp <- list()
#' #
#' #     regexp[[entry[1]]] <-
#' #       '(<bp:physicalEntity rdf:about="#compound:)(.+)(">)' # ChEBI ID
#' #     regexp[[entry[2]]] <-
#' #       '(<bp:NAME rdf:datatype .+>)(.*)(</bp:NAME>)' # compound common name used in equation
#' #     regexp[[entry[3]]] <-
#' #       '(<bp:COMMENT rdf:datatype = "http://www.w3.org/2001/XMLSchema#string">.* CHEBI:)(.*)(</bp:COMMENT>)'
#' #     regexp[[entry[4]]] <-
#' #       '</bp:physicalEntity>' # close
#' #
#' #     # pre-calculate condtions (TRUE/FALSE) for fast 'for' loop operation
#' #
#' #     condition <- list()
#' #
#' #     pb <-
#' #       progress::progress_bar$new(total = length(entry))
#' #
#' #     message("Matching entry...")
#' #     for (i in entry) {
#' #       pb$tick()
#' #       condition[[i]] <-
#' #         grepl(regexp[[i]], owl)
#' #     }
#' #
#' #     # parsing using regular expressions
#' #
#' #     message("Parsing owl to create list mapping chemical common name into ChEBI ID")
#' #
#' #     rhea_chemical <-
#' #       data.frame()
#' #     rhea_chemical_row <-
#' #       list()
#' #
#' #     pb <-
#' #       progress::progress_bar$new(total = length(owl))
#' #
#' #     for (i in 1:length(owl)) {
#' #       pb$tick()
#' #       for (j in 1:length(entry)) {
#' #         # Tag start. Please change the number in if statement depending on length of entry
#' #         if (j == 1 && condition[[entry[j]]][i]) {
#' #           rhea_chemical_row <- list()
#' #           value <- unlist(strsplit(owl[i], split = '"'))[2]
#' #           value <- sub('#', '', value)
#' #           rhea_chemical_row[[entry[j]]] <-
#' #             c(rhea_chemical_row[[entry[1]]], value)
#' #         }
#' #         else if (j %in% 2:3 && condition[[entry[j]]][i]) {
#' #           value <- sub(regexp[[entry[j]]], '\\2', owl[i])
#' #           rhea_chemical_row[[entry[j]]] <- value
#' #         }
#' #
#' #         # Tag close
#' #         else if (j == 4 && condition[[entry[j]]][i]) {
#' #           rhea_chemical_row <-
#' #             as.data.frame(rhea_chemical_row, stringsAsFactors = FALSE)
#' #           rhea_chemical <-
#' #             rbind.fill(rhea_chemical, rhea_chemical_row)
#' #           rhea_chemical_row <- list()
#' #         }
#' #       }
#' #     }
#' #
#' #     rhea_chemical <- trim(rhea_chemical)
#' #     rhea_chemical[, 'name'] <-
#' #       gsub(pattern = "&gt;",
#' #            replacement = ">",
#' #            x = rhea_chemical[, 'name'])
#' #     rhea_chemical[, 'name'] <-
#' #       gsub(pattern = "&lt;",
#' #            replacement = "<",
#' #            x = rhea_chemical[, 'name'])
#' #     rhea_chemical[, 'name'] <-
#' #       gsub(pattern = "&apos;",
#' #            replacement = "'",
#' #            x = rhea_chemical[, 'name'])
#' #
#' #     numberOfChemicals <-
#' #       nrow(rhea_chemical)
#' #
#' #     message("# of chemicals used in Rhea: ", numberOfChemicals)
#' #
#' #     ##### Building equation with the list #####
#' #
#' #     message("Building equation with the list...")
#' #
#' #     # Define direction
#' #
#' #     directionList <- list()
#' #     directionList[['undefined']] <- ' <?> '
#' #     directionList[['bidirectional']] <- ' <=> '
#' #     directionList[['right-to-left']] <- ' => '
#' #     directionList[['left-to-right']] <- ' => '
#' #
#' #     # Define regular expressions
#' #
#' #     regexp_coefficient <- "(^\\(?[0-9]*n?\\+?[0-9]*\\)? )(.+)"
#' #     regexp_localization <- "(.+)(\\([inout]+\\))"
#' #
#' #     # Conversion (split - conversion - paste)
#' #
#' #     equationWithChebi <- character(nrow(rhea_reaction))
#' #     equationParticipant <- character(nrow(rhea_reaction))
#' #     for (i in 1:nrow(rhea_reaction)) {
#' #       reactantsChebi <- c()
#' #       productsChebi <- c()
#' #       arrow <- directionList[[rhea_reaction[i, 'direction']]]
#' #       arrow2 <- sub('\\?', '\\\\?', arrow)
#' #       participants <-
#' #         unlist(strsplit(rhea_reaction[i, 'equationWithCommonName'],
#' #                         split = arrow2))
#' #       reactants <-
#' #         unlist(strsplit(participants[1], split = " \\+ "))
#' #       reactantsWOcoefficient <-
#' #         sub(regexp_coefficient, "\\2", reactants)
#' #       for (j in 1:length(reactants)) {
#' #         if (grepl(regexp_coefficient, reactants[j])) {
#' #           # process coefficeint
#' #           coefficient <-
#' #             sub(regexp_coefficient, "\\1", reactants[j])
#' #         } else {
#' #           coefficient <- ""
#' #         }
#' #         if (grepl(regexp_localization, reactants[j])) {
#' #           # process (in), (out)
#' #           localization <-
#' #             sub(regexp_localization, "\\2", x = reactants[j])
#' #         } else {
#' #           localization <- ""
#' #         }
#' #         reactantsWOcoefficient[j] <-
#' #           gsub("\\(in\\)", "", reactantsWOcoefficient[j])
#' #         reactantsWOcoefficient[j] <-
#' #           gsub("\\(out\\)", "", reactantsWOcoefficient[j])
#' #         chebi <-
#' #           rhea_chemical[rhea_chemical$name == reactantsWOcoefficient[j], 'chebiId']
#' #         if (length(chebi) > 1) {
#' #           reactants[j] <- "Unknown"
#' #         } else {
#' #           reactants[j] <- paste(coefficient, chebi, localization, sep = "")
#' #         }
#' #         reactantsChebi <- c(reactantsChebi, chebi)
#' #       }
#' #       products <-
#' #         unlist(strsplit(participants[2], split = " \\+ "))
#' #       productsWOcoefficient <-
#' #         sub(regexp_coefficient, "\\2", products)
#' #
#' #       for (k in 1:length(products)) {
#' #         if (grepl(regexp_coefficient, products[k])) {
#' #           # process coefficeint
#' #           coefficient <-
#' #             sub(regexp_coefficient, "\\1", products[k])
#' #         } else {
#' #           coefficient <- ""
#' #         }
#' #         if (grepl(regexp_localization, products[k])) {
#' #           # process (in), (out)
#' #           localization <-
#' #             sub(regexp_localization, "\\2", x = products[k])
#' #         } else {
#' #           localization <- ""
#' #         }
#' #         productsWOcoefficient[k] <-
#' #           gsub("\\(in\\)", "", productsWOcoefficient[k])
#' #         productsWOcoefficient[k] <-
#' #           gsub("\\(out\\)", "", productsWOcoefficient[k])
#' #         chebi <-
#' #           rhea_chemical[rhea_chemical$name == productsWOcoefficient[k], 'chebiId']
#' #         if (length(chebi) > 1) {
#' #           products[k] <- "Unknown"
#' #         } else {
#' #           products[k] <- paste(coefficient, chebi, localization, sep = '')
#' #         }
#' #         productsChebi <- c(productsChebi, chebi)
#' #       }
#' #       equationWithChebi[i] <-
#' #         paste(paste(reactants, collapse = " + "),
#' #               paste(products, collapse = " + "),
#' #               sep = arrow)
#' #       equationParticipant[i] <-
#' #         paste(c(reactantsChebi, productsChebi), collapse =
#' #                 ',')
#' #     }
#' #     result <- cbind(rhea_reaction,
#' #                     I(equationWithChebi),
#' #                     I(equationParticipant))
#' #     result[is.na(result)] <- ''
#' #     return(result)
#' #   }
