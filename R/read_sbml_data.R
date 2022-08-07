# read_sbml_data <-
#   function (file) {
#     edoc <- XML::xmlEventParse(file, handlers = sbml_handler(),
#                                ignoreBlanks = TRUE)
#     model <-
#       edoc$get_model()
#
#     doc <-
#       XML::xmlTreeParse(file, ignoreBlanks = TRUE)
#
#     model$htmlNotes <-
#       doc$doc$children$sbml[["model"]][["notes"]]
#
#     rules <-
#       doc$doc$children$sbml[["model"]][["listOfRules"]]
#
#     reactions <-
#       doc$doc$children$sbml[["model"]][["listOfReactions"]]
#
#     globa_parameters <-
#       names(model$globa_parameters)
#
#     rule_number <- length(rules)
#
#     if (rule_number > 0) {
#       for (i in seq_len(rule_number)) {
#         mathml <- rules[[i]][["math"]][[1]]
#         model$rules[[i]]$mathmlLaw = mathml
#         e <- mathml2R(mathml)
#         model$rules[[i]]$exprLaw <- e[[1]]
#         model$rules[[i]]$strLaw <- gsub(" ", "", toString(e[1]))
#         leaves <- getRuleLeaves(mathml)
#         r <-
#           model$rules[[i]]$inputs <-
#           setdiff(leaves, globa_parameters)
#         model$rules[[i]]$law = make_law(r, NULL, model$rules[[i]]$exprLaw)
#       }
#     }
#
#     reaction_number = length(reactions)
#     if (reaction_number > 0) {
#       for (i in 1:reaction_number) {
#         model$reactions[[i]]$mathmlLaw <-
#           reactions[[i]][["kineticLaw"]][["math"]][[1]]
#         e <-
#           mathml2R(reactions[[i]][["kineticLaw"]][["math"]][[1]])
#         model$reactions[[i]]$exprLaw = e[[1]]
#         model$reactions[[i]]$strLaw = gsub(" ", "", toString(e[1]))
#         r = model$reactions[[i]]$reactants
#         p = names(model$reactions[[i]]$parameters)
#         m = model$reactions[[i]]$modifiers
#         e = model$reactions[[i]]$exprLaw
#         model$reactions[[i]]$law <- make_law(c(r, m), p, e)
#       }
#     }
#     model
#   }
#
#
#
#
#
#
# sbml_handler <- function() {
#   sbml <- "x"
#   modelid <- "x"
#   lnotes <- NULL
#   compartments <- list()
#   reactLaws <- list()
#   species <- list()
#   rules <- list()
#   reactions <- list()
#   globa_parameters = list()
#   reactants = NULL
#   products = NULL
#   modifiers = NULL
#   currRxnID = NULL
#   parameters = NULL
#   parameterIDs = NULL
#   globalParameterIDs = NULL
#   notes = FALSE
#   reactant = FALSE
#   product = FALSE
#   law = FALSE
#   parameter = FALSE
#   math = FALSE
#   .startElement <- function(name, atts, ...) {
#     if (name == "sbml")
#       sbml <<- atts
#     if (name == "annotation")
#       print("skipping annotation")
#     if (name == "model") {
#       modelid <<- atts[["id"]]
#     }
#     if (name == "compartment")
#       if ("id" %in% names(atts))
#         compartments[[atts["id"]]] <<- atts
#     if (name == "species")
#       if ("id" %in% names(atts))
#         species[[atts["id"]]] <<- atts
#     if (name == "assignmentRule")
#       rules[[atts["variable"]]]$idOutput <<- atts[["variable"]]
#     if (name == "reaction") {
#       lstnames <- names(atts)
#       numitems <- length(lstnames)
#       nameslist <- list()
#       id <- "x"
#       reverse <- FALSE
#       name <- "x"
#       count <- 1
#       while (count <= numitems) {
#         switch(
#           lstnames[[count]],
#           id = {
#             id = atts[[count]]
#             nameslist[[length(nameslist) + 1]] <- "id"
#           },
#           reversible = {
#             reverse = as.logical(atts[[count]])
#             nameslist[[length(nameslist) + 1]] <- "reversible"
#           },
#           name = {
#             name = as.character(atts[[count]])
#             nameslist[[length(nameslist) + 1]] <- "name"
#           }
#         )
#         count <- count + 1
#       }
#       reactions[[atts["id"]]]$id <<- id
#       reactions[[atts["id"]]]$reversible <<- reverse
#       currRxnID <<- atts["id"]
#     }
#     if (name == "listOfReactants")
#       reactant <<- TRUE
#     if (name == "listOfProducts")
#       product <<- TRUE
#     if (name == "kineticLaw")
#       law <<- TRUE
#     if (name == "math")
#       math <<- TRUE
#     if ((name == "speciesReference") & reactant)
#       reactants <<- c(reactants, species = atts[["species"]])
#     if ((name == "speciesReference") & product)
#       products <<- c(products, species = atts[["species"]])
#     if (name == "modifierSpeciesReference")
#       modifiers <<- c(modifiers, species = atts[["species"]])
#     if ((name == "parameter") & law) {
#       parameterIDs <<- c(parameterIDs, atts[["id"]])
#       parameters <<- c(parameters, atts[["value"]])
#     }
#     if ((name == "parameter") & (!law)) {
#       globalParameterIDs <<- c(globalParameterIDs,
#                                atts[["id"]])
#       globa_parameters <<-
#         c(globa_parameters, as.numeric(atts[["value"]]))
#     }
#   }
#   .endElement <- function(name) {
#     if (name == "listOfReactants")
#       reactant <<- FALSE
#     if (name == "listOfProducts")
#       product <<- FALSE
#     if (name == "kineticLaw")
#       law <<- FALSE
#     if (name == "math")
#       math <<- FALSE
#     if ((name == "listOfParameters") & (!law))
#       names(globa_parameters) <<- globalParameterIDs
#     if (name == "reaction") {
#       names(reactants) <<- NULL
#       names(modifiers) <<- NULL
#       names(products) <<- NULL
#       reactions[[currRxnID]]$reactants <<- reactants
#       reactions[[currRxnID]]$modifiers <<- modifiers
#       reactions[[currRxnID]]$products <<- products
#       parameters <<- as.numeric(parameters)
#       names(parameters) <<- parameterIDs
#       reactions[[currRxnID]]$parameters <<- parameters
#       reactants <<- NULL
#       products <<- NULL
#       modifiers <<- NULL
#       parameters <<- NULL
#       parameterIDs <<- NULL
#     }
#   }
#   .text <- function(x, ...) {
#     if (!math)
#       lnotes <<- c(lnotes, x)
#   }
#   get_model <- function() {
#     fixComps = function(x) {
#       lstnames <- names(x)
#       count <- 1
#       numit <- length(lstnames)
#       id <- "x"
#       size <- 0
#       name <- "x"
#       nameslist <- list()
#       while (count <= numit) {
#         switch(lstnames[[count]],
#                id = {
#                  id = x[[count]]
#                  nameslist[[length(nameslist) + 1]] <- "id"
#                },
#                size = {
#                  size = as.numeric(x[[count]])
#                  nameslist[[length(nameslist) + 1]] <- "size"
#                },
#                name = {
#                  name = as.character(x[[count]])
#                  nameslist[[length(nameslist) + 1]] <- "name"
#                })
#         count = count + 1
#       }
#       if (numit == 2) {
#         if (id == "x")
#           id <- "default"
#         else if (name == "x")
#           name <- id
#         else if (size == "0")
#           size <- 1
#         lst = list(id, size, name)
#         names(lst) <- c("id", "size", "name")
#         lst
#       }
#       else if (numit == 3) {
#         lst = list(id, size, name)
#         names(lst) <- c("id", "size", "name")
#         lst
#       }
#     }
#     fixSpecies = function(x) {
#       numitems <- length(x)
#       lstnames <- names(x)
#       count <- 1
#       id <- "x"
#       ic <- 0
#       compart <- "def"
#       bc <- FALSE
#       name <- "def"
#       nameslist <- list()
#       while (count <= numitems) {
#         switch(
#           lstnames[[count]],
#           id = {
#             id <- x[[count]]
#             nameslist[[length(nameslist) + 1]] <- "id"
#           },
#           name = {
#             name <- x[[count]]
#             nameslist[[length(nameslist) + 1]] <- "name"
#           },
#           initialConcentration = {
#             ic <- as.numeric(x[[count]])
#             nameslist[[length(nameslist) + 1]] <- "ic"
#           },
#           compartment = {
#             compart <- as.character(x[[count]])
#             nameslist[[length(nameslist) + 1]] <- "compartment"
#           },
#           boundaryCondition = {
#             bc <- as.logical(x[[count]])
#             nameslist[[length(nameslist) + 1]] <- "bc"
#           }
#         )
#         count = count + 1
#       }
#       lst = list(id, as.numeric(ic), compart, as.logical(bc))
#       names(lst) <- c("id", "ic", "compartment", "bc")
#       lst
#     }
#     fixParams = function(x) {
#       numitems <- length(x)
#       lstnames <- names(x)
#       count <- 1
#       id <- "x"
#       value <- 0
#       name <- "def"
#       constant <- FALSE
#       nameslist <- list()
#       while (count <= numitems) {
#         switch(
#           lstnames[[count]],
#           id = {
#             id <- x[[count]]
#             nameslist[[length(nameslist) + 1]] <- "id"
#           },
#           name = {
#             name <- x[[count]]
#             nameslist[[length(nameslist) + 1]] <- "name"
#           },
#           value = {
#             value <- as.numeric(x[[count]])
#             nameslist[[length(nameslist) + 1]] <- "value"
#           },
#           constant = {
#             constant <- as.logical(x[[count]])
#             nameslist[[length(nameslist) + 1]] <- "constant"
#           }
#         )
#         count = count + 1
#       }
#       lst = list(id, as.numeric(value))
#       names(lst) <- c("id", "value")
#       lst
#     }
#     compartments = sapply(compartments, fixComps, simplify = FALSE)
#     species = sapply(species, fixSpecies, simplify = FALSE)
#     list(
#       sbml = sbml,
#       id = modelid[[1]],
#       notes = lnotes,
#       compartments = compartments,
#       species = species,
#       globa_parameters = globa_parameters,
#       rules = rules,
#       reactions = reactions
#     )
#   }
#   list(
#     .startElement = .startElement,
#     .endElement = .endElement,
#     .text = .text,
#     get_model = get_model
#   )
# }
#
#
# mathml2R <-
#   function(children) {
#     UseMethod("mathml2R", children)
#   }
#
#
# mathml2R.XMLDocument <- function(children) {
#   return(mathml2R(children$doc$children))
# }
#
# mathml2R.default <-
#   function(children) {
#     if (is.null(children)) {
#       return(NA)
#     }
#     expr <- expression()
#     n = length(children)
#     for (i in 1:n)
#       expr = c(expr, mathml2R(children[[i]]))
#     if (n > 3) {
#       if (expr[[1]] == "*")
#         expr[[1]] = as.name("prod")
#       if (expr[[1]] == "+")
#         expr[[1]] = as.name("sum")
#     }
#     return(expr)
#   }
# mathml2R.XMLNode <- function(node) {
#   nm <- xmlName(node)
#   if (nm == "power" || nm == "divide" || nm == "times" ||
#       nm == "plus" || nm == "minus") {
#     op <- switch(
#       nm,
#       power = "^",
#       divide = "/",
#       times = "*",
#       plus = "+",
#       minus = "-"
#     )
#     val <- as.name(op)
#   }
#   else if ((nm == "ci") | (nm == "cn")) {
#     if (nm == "ci")
#       val <- as.name(node$children[[1]]$value)
#     if (nm == "cn")
#       val <- as.numeric(node$children[[1]]$value)
#   }
#   else if (nm == "apply") {
#     val <- mathml2R(node$children)
#     mode(val) <- "call"
#   }
#   else {
#     cat("error: nm =", nm, " not in set!\n")
#   }
#   return(as.expression(val))
# }
#
#
# make_law <-
#   function(listofSpecies,
#            listofParams,
#            e,
#            compartments = NULL) {
#     attach(compartments)
#     # takes reactant list of Species, parameter list and rate law R expression e
#     # and makes a reaction rate law function out of them.
#     lawTempl = function(listofSpecies, listofParams = NULL) {
#
#     }
#     i = 2
#     for (j in seq(along = listofParams)) {
#       body(lawTempl)[[i]] <-
#         call("=",
#              as.name(listofParams[j]),
#              call("[", as.name("listofParams"), listofParams[j]))
#       i = i + 1
#     }
#     for (j in seq(along = listofSpecies)) {
#       body(lawTempl)[[i]] <-
#         call("=",
#              as.name(listofSpecies[j]),
#              call("[", as.name("listofSpecies"), listofSpecies[j]))
#       i = i + 1
#     }
#     body(lawTempl)[[i]] <- e
#     lawTempl
#   }
