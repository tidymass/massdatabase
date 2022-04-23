#' @title Convert species to source
#' @description Convert species to source
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param x source
#' @param match_table match_table
#' @return A data.frame
#' @importFrom magrittr %>%
#' @importFrom dplyr case_when everything select filter
#' @importFrom purrr map map2 walk
#' @importFrom crayon green
#' @export

convert_species2source <-
  function(x,
           match_table) {
    if (is.na(x)) {
      return(
        data.frame(
          From_human = "No",
          From_animal = "No",
          From_microbiota = "No",
          From_archaea = "No",
          From_bacteria = "No",
          From_fungi = "No",
          From_food = "No",
          From_plant = "No",
          From_drug = "No",
          From_environment = "No",
          From_eukaryota = "No",
          From_other = "Yes"
        )
      )
    }

    x <-
      stringr::str_split(x, "\\{\\}")[[1]]

    match_table <-
      as.data.frame(match_table)
    x <-
      as.character(match_table[, 2])[match(x, as.character(match_table[, 1]))]
    x <- x[!is.na(x)]

    if (length(x) == 0) {
      return(
        data.frame(
          From_human = "No",
          From_animal = "No",
          From_microbiota = "No",
          From_archaea = "No",
          From_bacteria = "No",
          From_fungi = "No",
          From_food = "No",
          From_plant = "No",
          From_drug = "No",
          From_environment = "No",
          From_eukaryota = "No",
          From_virus = "No",
          From_other = "Yes"
        )
      )
    }

    ###
    From_human = "No"
    From_animal = "No"
    From_microbiota = "No"
    From_archaea = "No"
    From_bacteria = "No"
    From_fungi = "No"
    From_food = "No"
    From_plant = "No"
    From_drug = "No"
    From_environment = "No"
    From_eukaryota = "No"
    From_virus = "No"
    From_other = "Yes"

    if (any(x == "Human")) {
      From_human <- "Yes"
      From_other = "No"
    }

    if (any(x == "Animalia")) {
      From_animal <- "Yes"
      From_other = "No"
    }

    if (any(x == "Plantae") |
        any(x == "Archaeplastida") | any(x == "Viridiplantae")) {
      From_plant <- "Yes"
      From_other = "No"
    }

    if (any(x == "Bacteria")) {
      From_microbiota <- "Yes"
      From_bacteria <- "Yes"
      From_other = "No"
    }

    if (any(x == "Fungi")) {
      From_microbiota <- "Yes"
      From_fungi <- "Yes"
      From_other = "No"
    }

    if (any(x == "Archaea")) {
      From_microbiota <- "Yes"
      From_archaea <- "Yes"
      From_other = "No"
    }


    if (any(x == "Eukaryota")) {
      From_eukaryota <- "Yes"
      From_other = "No"
    }

    if (any(x == "Food")) {
      From_food <- "Yes"
      From_other = "No"
    }

    if (any(x == "Environment")) {
      From_environment <- "Yes"
      From_other = "No"
    }

    if (any(x == "Virus")) {
      From_virus <- "Yes"
      From_other = "No"
    }

    if (any(x == "Drug")) {
      From_drug <- "Yes"
      From_other = "No"
    }


    if (any(x == "Food_plant")) {
      From_food <- "Yes"
      From_plant <- "Yes"
      From_other = "No"
    }

    data.frame(
      From_human = From_human,
      From_animal = From_animal,
      From_microbiota = From_microbiota,
      From_archaea = From_archaea,
      From_bacteria = From_bacteria,
      From_fungi = From_fungi,
      From_food = From_food,
      From_plant = From_plant,
      From_drug = From_drug,
      From_environment = From_environment,
      From_eukaryota = From_eukaryota,
      From_virus = From_virus,
      From_other = From_other
    )
  }



lipid_class_table <-
  data.frame(
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
    url = c(
      "ALL",
      "NAG",
      "BBA",
      "DFA",
      "DLL",
      "DLD",
      "DLB",
      "XPR",
      "EEL",
      "VCA",
      "VCQ",
      "VVA",
      "VVD",
      "VVE",
      "VVF",
      "VVK",
      "GSG",
      "GCG",
      "IIP",
      "OPO",
      "ALA",
      "CLS",
      "TLP",
      "MMA",
      "PGP",
      "PSP",
      "SST",
      "WWA"
    )
  )
