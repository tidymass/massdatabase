#' Retrieve Metabolite Information from Wikidata
#'
#' This function queries Wikidata for detailed information about a given metabolite using its Wikidata ID.
#'
#' @param metabolite_id A character string specifying the Wikidata ID of the metabolite (default: "Q18216").
#' @return A list containing detailed metabolite information, including:
#'   \itemize{
#'     \item \code{metabolite_id} - The Wikidata ID of the metabolite.
#'     \item \code{metabolite_name} - The name of the metabolite.
#'     \item \code{metabolite_description} - A brief description of the metabolite.
#'     \item \code{aliases} - Alternative names for the metabolite.
#'     \item \code{wikipedia_url} - The Wikipedia URL for the metabolite.
#'     \item \code{chemical_formula} - The chemical formula of the metabolite.
#'     \item \code{smiles} - The SMILES notation of the metabolite.
#'     \item \code{inchi} - The InChI representation of the metabolite.
#'     \item \code{inchikey} - The InChIKey of the metabolite.
#'     \item \code{cas_id} - The CAS Registry Number.
#'     \item \code{pubchem_cid} - The PubChem Compound ID.
#'     \item \code{chebi_id} - The ChEBI ID.
#'     \item \code{chemspider_id} - The ChemSpider ID.
#'     \item \code{chembl_id} - The ChEMBL ID.
#'     \item \code{kegg_id} - The KEGG Compound ID.
#'     \item \code{hmdb_id} - The Human Metabolome Database (HMDB) ID.
#'     \item \code{lipidmaps_id} - The LIPID MAPS ID.
#'     \item \code{drugbank_id} - The DrugBank ID.
#'     \item \code{rxnorm_id} - The RxNorm ID.
#'     \item \code{ndf_rt_id} - The NDF-RT ID.
#'     \item \code{mesh_descriptor_id} - The MeSH Descriptor ID.
#'     \item \code{unii} - The Unique Ingredient Identifier (UNII).
#'     \item \code{guide_to_pharmacology_ligand_id} - The Guide to Pharmacology Ligand ID.
#'     \item \code{atc_code} - The Anatomical Therapeutic Chemical (ATC) classification code.
#'     \item \code{echa_substance_infocard_id} - The ECHA Substance Infocard ID.
#'     \item \code{chembl_id} - The ChEMBL ID.
#'     \item \code{freebase_id} - The Freebase ID.
#'     \item \code{reaxys} - The Reaxys Registry Number.
#'     \item \code{zvg_number} - The ZVG number.
#'     \item \code{probes_and_drugs_id} - The Probes and Drugs ID.
#'     \item \code{australian_register_of_therapeutic_goods_id} - The Australian Register of Therapeutic Goods ID.
#'     \item \code{other_identifiers} - Various other database identifiers (e.g., SPLASH, DSSTox, PDB Ligand ID, etc.).
#'   }
#' @export
#' @examples
#' request_wikidata_metabolite("Q18216") # Fetches metabolite data for a given Wikidata ID.

request_wikidata_metabolite <-
  function(metabolite_id = "Q18216") {
    url <-
      paste0("https://www.wikidata.org/wiki/Special:EntityData/",
             metabolite_id,
             ".json")

    # Read the JSON file
    result <-
      tryCatch(
        expr = jsonlite::fromJSON(url),
        error = function(e) {
          message("Error: ", e)
          return(NULL)
        }
      )

    result <-
      result$entities[[metabolite_id]]

    # Extract metabolite information
    metabolite_id <- result$id
    if (is.null(metabolite_id)) {
      metabolite_name <- NA
    }
    metabolite_name <- result$labels$en$value
    if (is.null(metabolite_name)) {
      metabolite_name <- NA
    }
    metabolite_description <- result$descriptions$en$value
    if (is.null(metabolite_description)) {
      metabolite_description <- NA
    }
    aliases <-
      result$aliases$en$value
    if (is.null(aliases)) {
      aliases <- NA
    }

    wikipedia_url <-
      result$sitelinks$enwiki$url
    if (is.null(wikipedia_url)) {
      wikipedia_url <- NA
    }

    ###chemical formula
    chemical_formula <-
      result$claims$P274$mainsnak$datavalue$value
    if (is.null(chemical_formula)) {
      chemical_formula <- NA
    }

    ###SMILES
    smiles <-
      result$claims$P233$mainsnak$datavalue$value
    if (is.null(smiles)) {
      smiles <- NA
    }

    ###Identifiers
    ###GND ID
    gnd_id <-
      result$claims$P227$mainsnak$datavalue$value
    if (is.null(gnd_id)) {
      gnd_id <- NA
    }

    ###J9U ID
    j9u_id <-
      result$claims$P8189$mainsnak$datavalue$value
    if (is.null(j9u_id)) {
      j9u_id <- NA
    }

    ##Library of Congress authority ID
    library_of_congress_authority_id <-
      result$claims$P244$mainsnak$datavalue$value
    if (is.null(library_of_congress_authority_id)) {
      library_of_congress_authority_id <- NA
    }

    ###NDL Authority ID
    ndl_authority_id <-
      result$claims$P349$mainsnak$datavalue$value
    if (is.null(ndl_authority_id)) {
      ndl_authority_id <- NA
    }

    ###InChI
    inchi <-
      result$claims$P234$mainsnak$datavalue$value
    if (is.null(inchi)) {
      inchi <- NA
    }

    ###InChIKey
    inchikey <-
      result$claims$P235$mainsnak$datavalue$value
    if (is.null(inchikey)) {
      inchikey <- NA
    }

    ###CAS Registry Number
    cas_id <-
      result$claims$P231$mainsnak$datavalue$value
    if (is.null(cas_id)) {
      cas_id <- NA
    }

    ###ChemSpider ID
    chemspider_id <-
      result$claims$P661$mainsnak$datavalue$value
    if (is.null(chemspider_id)) {
      chemspider_id <- NA
    }

    ##PubChem CID
    pubchem_cid <-
      result$claims$P662$mainsnak$datavalue$value
    if (is.null(pubchem_cid)) {
      pubchem_cid <- NA
    }

    ###Reaxys registry number
    reaxys <-
      result$claims$P1579$mainsnak$datavalue$value
    if (is.null(reaxys)) {
      reaxys <- NA
    }

    ##Gmelin number
    gmelin_number <-
      result$claims$P1578$mainsnak$datavalue$value
    if (is.null(gmelin_number)) {
      gmelin_number <- NA
    }

    ####ChEBI ID
    chebi_id <-
      result$claims$P683$mainsnak$datavalue$value
    if (is.null(chebi_id)) {
      chebi_id <- NA
    }

    ###ChEMBL ID
    chembl_id <-
      result$claims$P592$mainsnak$datavalue$value
    if (is.null(chembl_id)) {
      chembl_id <- NA
    }


    ###UniChem compound ID
    unichem_id <-
      result$claims$P11089$mainsnak$datavalue$value
    if (is.null(unichem_id)) {
      unichem_id <- NA
    }

    ###MassBank
    massbank_id <-
      result$claims$P6689$mainsnak$datavalue$value
    if (is.null(massbank_id)) {
      massbank_id <- NA
    }

    ##SPLASH
    splash <-
      result$claims$P4964$mainsnak$datavalue$value
    if (is.null(splash)) {
      splash <- NA
    }

    ###ICSC
    icsc_id <-
      result$claims$P5220$mainsnak$datavalue$value

    if (is.null(icsc_id)) {
      icsc_id <- NA
    }

    ####CAMEO Chemicals ID
    cameo_chemicals_id <-
      result$claims$P11931$mainsnak$datavalue$value

    if (is.null(cameo_chemicals_id)) {
      cameo_chemicals_id <- NA
    }

    ###RTECS number
    rtecs_number <-
      result$claims$P657$mainsnak$datavalue$value
    if (is.null(rtecs_number)) {
      rtecs_number <- NA
    }

    ###ZVG number
    zvg_number <-
      result$claims$P679$mainsnak$datavalue$value
    if (is.null(zvg_number)) {
      zvg_number <- NA
    }

    ###HSDB ID
    hsdb_id <-
      result$claims$P2062$mainsnak$datavalue$value
    if (is.null(hsdb_id)) {
      hsdb_id <- NA
    }

    ###DSSTox substance ID
    dsstox_substance_id <-
      result$claims$P3117$mainsnak$datavalue$value
    if (is.null(dsstox_substance_id)) {
      dsstox_substance_id <- NA
    }

    ###DSSTOX compound identifier
    dsstox_compound_id <-
      result$claims$P8494$mainsnak$datavalue$value
    if (is.null(dsstox_compound_id)) {
      dsstox_compound_id <- NA
    }

    ###NSC number
    nsc_number <-
      result$claims$P2840$mainsnak$datavalue$value

    if (is.null(nsc_number)) {
      nsc_number <- NA
    }

    ###EC number
    ec_number <-
      result$claims$P232$mainsnak$datavalue$value
    if (is.null(ec_number)) {
      ec_number <- NA
    }

    ###ECHA Substance Infocard ID
    echa_substance_infocard_id <-
      result$claims$P2566$mainsnak$datavalue$value
    if (is.null(echa_substance_infocard_id)) {
      echa_substance_infocard_id <- NA
    }

    ####FL number
    fl_number <-
      result$claims$P9066$mainsnak$datavalue$value
    if (is.null(fl_number)) {
      fl_number <- NA
    }

    ###CosIng number
    cosing_number <-
      result$claims$P3073$mainsnak$datavalue$value
    if (is.null(cosing_number)) {
      cosing_number <- NA
    }

    ###Cannabis Database ID
    cannabis_database_id <-
      result$claims$P11160$mainsnak$datavalue$value
    if (is.null(cannabis_database_id)) {
      cannabis_database_id <- NA
    }

    ####Nikkaji ID
    nikkaji_id <-
      result$claims$P2085$mainsnak$datavalue$value
    if (is.null(nikkaji_id)) {
      nikkaji_id <- NA
    }

    ####Human Metabolome Database ID
    hmdb_id <-
      result$claims$P2057$mainsnak$datavalue$value
    if (is.null(hmdb_id)) {
      hmdb_id <- NA
    }

    ####KNApSAcK ID
    knapsack_id <-
      result$claims$P2064$mainsnak$datavalue$value
    if (is.null(knapsack_id)) {
      knapsack_id <- NA
    }

    ####UNII
    unii <-
      result$claims$P652$mainsnak$datavalue$value

    if (is.null(unii)) {
      unii <- NA
    }

    ####Probes And Drugs ID
    probes_and_drugs_id <-
      result$claims$P11199$mainsnak$datavalue$value
    if (is.null(probes_and_drugs_id)) {
      probes_and_drugs_id <- NA
    }

    ###KEGG ID
    kegg_id <-
      result$claims$P665$mainsnak$datavalue$value
    if (is.null(kegg_id)) {
      kegg_id <- NA
    }

    ###LIPID MAPS ID
    lipidmaps_id <-
      result$claims$P2063$mainsnak$datavalue$value
    if (is.null(lipidmaps_id)) {
      lipidmaps_id <- NA
    }

    ###Brockhaus Enzyklopädie online ID
    brockhaus_enzyklopadie_online_id <-
      result$claims$P5019$mainsnak$datavalue$value

    if (is.null(brockhaus_enzyklopadie_online_id)) {
      brockhaus_enzyklopadie_online_id <- NA
    }

    ###CA PROP 65 ID
    ca_prop_65_id <-
      result$claims$P7524$mainsnak$datavalue$value
    if (is.null(ca_prop_65_id)) {
      ca_prop_65_id <- NA
    }

    ###Encyclopædia Britannica Online ID
    encyclopaedia_britannica_online_id <-
      result$claims$P1417$mainsnak$datavalue$value
    if (is.null(encyclopaedia_britannica_online_id)) {
      encyclopaedia_britannica_online_id <- NA
    }

    ###Freebase ID
    freebase_id <-
      result$claims$P646$mainsnak$datavalue$value
    if (is.null(freebase_id)) {
      freebase_id <- NA
    }

    ###HCIS ID
    hcis_id <-
      result$claims$P7025$mainsnak$datavalue$value
    if (is.null(hcis_id)) {
      hcis_id <- NA
    }

    ####J-GLOBAL ID
    j_global_id <-
      result$claims$P7783$mainsnak$datavalue$value
    if (is.null(j_global_id)) {
      j_global_id <- NA
    }

    ###JECFA number
    jecfa_number <-
      result$claims$P9557$mainsnak$datavalue$value
    if (is.null(jecfa_number)) {
      jecfa_number <- NA
    }

    ###Microsoft Academic ID
    microsoft_academic_id <-
      result$claims$P6366$mainsnak$datavalue$value

    if (is.null(microsoft_academic_id)) {
      microsoft_academic_id <- NA
    }

    ###OpenAlex ID
    openalex_id <-
      result$claims$P10283$mainsnak$datavalue$value
    if (is.null(openalex_id)) {
      openalex_id <- NA
    }

    ####Quora topic ID
    quora_topic_id <-
      result$claims$P3417$mainsnak$datavalue$value
    if (is.null(quora_topic_id)) {
      quora_topic_id <- NA
    }

    ###ScienceDirect topic ID
    sciencedirect_topic_id <-
      result$claims$P10376$mainsnak$datavalue$value

    if (is.null(sciencedirect_topic_id)) {
      sciencedirect_topic_id <- NA
    }

    ###SureChEMBL ID
    surechembl_id <-
      result$claims$P2877$mainsnak$datavalue$value
    if (is.null(surechembl_id)) {
      surechembl_id <- NA
    }

    ###CCDC Number
    ccdc_number <-
      result$claims$P6852$mainsnak$datavalue$value
    if (is.null(ccdc_number)) {
      ccdc_number <- NA
    }

    ###CSD Refcode
    csd_refcode <-
      result$claims$P11375$mainsnak$datavalue$value
    if (is.null(csd_refcode)) {
      csd_refcode <- NA
    }

    ####ICSC ID
    icsc_id <-
      result$claims$P5220$mainsnak$datavalue$value
    if (is.null(icsc_id)) {
      icsc_id <- NA
    }

    ###LiverTox ID
    livertox_id <-
      result$claims$P7830$mainsnak$datavalue$value
    if (is.null(livertox_id)) {
      livertox_id <- NA
    }

    ###OSHA Occupational Chemical Database ID
    osha_occupational_chemical_database_id <-
      result$claims$P12594$mainsnak$datavalue$value
    if (is.null(osha_occupational_chemical_database_id)) {
      osha_occupational_chemical_database_id <- NA
    }

    ###PesticideInfo chemical ID
    pesticideinfo_chemical_id <-
      result$claims$P11949$mainsnak$datavalue$value
    if (is.null(pesticideinfo_chemical_id)) {
      pesticideinfo_chemical_id <- NA
    }

    ###ATC code
    atc_code <-
      result$claims$P267$mainsnak$datavalue$value
    if (is.null(atc_code)) {
      atc_code <- NA
    }

    ###World Health Organisation international non-proprietary name numeric ID
    who_inn_numeric_id <-
      result$claims$P3350$mainsnak$datavalue$value
    if (is.null(who_inn_numeric_id)) {
      who_inn_numeric_id <- NA
    }

    ###Guide to Pharmacology Ligand ID
    guide_to_pharmacology_ligand_id <-
      result$claims$P595$mainsnak$datavalue$value
    if (is.null(guide_to_pharmacology_ligand_id)) {
      guide_to_pharmacology_ligand_id <- NA
    }

    ###MeSH descriptor ID
    mesh_descriptor_id <-
      result$claims$P486$mainsnak$datavalue$value
    if (is.null(mesh_descriptor_id)) {
      mesh_descriptor_id <- NA
    }

    ###MeSH tree code
    mesh_tree_code <-
      result$claims$P672$mainsnak$datavalue$value
    if (is.null(mesh_tree_code)) {
      mesh_tree_code <- NA
    }

    ###RxNorm ID
    rxnorm_id <-
      result$claims$P3345$mainsnak$datavalue$value
    if (is.null(rxnorm_id)) {
      rxnorm_id <- NA
    }

    ###DrugBank ID
    drugbank_id <-
      result$claims$P715$mainsnak$datavalue$value
    if (is.null(drugbank_id)) {
      drugbank_id <- NA
    }

    ###MedlinePlus ID
    medlineplus_id <-
      result$claims$P604$mainsnak$datavalue$value
    if (is.null(medlineplus_id)) {
      medlineplus_id <- NA
    }

    ###MedlinePlus drug identifier
    medlineplus_drug_identifier <-
      result$claims$P10245$mainsnak$datavalue$value
    if (is.null(medlineplus_drug_identifier)) {
      medlineplus_drug_identifier <- NA
    }

    ###PatientsLikeMe treatment ID
    patientslikeme_treatment_id <-
      result$claims$P4235$mainsnak$datavalue$value
    if (is.null(patientslikeme_treatment_id)) {
      patientslikeme_treatment_id <- NA
    }

    ###NDF-RT ID
    ndf_rt_id <-
      result$claims$P2115$mainsnak$datavalue$value
    if (is.null(ndf_rt_id)) {
      ndf_rt_id <- NA
    }

    ##Australian Register of Therapeutic Goods ID
    australian_register_of_therapeutic_goods_id <-
      result$claims$P3550$mainsnak$datavalue$value
    if (is.null(australian_register_of_therapeutic_goods_id)) {
      australian_register_of_therapeutic_goods_id <- NA
    }

    ###DrugCentral ID
    drugcentral_id <-
      result$claims$P11198$mainsnak$datavalue$value
    if (is.null(drugcentral_id)) {
      drugcentral_id <- NA
    }

    ###IEDB Epitope ID
    iedb_epitope_id <-
      result$claims$P4168$mainsnak$datavalue$value
    if (is.null(iedb_epitope_id)) {
      iedb_epitope_id <- NA
    }

    ###PDB ligand ID
    pdb_ligand_id <-
      result$claims$P3636$mainsnak$datavalue$value
    if (is.null(pdb_ligand_id)) {
      pdb_ligand_id <- NA
    }

    ###BBC Things ID
    bbc_things_id <-
      result$claims$P1617$mainsnak$datavalue$value
    if (is.null(bbc_things_id)) {
      bbc_things_id <- NA
    }

    ###BNCF Thesaurus ID
    bncf_thesaurus_id <-
      result$claims$P508$mainsnak$datavalue$value
    if (is.null(bncf_thesaurus_id)) {
      bncf_thesaurus_id <- NA
    }

    ###DeCS ID
    decs_id <-
      result$claims$P9272$mainsnak$datavalue$value
    if (is.null(decs_id)) {
      decs_id <- NA
    }

    ####Encyclopædia Universalis ID
    encyclopaedia_universalis_id <-
      result$claims$P3219$mainsnak$datavalue$value
    if (is.null(encyclopaedia_universalis_id)) {
      encyclopaedia_universalis_id <- NA
    }

    ###Encyclopedia of China (Third Edition) ID
    encyclopedia_of_china_third_edition_id <-
      result$claims$P10565$mainsnak$datavalue$value
    if (is.null(encyclopedia_of_china_third_edition_id)) {
      encyclopedia_of_china_third_edition_id <- NA
    }

    ###Gran Enciclopèdia Catalana ID
    gran_enciclopedia_catalana_id <-
      result$claims$P12385$mainsnak$datavalue$value
    if (is.null(gran_enciclopedia_catalana_id)) {
      gran_enciclopedia_catalana_id <- NA
    }

    ###LEM ID
    lem_id <-
      result$claims$P902$mainsnak$datavalue$value
    if (is.null(lem_id)) {
      lem_id <- NA
    }

    ###Lex ID
    lex_id <-
      result$claims$P8313$mainsnak$datavalue$value
    if (is.null(lex_id)) {
      lex_id <- NA
    }

    ###NALT ID
    nalt_id <-
      result$claims$P2004$mainsnak$datavalue$value
    if (is.null(nalt_id)) {
      nalt_id <- NA
    }

    ###NE.se ID
    ne_se_id <-
      result$claims$P3222$mainsnak$datavalue$value
    if (is.null(ne_se_id)) {
      ne_se_id <- NA
    }

    ###NHS Health A to Z ID
    nhs_health_a_to_z_id <-
      result$claims$P7995$mainsnak$datavalue$value
    if (is.null(nhs_health_a_to_z_id)) {
      nhs_health_a_to_z_id <- NA
    }

    ###OmegaWiki Defined Meaning
    omegawiki_defined_meaning <-
      result$claims$P1245$mainsnak$datavalue$value
    if (is.null(omegawiki_defined_meaning)) {
      omegawiki_defined_meaning <- NA
    }

    ##Online PWN Encyclopedia ID
    online_pwn_encyclopedia_id <-
      result$claims$P7305$mainsnak$datavalue$value
    if (is.null(online_pwn_encyclopedia_id)) {
      online_pwn_encyclopedia_id <- NA
    }

    ####PDB structure ID
    pdb_structure_id <-
      result$claims$P638$mainsnak$datavalue$value

    if (is.null(pdb_structure_id)) {
      pdb_structure_id <- NA
    }

    return_result <-
      list(
        metabolite_id = metabolite_id,
        metabolite_name = metabolite_name,
        metabolite_description = metabolite_description,
        aliases = aliases,
        wikipedia_url = wikipedia_url,
        chemical_formula = chemical_formula,
        smiles = smiles,
        gnd_id = gnd_id,
        j9u_id = j9u_id,
        library_of_congress_authority_id = library_of_congress_authority_id,
        ndl_authority_id = ndl_authority_id,
        inchi = inchi,
        inchikey = inchikey,
        cas_id = cas_id,
        chemspider_id = chemspider_id,
        pubchem_cid = pubchem_cid,
        reaxys = reaxys,
        gmelin_number = gmelin_number,
        chebi_id = chebi_id,
        chembl_id = chembl_id,
        unichem_id = unichem_id,
        massbank_id = massbank_id,
        splash = splash,
        icsc_id = icsc_id,
        cameo_chemicals_id = cameo_chemicals_id,
        rtecs_number = rtecs_number,
        zvg_number = zvg_number,
        hsdb_id = hsdb_id,
        dsstox_substance_id = dsstox_substance_id,
        dsstox_compound_id = dsstox_compound_id,
        nsc_number = nsc_number,
        ec_number = ec_number,
        echa_substance_infocard_id = echa_substance_infocard_id,
        fl_number = fl_number,
        cosing_number = cosing_number,
        cannabis_database_id = cannabis_database_id,
        nikkaji_id = nikkaji_id,
        hmdb_id = hmdb_id,
        knapsack_id = knapsack_id,
        unii = unii,
        probes_and_drugs_id = probes_and_drugs_id,
        kegg_id = kegg_id,
        lipidmaps_id = lipidmaps_id,
        brockhaus_enzyklopadie_online_id = brockhaus_enzyklopadie_online_id,
        ca_prop_65_id = ca_prop_65_id,
        encyclopaedia_britannica_online_id = encyclopaedia_britannica_online_id,
        freebase_id = freebase_id,
        hcis_id = hcis_id,
        j_global_id = j_global_id,
        jecfa_number = jecfa_number,
        microsoft_academic_id = microsoft_academic_id,
        openalex_id = openalex_id,
        quora_topic_id = quora_topic_id,
        sciencedirect_topic_id = sciencedirect_topic_id,
        surechembl_id = surechembl_id,
        ccdc_number = ccdc_number,
        csd_refcode = csd_refcode,
        icsc_id = icsc_id,
        livertox_id = livertox_id,
        osha_occupational_chemical_database_id = osha_occupational_chemical_database_id,
        pesticideinfo_chemical_id = pesticideinfo_chemical_id,
        atc_code = atc_code,
        who_inn_numeric_id = who_inn_numeric_id,
        guide_to_pharmacology_ligand_id = guide_to_pharmacology_ligand_id,
        mesh_descriptor_id = mesh_descriptor_id,
        mesh_tree_code = mesh_tree_code,
        rxnorm_id = rxnorm_id,
        drugbank_id = drugbank_id,
        medlineplus_id = medlineplus_id,
        medlineplus_drug_identifier = medlineplus_drug_identifier,
        patientslikeme_treatment_id = patientslikeme_treatment_id,
        ndf_rt_id = ndf_rt_id,
        australian_register_of_therapeutic_goods_id = australian_register_of_therapeutic_goods_id,
        drugcentral_id = drugcentral_id,
        iedb_epitope_id = iedb_epitope_id,
        pdb_ligand_id = pdb_ligand_id,
        bbc_things_id = bbc_things_id,
        bncf_thesaurus_id = bncf_thesaurus_id,
        decs_id = decs_id,
        encyclopaedia_universalis_id = encyclopaedia_universalis_id,
        encyclopedia_of_china_third_edition_id = encyclopedia_of_china_third_edition_id,
        gran_enciclopedia_catalana_id = gran_enciclopedia_catalana_id,
        lem_id = lem_id,
        lex_id = lex_id,
        nalt_id = nalt_id,
        ne_se_id = ne_se_id,
        nhs_health_a_to_z_id = nhs_health_a_to_z_id,
        omegawiki_defined_meaning = omegawiki_defined_meaning,
        online_pwn_encyclopedia_id = online_pwn_encyclopedia_id,
        pdb_structure_id = pdb_structure_id
      )

    return(return_result)

  }
