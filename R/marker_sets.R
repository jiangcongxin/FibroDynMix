#' Get built-in FibroDynMix marker priors
#'
#' Returns lightweight fibroblast state marker priors for common species and
#' contexts. These priors are intended as weak starting points and should be
#' reviewed against tissue-specific marker expression before formal use.
#'
#' @param species Species name: `human` or `mouse`.
#' @param context Biological context. Supported values include `generic`,
#'   `skin`, `scar`, and `caf`.
#'
#' @return A named list of marker genes per FibroDynMix state.
#' @export
get_fibrodynmix_markers <- function(species = c("human", "mouse"),
                                    context = c("generic", "skin", "scar", "caf")) {
  species <- match.arg(species)
  context <- match.arg(context)

  markers <- switch(
    species,
    human = human_fibrodynmix_markers(context),
    mouse = mouse_fibrodynmix_markers(context)
  )
  class(markers) <- c("FibroDynMixMarkerSet", class(markers))
  attr(markers, "species") <- species
  attr(markers, "context") <- context
  markers
}

human_fibrodynmix_markers <- function(context) {
  base <- list(
    resident = c("DCN", "LUM", "COL14A1", "PDGFRA", "PI16"),
    inflammatory = c("IL6", "CXCL12", "CXCL14", "CCL2", "CXCL2"),
    myofibroblast = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1"),
    `ECM-remodeling` = c("COL1A1", "COL1A2", "FN1", "POSTN", "MMP2"),
    `antigen-presenting` = c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1"),
    `IFN-stress` = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1")
  )
  if (context %in% c("skin", "scar")) {
    base$resident <- unique(c(base$resident, "APOD", "FBLN1"))
    base$myofibroblast <- unique(c(base$myofibroblast, "PDGFRB", "MCAM"))
    base$`ECM-remodeling` <- unique(c(base$`ECM-remodeling`, "COL3A1", "COMP", "THBS2"))
  }
  if (context == "scar") {
    base$inflammatory <- unique(c(base$inflammatory, "CXCL8", "PTGS2"))
    base$`ECM-remodeling` <- unique(c(base$`ECM-remodeling`, "TNC", "LOX"))
  }
  if (context == "caf") {
    base$myofibroblast <- unique(c(base$myofibroblast, "FAP", "PDPN"))
    base$`ECM-remodeling` <- unique(c(base$`ECM-remodeling`, "COL11A1", "MMP11"))
    base$inflammatory <- unique(c(base$inflammatory, "CXCL1", "CXCL8"))
  }
  base
}

mouse_fibrodynmix_markers <- function(context) {
  base <- list(
    resident = c("Dcn", "Lum", "Col14a1", "Pdgfra", "Pi16"),
    inflammatory = c("Il6", "Cxcl12", "Cxcl14", "Ccl2", "Cxcl2"),
    myofibroblast = c("Acta2", "Tagln", "Myl9", "Tpm2", "Cnn1"),
    `ECM-remodeling` = c("Col1a1", "Col1a2", "Fn1", "Postn", "Mmp2"),
    `antigen-presenting` = c("H2-Aa", "H2-Ab1", "Cd74", "H2-Eb1", "H2-DMa"),
    `IFN-stress` = c("Isg15", "Ifit1", "Ifit3", "Mx1", "Oas1")
  )
  if (context %in% c("skin", "scar")) {
    base$resident <- unique(c(base$resident, "Apod", "Fbln1"))
    base$myofibroblast <- unique(c(base$myofibroblast, "Pdgfrb", "Mcam"))
    base$`ECM-remodeling` <- unique(c(base$`ECM-remodeling`, "Col3a1", "Comp", "Thbs2"))
  }
  if (context == "scar") {
    base$inflammatory <- unique(c(base$inflammatory, "Cxcl1", "Ptgs2"))
    base$`ECM-remodeling` <- unique(c(base$`ECM-remodeling`, "Tnc", "Lox"))
  }
  if (context == "caf") {
    base$myofibroblast <- unique(c(base$myofibroblast, "Fap", "Pdpn"))
    base$`ECM-remodeling` <- unique(c(base$`ECM-remodeling`, "Col11a1", "Mmp11"))
    base$inflammatory <- unique(c(base$inflammatory, "Cxcl1", "Cxcl5"))
  }
  base
}
