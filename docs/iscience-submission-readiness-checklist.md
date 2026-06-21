# iScience Submission Readiness Checklist

Date: 2026-06-15

## Package Location

Review package:
`submission/iscience_review_package/`

## Official Guidance Checked

- iScience author information and journal policies.
- iScience article-type guidance, including concise article length and STAR
  Methods placement.
- Cell Press graphical abstract guidance: single-panel visual, sparse labels,
  distinct from main-text figures, and readable as a square artwork.
- Cell Press/Elsevier Key Resources Table guidance: essential datasets,
  software, source data, and identifiers should be listed and also described in
  STAR Methods.
- Cell Press declaration-of-interests guidance: declaration form and manuscript
  statement must be consistent.

## Manuscript Structure

| Item | Status | Notes |
|---|---|---|
| Title | Ready for review | Current title: "FibroDynMix models fibroblast state plasticity as raw-count latent mixtures" |
| Summary | Ready for review | 126 words in the review package; no citations |
| Highlights | Ready for review | 4 bullets; should be checked against final character-limit screen |
| Graphical Abstract section | Ready for review | Points to separate 1200 x 1200 px artwork draft |
| Introduction | Ready for review | Method gap and bounded evidence sequence are stated |
| Results | Ready for review | Claim-style headings and figure sequence are aligned |
| Discussion | Ready for review | Four result-first paragraphs |
| Limitations of the study | Ready for review | Diagnostic, causal, therapeutic, atlas-level, and lineage boundaries retained |
| STAR Methods | Ready for review | KRT, resource availability, public dataset details, public GitHub/Zenodo links, method details, statistics, and quality gates are present |
| References | Ready for review | 38 references, all cited |

## Front Matter and Submission Metadata

| Item | Status | Notes |
|---|---|---|
| eTOC blurb | Ready for review | Primary version 35 words; backup version 32 words |
| Graphical abstract artwork | Draft ready | SVG source plus PNG/PDF/TIFF/JPG exports in review package |
| Cover letter | Draft ready | Full draft plus fit paragraph included; requires author/corresponding-author details |
| Keywords | Draft ready | 6 recommended keywords plus 5-keyword backup set |
| Title page | Template ready | Requires author list, affiliations, corresponding author, and ORCID IDs |
| Funding statement | Template ready | Requires author confirmation and grant details |
| Author contributions | Template ready | Requires final author list and contribution roles |
| Declaration of interests | Template ready | Requires author form and exact manuscript statement |
| Generative AI / writing-assistance statement | Template ready | Use only if required by final journal workflow or institutional policy |

## Figures and Source Data

| Item | Status | Notes |
|---|---|---|
| Main Figure 1-6 packages | Ready for review | Exports, legends, manifests, source data, palettes, and QC files copied |
| Supplementary Figure S1-S4 packages | Ready for review | Supplementary packages copied with legends, manifests, and source data |
| Graphical abstract exports | Draft ready | 1200 x 1200 px PNG/JPG/TIFF plus PDF and SVG |
| Figure claim-result-source alignment | Ready for review | No open issues in current alignment report |
| Source-data manifests | Ready for review | Panel-level source-data manifests included for all figure packages |
| Final upload naming | Pending | Rename after journal file-format decision |

## Data, Code, and Reproducibility

| Item | Status | Notes |
|---|---|---|
| Data availability statement | Ready for review | Public accessions/DOI listed; generated summaries, figure source data, manifests, and QC files archived with GitHub/Zenodo release |
| Code availability statement | Ready for review | R package source and analysis scripts public at GitHub and archived at Zenodo DOI 10.5281/zenodo.20787528 |
| Key resources table | Ready for review | Public datasets, package source, R version, conventions, and source-data outputs listed |
| Runtime information | Ready for review | R 4.6.0 and runtime records stated |
| Project integrity check | Passed | 7 figure packages and required analysis outputs checked |
| Unit tests | Passed | 359 tests passed |
| R CMD check | Passed | `R CMD check --no-manual --no-build-vignettes` status OK on 2026-06-15 |

## Claim Boundary

| Claim area | Boundary |
|---|---|
| Fibroblast states | State mixtures are inferred from raw counts with weak marker orientation; marker genes do not directly define final weights |
| Simulation | Simulation truth supports latent-state recovery claims |
| Public data | Public datasets support execution, transfer, marker-gradient recovery, and directional scar-context evidence |
| Scar module | Directional effect-size evidence only; not FDR-significant |
| Transition flow/FPI | Cross-sectional state-composition summary; not observed lineage tracing |
| Clinical relevance | No diagnostic, therapeutic, causal disease-mechanism, or completed human atlas claim |

## Final Freeze Gate

Before journal upload, rerun:

```bash
Rscript scripts/check_project_integrity.R
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual --no-build-vignettes FibroDynMix_0.0.0.9000.tar.gz
```

Latest final-freeze run on 2026-06-15:

- `Rscript scripts/check_project_integrity.R`: passed.
- `Rscript -e 'testthat::test_local()'`: 359 passed, 0 failed, 0 warnings, 0 skipped.
- `R CMD build .`: passed after excluding `submission/` from the R package tarball.
- `R CMD check --no-manual --no-build-vignettes FibroDynMix_0.0.0.9000.tar.gz`: status OK.
