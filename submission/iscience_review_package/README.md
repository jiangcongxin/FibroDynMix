# iScience Review Package

This package gathers the current iScience-facing manuscript materials for
internal review before journal upload. It is a review package, not yet a final
Editorial Manager upload bundle.

## Package Contents

| Folder | Contents | Review purpose |
|---|---|---|
| `manuscript/` | Current iScience manuscript draft | Main text, Summary, Highlights, STAR Methods, references |
| `frontmatter/` | eTOC, graphical abstract plan, cover-letter fit text | Submission metadata and front matter review |
| `cover_letter/` | Full cover letter draft, fit paragraph, and one-sentence fit statement | Editorial-fit language |
| `figures/main/` | Figure 1-6 full packages | Main figures, exports, legends, source data, manifests, QC |
| `figures/supplementary/` | Supplementary Figure S1-S4 packages | Supplementary figures, exports, legends, source data, manifests |
| `source_data_and_manifests/` | Claim-source alignment and reference gap documents | Claim boundary and source-data audit |
| `quality_control/` | Review-package QC outputs | Package inventory and validation summaries |

## Current Manuscript Spine

- Main figures: Figure 1-6.
- Supplementary figures: S1-S4.
- References: 38, all currently cited in the manuscript.
- Summary length: 126 words in the packaged manuscript.
- eTOC blurb: 35-word primary version, 32-word backup version.

## Remaining Human Inputs Before Submission

- Final author list, affiliations, and corresponding author.
- ORCID IDs, especially for corresponding author.
- Funding statement values and grant numbers.
- Author contributions using the chosen contribution taxonomy.
- Declaration of interests form and manuscript statement.
- Final repository or archive links for code and data, if the local workspace is
  not the final public access point.
- Final journal-upload figure naming and file-format preferences.

## Quality Gate

The source workspace passed the current checks after package assembly:

- Citation/reference consistency: all 38 references cited; no unused references.
- Project integrity: passed for 7 required figure packages and analysis outputs.
- Local unit tests: 359 passing tests.
- `.DS_Store`: none found.
- `R CMD build .`: passed after excluding `submission/` from the R package tarball.
- `R CMD check --no-manual --no-build-vignettes`: status OK.

## Added Submission Templates

- `frontmatter/title_page_template.md`
- `frontmatter/keywords.md`
- `frontmatter/declarations_template.md`
- `cover_letter/cover_letter_full_draft.md`
