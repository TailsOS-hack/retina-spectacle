# Retina Spectacle

An independent adaptation of the Spectacle interface for public retinal single-cell research. The website runs on Vercel. Standard GitHub Actions runners retrieve repository matrices, process datasets and run requested analyses. Versioned release assets hold the processed public research data. No paid runner, paid storage upgrade or payment method is required; provider free-tier limits and queue availability apply.

## Data and methods

The source review is [PMC11214886](https://pmc.ncbi.nlm.nih.gov/articles/PMC11214886/). Each dataset retains its source accession and downloaded-file provenance. The 17 migrated studies contain 470,416 retained cells. GSE209872 is partial because one deposited control matrix was corrupt. Reference studies and mixed tissues are labeled; these are not all direct diabetic-retinal cohorts.

Raw counts undergo documented gene/cell filtering, CP10K normalization, log1p, variable-gene selection, scaled PCA, UMAP, t-SNE and Leiden clustering. GSE150703 retains its verified deposited normalized scale. Differential expression is an exploratory cell-level Wilcoxon test with BH correction; it does not replace independent biological replication.

## Interface

Data Select follows the original Spectacle study table. Dashboard pairs dimensionality reduction / heatmap / group 1 with violin plots / reclustering / group 2. Regional anatomical expression, sc-eQTL and private uploads require further datasets/workflows and are clearly marked unavailable. This project is not affiliated with Spectacle or the University of Iowa. Original inspiration: https://singlecell-eye.org/app/spectacle/ and https://github.com/drewvoigt10/cellcuratoR . Font Awesome is distributed under its included license.

## Operation

Visitors request shared repository checks at most once per 30-minute window. A GitHub workflow imports supported new or changed sources. Existing immutable versions remain readable during processing. Owner authentication is required for analysis jobs and manual retries. Public API endpoints never disclose server credentials. Analysis requests refer only to public cell indices.

Vercel environment: `GH_TOKEN` (repository workflow/content access), `ADMIN_KEY` (analysis access key), `GITHUB_REPOSITORY` (owner/name). Browser configuration contains only `/api`. Keep credentials out of commits. GitHub Actions uses its built-in repository token. All workflows use standard free `ubuntu-24.04` runners; no larger runners or paid storage addons are configured.
