# Retina Spectacle

A research workspace for public diabetic-retinopathy and related single-cell studies. The interface is designed around selecting cells, comparing recorded conditions and making expression plots. It is independent of the original Spectacle project.

The website runs on Vercel. Standard GitHub Actions runners retrieve repository matrices and publish processed, versioned data as public release assets. Expression matrices are downloaded when a study is opened, not bundled with the website. No paid runner, storage upgrade or payment method is required; provider free-tier limits and queue availability apply.

## Research workspace

- Pool sample libraries by recorded condition within one study, or group by sample, provisional cell type, cluster, cell type × condition, custom sample groups, or cell selections A/B.
- Include and exclude plot groups. Intersect sample, tissue, cell type, cluster and selection filters across plots and exports. Clearing every option produces an empty plot, not a fallback to all cells.
- Draw box-and-whisker plots with every included cell, including zeros. Switch to a violin distribution or hide individual points. Boxes show Q1–Q3 and the median; whiskers reach observed values within 1.5× IQR.
- Draw dot plots and heatmaps for up to 12 genes. Dot area encodes detection percentage; color encodes the mean log-normalized expression including zeros.
- Select cells with a lasso or exact barcode/index lookup on UMAP/t-SNE. A and B remain disjoint. Shift-lasso adds cells; ordinary lasso replaces that selection.
- Compare up to 12 selected genes in the browser. A worker calculates two-sided asymptotic Mann–Whitney tests with tie/continuity correction and BH adjustment across the genes tested. There are no new-job credentials or paid services involved in this workflow.
- Export plotted cells, summaries, selected-gene comparisons, gene panels and cell selections as CSV; plots as PNG/SVG. Copy a workspace link that includes the study version, gene, filters, grouping and cell selections.
- Browse existing saved differential-expression results and their Excel workbook; reopen the existing saved reclustered subset.

## Data and scientific limits

The source review is [PMC11214886](https://pmc.ncbi.nlm.nih.gov/articles/PMC11214886/). The current 17 ready studies contain 470,416 retained cells. Related models and reference studies are separate from the default diabetes/DR collection. GSE209872 is partial because one deposited control matrix was corrupt.

Recorded conditions are not interchangeable: GSE248284 contains type 1 diabetes with versus without DR, not healthy controls. GSE245561 contains PDR membranes only. OIR is a retinal disease model, not diabetes. GSE204880 defaults to the imported mouse retina samples and keeps kidney separate; its broader GEO series also describes human bulk data. Each workspace exposes its design and condition-label evidence.

Raw counts undergo documented QC, CP10K normalization, log1p, variable-gene selection, scaled PCA, UMAP, t-SNE and Leiden clustering. GSE150703 retains its documented deposited normalized scale. Marker-panel cell types are provisional, not validated author annotations.

Comparisons stay within one study. Cell-level tests do not account for donor replication, batch effects or covariates. Library counts are not necessarily donor counts. GSE178121 has one pooled library per condition. The new comparison corrects only across the selected genes, not the whole transcriptome. Fold change uses log2((mean expm1(expression) in A + 1) / (mean expm1(expression) in B + 1)).

## Operation

Source checks are scheduled every 30 minutes using GitHub Actions. The site displays the last successful check, and offers a newly completed dataset version without changing an open comparison underneath the user. Supported sources are imported automatically; restricted and unsupported sources remain explicitly unavailable. See [OPERATIONS.md](OPERATIONS.md).

Public browsing, in-browser selected-gene comparisons and existing saved results work with no Vercel secrets. The separate backend job API still requires private server configuration for new genome-wide analyses or reclustering; it is not exposed as an available website action. Unfinished upload, regional-expression and sc-eQTL controls are absent.

## Verification

Run `node --test test/research.test.mjs` for 35 cohort, filtering, plotting and statistical checks, including expected results generated independently with SciPy. Run `node test-gateway.cjs` for 21 API/concurrency checks. Public fixtures contain metadata and selected real gene vectors, not entire expression matrices.

Browser validation covers real control/diabetes groups, human DR/non-DR groups, the PDR-only study, mixed-tissue defaults, all-cell plots, empty filters, custom groups, map selections, sharing, exports and responsive layouts. A 150,504-cell study exercises the largest available cohort.

Source inspiration: [Spectacle](https://singlecell-eye.org/app/spectacle/) and [cellcuratoR](https://github.com/drewvoigt10/cellcuratoR). This project is not affiliated with Spectacle or the University of Iowa.
