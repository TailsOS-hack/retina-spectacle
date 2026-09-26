# Operating Retina Spectacle

The website is https://retina-spectacle.vercel.app. Dataset versions and completed analyses are public GitHub release assets. The frontend contains no bundled expression matrices.

## Source checks

`process.yml` is active on the default `main` branch and schedules checks at minutes 17 and 47 of each hour. GitHub may delay or skip scheduled starts under load. The website displays the last successful source-check timestamp, which is the relevant freshness measure.

An owner can also open the repository's Actions tab, choose **Process public retinal data**, select **Run workflow**, and use `kind=refresh`. No payment method is required. Check the run and the resulting `state/catalog.json` rather than assuming a scheduled time proves a successful import.

Every check re-reads the motivating review, screens current GEO search results, checks existing sources, and imports up to three supported new or changed studies. Repeated failed sources wait a day before another automatic attempt. Completed versions remain readable. The job's optional `target` accession forces a retry. Repository indexing and headers do not guarantee discovery of every relevant publication or unannounced byte change.

## Website comparisons and backend jobs

The public deployment has no account credential installed. Selected-gene comparisons run in a browser worker on the filtered public expression vectors, without a server job or owner code. Public browsing, plots, CSV/image exports and completed saved analyses also work without authentication.

New genome-wide analyses or reclustering through the backend API still require separately authorized server credentials. Those submission controls are absent from the public UI. Optional server configuration uses a repository-limited credential as `GH_TOKEN`, a private `ADMIN_KEY`, and `GITHUB_REPOSITORY`; none belongs in browser files or commits. Visitor-triggered source scans also require this optional setup; the existing scheduled source checks do not.

## Results and recovery

Analysis requests use public cell indices and a dataset revision. Identical requests reuse the same job. Jobs store their status under `state/jobs` and immutable results in `job-<id>` releases. If a worker stops after uploading results, rerunning that job recovers the existing immutable files.

The gateway rejects unpublished dataset revisions, invalid cell indices, overlapping comparisons, and subsets outside its size limits. Controlled-access repositories, unsupported formats and raw-read-only records stay explicitly unavailable. Do not remove source provenance or label reference studies as diabetic retinal cohorts.

## Limits

This independently implements expression workflows in a comparison-focused interface. Regional anatomical expression, sc-eQTL, private uploads, cross-study integration and donor-aware tests are not included. Cell labels are provisional. Browser comparisons correct p values across at most 12 selected genes and are not donor-aware or whole-transcriptome tests. Free runner availability and provider quotas apply; a passing bounded test is not an unlimited-capacity guarantee.

Provider references: [GitHub scheduled events](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#schedule), [standard hosted runners](https://docs.github.com/en/actions/reference/runners/github-hosted-runners), [release assets](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases), and [Vercel function limits](https://vercel.com/docs/functions/limitations).
