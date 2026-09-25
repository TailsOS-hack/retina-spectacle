# Operating Retina Spectacle

The website is https://retina-spectacle.vercel.app. Dataset versions and completed analyses are public GitHub release assets. The frontend contains no bundled expression matrices.

## Source checks

`process.yml` is active on the default `main` branch and schedules checks at minutes 17 and 47 of each hour. GitHub may delay or skip scheduled starts under load. The website displays the last successful source-check timestamp, which is the relevant freshness measure.

An owner can also open the repository's Actions tab, choose **Process public retinal data**, select **Run workflow**, and use `kind=refresh`. No payment method is required. Check the run and the resulting `state/catalog.json` rather than assuming a scheduled time proves a successful import.

Every check re-reads the motivating review, screens current GEO search results, checks existing sources, and imports up to three supported new or changed studies. Repeated failed sources wait a day before another automatic attempt. Completed versions remain readable. The job's optional `target` accession forces a retry. Repository indexing and headers do not guarantee discovery of every relevant publication or unannounced byte change.

## New website analyses

The public-read deployment intentionally has no account credential installed. Public browsing and completed results work without authentication. New website analysis submissions and visitor-triggered scans remain disabled until the owner authorizes the Vercel server-secret setup.

Use a token limited to this repository when possible. The API needs GitHub Contents and Actions write permissions. Set `GH_TOKEN`, `ADMIN_KEY` and `GITHUB_REPOSITORY` as server environment variables, then redeploy. The browser receives only `/api`, never these server values. Owners enter their analysis access code in the sidebar; it is retained for that browser session only.

## Results and recovery

Analysis requests use public cell indices and a dataset revision. Identical requests reuse the same job. Jobs store their status under `state/jobs` and immutable results in `job-<id>` releases. If a worker stops after uploading results, rerunning that job recovers the existing immutable files.

The gateway rejects unpublished dataset revisions, invalid cell indices, overlapping comparisons, and subsets outside its size limits. Controlled-access repositories, unsupported formats and raw-read-only records stay explicitly unavailable. Do not remove source provenance or label reference studies as diabetic retinal cohorts.

## Limits

This independently implements the core expression workflows and Spectacle-style interface. Regional anatomical expression, sc-eQTL, private uploads, cross-study integration and donor-aware tests are not included. Cell labels are provisional. Free runner availability and provider quotas apply; a passing bounded test is not an unlimited-capacity guarantee.

Provider references: [GitHub scheduled events](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#schedule), [standard hosted runners](https://docs.github.com/en/actions/reference/runners/github-hosted-runners), [release assets](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases), and [Vercel function limits](https://vercel.com/docs/functions/limitations).
