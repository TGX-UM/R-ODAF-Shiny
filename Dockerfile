# R-ODAF Shiny -- reproducible build
#
# Every layer that can move is pinned. The versions below were derived from the
# image running in production (built 2026-04-24, R 4.5.3, 313 packages) and each
# pin was verified against that inventory -- see HANDOVER.md for the evidence.
#
# Three independent sources have to be pinned, or a rebuild silently changes the
# differential-expression results this app reports:
#   1. the base image  -> R itself + the 14 base/recommended packages
#   2. CRAN            -> 271 packages
#   3. Bioconductor    -> 28 packages, incl. DESeq2 and edgeR
#
# Do NOT replace this with `renv::restore()`. The renv.lock committed in this
# repo is stale (R 4.3.3 / Bioconductor 3.17 / DESeq2 1.40.2) and restoring it
# would roll production back two Bioconductor releases.

# 1. Base image, pinned by digest.
#    Tag alone is not enough: rocker/shiny:4.5.3 was rebuilt 2026-06-24, so the
#    tag is a moving target even though the R version in its name is not.
FROM rocker/shiny:4.5.3@sha256:5d2f546daf4aa2cbd0aa5f979d5141c066e9b1d3195fe411f190a3f74c6269d4

# 2. CRAN, pinned to the Posit Package Manager snapshot matching the production
#    build date. The __linux__/noble path serves prebuilt binaries, which is both
#    faster and what production actually installed.
ENV CRAN_SNAPSHOT=https://p3m.dev/cran/__linux__/noble/2026-04-24

# 3. Bioconductor 3.22. That release is archived, so bioconductor.org redirects
#    to the archive -- which is exactly what we want, because an archived release
#    no longer moves. Set as plain repos rather than going through
#    BiocManager::install(version=), so no version negotiation can override us.
ENV BIOC_RELEASE=3.22

# Compile parallelism. Defaults to 0 = "use every core", which is right for a
# GitHub runner. Set it lower when building on a shared cluster node so the
# Bioconductor source compile does not starve the services running alongside it
# (e.g. --build-arg NCPUS=4 on tgx2).
ARG NCPUS=0

RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    build-essential \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libglpk40 \
    libproj-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Point every install at the pinned repos for the rest of the build.
#
# NOTE: this APPENDS -- it must never overwrite the file. rocker's Rprofile.site
# also sets `HTTPUserAgent`, which is what tells Posit Package Manager to serve
# precompiled Linux *binaries* rather than source tarballs
# (https://docs.rstudio.com/rspm/admin/serving-binaries/). Clobber it and every
# package builds from source, at which point installs start failing on system
# libraries that are not in the apt list above -- observed as `usethis` failing
# with "unable to load shared object fs.so: libuv.so.1: cannot open", taking
# devtools, sf, magick, units, s2, rsvg, gt and pkgdown down with it.
# Production installs these as binaries, so it never hits that path.
# A later options() call wins, so appending still overrides the repos.
#
# Only CRAN goes in the global repos. The Bioconductor repos are passed
# explicitly to the one install call that needs them (below) -- deliberately, so
# this reproduces production's *repo visibility at each step*. Production ran its
# CRAN install with only CRAN configured, so `Suggests` that live in
# Bioconductor simply could not resolve and were skipped. Making Bioconductor
# globally visible instead pulls in extras production does not have
# (BiocStyle, bookdown were observed).
#
# The timeout is not optional. Bioconductor 3.22 is archived, so its tarballs are
# served from a slow OSN mirror (mghp.osn.xsede.org), and R's default 60s download
# timeout is not enough for the large annotation packages -- `org.Hs.eg.db` (~180 MB)
# timed out and was silently skipped. See the warning-not-error note below.
RUN printf '%s\n' \
    'options(repos = c(CRAN = Sys.getenv("CRAN_SNAPSHOT")))' \
    'options(timeout = 1800)' \
    "options(Ncpus = if (${NCPUS} > 0) ${NCPUS} else max(1L, parallel::detectCores()))" \
    >> /usr/local/lib/R/etc/Rprofile.site \
    && R -q -e 'stopifnot(nzchar(getOption("HTTPUserAgent"))); cat("binary-serving user agent preserved\n"); print(getOption("repos"))'

# CRAN packages (binaries from the pinned p3m snapshot).
RUN R -q -e 'install.packages(c("DT","data.table","dplyr","tidyr", \
    "pheatmap","ggplot2","cowplot","corrplot","plotly","stringr", \
    "BiocManager","devtools","enrichR"), dependencies=TRUE)' \
    && R -q -e 'stopifnot(all(c("DT","data.table","dplyr","tidyr","pheatmap","ggplot2","cowplot","corrplot","plotly","stringr","BiocManager","devtools","enrichR") %in% rownames(installed.packages())))'

# R ships Matrix, nlme and rpart with the base image, but production upgraded all
# three from CRAN during its dependency resolution. Request them explicitly so the
# pinned snapshot supplies production's versions (1.7-5 / 3.1-169 / 4.1.27) rather
# than leaving the older bundled ones in place. Done before Bioconductor, so
# anything compiled afterwards links against the same Matrix production had.
#
# `lib = .Library` matters: the default target is site-library, which would leave
# the older bundled copies shadowed in the base library rather than replaced.
# R would still load the newer ones (site-library comes first in .libPaths()), so
# behaviour is identical either way -- but `installed.packages()` then reports two
# rows per package, and a naive inventory diff picks up the stale one and reports
# a drift that does not exist. Overwriting in place matches production and keeps
# the inventory unambiguous.
RUN R -q -e 'install.packages(c("Matrix","nlme","rpart"), lib = .Library)'

# Bioconductor packages, from the pinned (archived) release. These compile from
# source, which is what production did -- its Bioc packages all show Built=4.5.3.
# Repos are passed per-call rather than set globally; see the Rprofile note above.
RUN R -q -e 'bioc <- c( \
      BioCsoft = paste0("https://bioconductor.org/packages/", Sys.getenv("BIOC_RELEASE"), "/bioc"), \
      BioCann  = paste0("https://bioconductor.org/packages/", Sys.getenv("BIOC_RELEASE"), "/data/annotation")); \
    install.packages(c("DESeq2","edgeR","AnnotationDbi","org.Hs.eg.db","BiocVersion"), repos = c(getOption("repos"), bioc)); \
    install.packages("PCAtools", repos = c(getOption("repos"), bioc))'

# Fail the build loudly if the packages that determine the DE results drifted,
# rather than shipping a subtly different image. Update these expectations only
# together with a deliberate, validated version bump.
#
# This guard is load-bearing, not belt-and-braces. `install.packages()` reports a
# failed download as a *warning*, so the R process — and therefore the RUN step —
# exits 0 even when a package was never installed. That is exactly what happened
# on a validation build: `org.Hs.eg.db` timed out against the archived
# Bioconductor mirror, the build reported success, and only this check caught that
# the image was missing the package the app needs for human gene symbols.
# Never remove it, and never let it degrade into a warning.
RUN R -q -e 'exp <- c(DESeq2="1.50.2", edgeR="4.8.2", PCAtools="2.22.4", \
      org.Hs.eg.db="3.22.0", AnnotationDbi="1.72.0"); \
    got <- vapply(names(exp), function(p) as.character(packageVersion(p)), ""); \
    bad <- got != exp; \
    if (any(bad)) stop("pinned package drift: ", \
      paste0(names(exp)[bad], " expected ", exp[bad], " got ", got[bad], collapse="; ")); \
    cat("pinned versions OK\n"); print(data.frame(expected=exp, got=got))'

RUN rm -rf /srv/shiny-server/*
WORKDIR /srv/shiny-server/
COPY ./ ./

EXPOSE 3838
