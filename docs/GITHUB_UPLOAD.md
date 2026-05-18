# GitHub upload instructions

## First-time upload

```bash
cd /path/to/rnaseq2tracksP

git init
git add .
git commit -m "Initial commit: rnaseq2tracks v5.0"
git branch -M main

# Add remote (create repo first on github.com)
git remote add origin https://github.com/MichalGd/rnaseq2tracksP.git
git push -u origin main

# Tag release
git tag -a v5.0 -m "rnaseq2tracks v5.0 — enrichment analysis"
git push origin v5.0
```

## Updating an existing repo

```bash
cd /path/to/rnaseq2tracksP

git add .
git commit -m "v5.0 — gene set enrichment analysis"
git push origin main

git tag -a v5.0 -m "v5.0"
git push origin v5.0
```

## Verify .gitignore

Ensure these are NOT committed:
```bash
git status | grep -E "config\.conf|samplesheet\.csv|contrasts\.csv|\.RData|_output/"
# Should return empty
```

## Update MANIFEST

After any file changes, regenerate the SHA256 manifest:

```bash
find scripts/ docs/ config/ examples/ tests/ README.md CHANGELOG.md \
  LICENSE CITATION.cff FILE_TREE.txt environment.yml -type f \
  | sort | xargs sha256sum > MANIFEST.sha256
git add MANIFEST.sha256 && git commit -m "Update MANIFEST for v5.0"
```

## Update ORCID

Edit `CITATION.cff` and replace `0000-0000-0000-0000` with your real ORCID before making the repository public.

## Config templates

`config/config.conf`, `config/samplesheet.csv`, and `config/contrasts.csv` are excluded
by `.gitignore` because they contain server-specific paths. The corresponding template
files in `config/` (e.g., `config_template.conf`) and the examples in `examples/` are
committed and serve as documentation for users setting up the pipeline.
