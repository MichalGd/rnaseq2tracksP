# GitHub upload instructions

## First-time upload

```bash
cd /path/to/rnaseq2tracks

# Initialize
git init
git add .
git commit -m "Initial commit: rnaseq2tracks v4.0"
git branch -M main

# Add remote (create repo first on github.com)
git remote add origin https://github.com/MichalGd/rnaseq2tracks.git
git push -u origin main

# Tag release
git tag -a v4.0 -m "rnaseq2tracks v4.0 — RNA-seq QC module"
git push origin v4.0
```

## Verify .gitignore

Ensure these are NOT committed:
```bash
git status | grep -E "config.conf|samplesheet.csv|contrasts.csv|\.RData|output/"
# Should return empty
```

## Update MANIFEST

```bash
find scripts/ docs/ config/ examples/ tests/ README.md CHANGELOG.md \
  LICENSE CITATION.cff environment.yml -type f | sort | xargs sha256sum > MANIFEST.sha256
git add MANIFEST.sha256 && git commit -m "Update MANIFEST"
```

## Update ORCID

Edit `CITATION.cff` and replace `0000-0000-0000-0000` with your ORCID before publishing.
