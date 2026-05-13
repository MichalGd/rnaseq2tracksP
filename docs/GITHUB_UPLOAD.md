# GitHub upload

```bash
cd /path/to/rnaseq2tracks
git init && git add . && git commit -m "Initial commit: rnaseq2tracks v1.0"
git branch -M main
git remote add origin https://github.com/MichalGd/rnaseq2tracks.git
git push -u origin main
git tag -a v1.0 -m "v1.0" && git push origin v1.0
```

Verify excluded files (should NOT appear in `git status`):
- `config/config.conf`
- `config/samplesheet.csv`
- Any output directories
