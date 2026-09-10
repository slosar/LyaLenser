# Repo-level memory (caveats and portability notes)

- Machine: this work started on a Linux box where `~/work -> ~/Dropbox/work` (symlink). The repo lives in Dropbox, so avoid large data files here; put data under a non-synced path and record it below.
- Python: anaconda 3.11 at `/home/anze/anaconda3/bin/python3` with numpy, scipy, matplotlib, healpy, astropy. No `.venv` yet.
- Mathematica: available both via the `mathematica` MCP server and `/usr/bin/wolframscript`. MCP shows only printed output.
- Codex MCP: failed to connect on 2026-09-10 (CONNECTION_CLOSED). Reviews are queued in PROGRESS.md.
- Data: no DESI DR1 deltas or CMB lensing maps are present locally. Expected source is NERSC; paths to be recorded here when known.
- GitHub: private repo `slosar/LyaLenser` (gh authenticated as slosar).
