# Implementation round 3 (gpt-6-astra) final summary — 2026-09-11 (verbatim task output)

Implemented the iteration-3 changes and rebuilt validation from the fresh ensemble. **36 tests pass; 20/30 acceptance gates pass, 10 fail. Stage B remains blocked.**

Numerical choices and [GATES.md](/home/anze/Dropbox/work/LyaLenser/GATES.md) were frozen before validation. The prescribed precision rule triggered extension through seed 59; no subsequent tuning occurred.

| Acceptance gate | Result |
|---|---|
| Truth normalization slope | PASS |
| CMB normalization slope | FAIL |
| Matched normalization slope | FAIL |
| Fixed versus refitted baseline | PASS |
| Varying-template mean field: truth | PASS |
| Fixed-template mean field: truth | PASS |
| Varying-template mean field: CMB | PASS |
| Fixed-template mean field: CMB | PASS |
| Varying-template mean field: matched | FAIL |
| Fixed-template mean field: matched | PASS |
| Absolute stochastic recovery, two-SEM consistency | PASS |
| Response-only CMB prediction agreement | PASS |
| Response-only matched prediction agreement | FAIL |
| Absolute deprojected response-only null | FAIL |
| Absolute combined deprojected recovery | FAIL |
| Six-bin Hotelling shape | PASS |
| Covariance calibration | FAIL |
| Outside-box and total spectra | PASS |
| Analytic derivative convergence | PASS |
| Discrete-grid covariance quadrature | PASS |
| Unit and regression tests | PASS |
| Xi interpolation convergence | PASS |
| Dense physical normalization, g on | FAIL |
| Data versus analytic baseline | FAIL |
| First-moment correction | PASS |
| Dense physical normalization, g off | FAIL |
| Injection bookkeeping, including joint curl-error audit | PASS |
| Flag and common-realization margin diagnostics | PASS |
| 100 random-template diagnostics | PASS |
| Benchmark and memory | PASS |

The main unresolved absolute results are:

- Combined recovery: **−1.6339 ± 0.7132 SEM**, versus target 1; achieved 95% residual bound **4.0764 A**, exceeding 0.3 A.
- Deprojected response-only null: **−2.9660 ± 0.9563 SEM**; achieved bound **4.9004 A**.
- Truth slope: **0.99918 ± 0.02737 SEM** passes the frozen point tolerance, but its **0.05618** residual bound does not establish 3% precision at 95%.
- Covariance ratio: **0.695932**, below the required 0.7.
- Dense normalization slopes remain approximately **1.07**; measured/analytic slope ratio is **1.19362**.

An additional audit found reused RNG-stream prefixes between density generation and CMB noise. Its impact is unresolved; it was documented without changing the frozen campaign. Bias origins remain unisolated.

Validation took **5 h 30 min 23 s**: ensemble **4 h 5 min 11 s**, controls **1 h 24 min**. Tests took **50.9 s**. Validation peak RSS was **27.76 GB**, using **24 Numba threads**. An earlier interrupted development run has an unknown terminal peak, documented in the notes.

Full measurements, SEMs, bounds, and diagnoses are in [the validation report](/home/anze/Dropbox/work/LyaLenser/report/mock_validation.md), [JSON](/home/anze/Dropbox/work/LyaLenser/report/mock_validation.json), and [iteration-3 notes](/home/anze/Dropbox/work/LyaLenser/code/pipeline/NOTES.md:193). Disk reconstruction verified **2,229 fits** and reproduced all acceptance rows exactly. All four figures were visually checked.

:codex-file-citation{path="/home/anze/Dropbox/work/LyaLenser/report/figures/mock_physical_normalisation.pdf" purpose="output"}
:codex-file-citation{path="/home/anze/Dropbox/work/LyaLenser/report/figures/mock_combined_ensemble.pdf" purpose="output"}
:codex-file-citation{path="/home/anze/Dropbox/work/LyaLenser/report/figures/mock_shape.pdf" purpose="output"}
:codex-file-citation{path="/home/anze/Dropbox/work/LyaLenser/report/figures/mock_spectra.pdf" purpose="output"}

Protected files are unchanged, `git diff --check` passes, and no commit was made.
