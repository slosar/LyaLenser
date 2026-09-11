# Iteration 4 campaign on BNL RACF

Do not run scale-1 commands on the login node. The orchestrating session submits
these commands as separate HTCondor jobs. All numerical/source/GATES changes
must precede the five development jobs. No validation-seed tuning is permitted.

Every job runs from the repository's `code/pipeline` directory with:

```bash
export LYALENSER_DATA=/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser
export NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
export MPLCONFIGDIR=/tmp/lyalenser-mpl-${USER}
export LYALENSER_REPO=/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser
export LYALENSER_PYTHON=/gpfs/mnt/gpfs02/astro/workarea/anze/envs/lyalenser/bin/python
cd "$LYALENSER_REPO/code/pipeline"
```

Use the following ordered commands. Each line inside a loop denotes an
independent job command, **not** a serial loop to run on the login node.
All commands use `--scale 1 --mock-root "$LYALENSER_DATA/mocks/iteration4"`.
Request shared filesystem access (`should_transfer_files = NO`, `getenv = True`).
Request 20 GB disk for sparse/dev jobs and 40 GB for dense/control jobs. Aggregate
campaign storage should have at least 150 GB free; preserve completion markers
and all artifacts. Log paths must be outside each phase's product directory.

1. Five independent development jobs: **8 CPUs, 32 GB RAM each**.

```bash
for seed in 100 101 102 103 104; do
  "$LYALENSER_PYTHON" run_mock_validation.py --phase dev-seed --seed "$seed" --scale 1 --mock-root "$LYALENSER_DATA/mocks/iteration4"
done
```

2. After all five succeed, freeze: **8 CPUs, 32 GB RAM**.

```bash
"$LYALENSER_PYTHON" run_mock_validation.py --phase freeze --scale 1 --mock-root "$LYALENSER_DATA/mocks/iteration4"
```

This reads each completed development product, checks provenance and file hashes,
selects the declared smoothing width and copies the saved covariance table. It
writes `freeze/development.json`, `freeze/frozen.json`, and `freeze/complete.json`.
It does not generate a mock. Never copy the smoke freeze into this campaign.

3. Independent sparse seeds: **8 CPUs, 32 GB RAM each**.

```bash
for seed in $(seq 0 59); do
  "$LYALENSER_PYTHON" run_mock_validation.py --phase seed --seed "$seed" --variant sparse --scale 1 --mock-root "$LYALENSER_DATA/mocks/iteration4"
done
```

Seeds 0–19 are recovery, 20–39 null, 40–59 prospective extension. Computing the
extension in advance is permitted; collect includes it only if the frozen
stopping criterion on the initial forty triggers. Alternatively submit 0–39,
run collect after dense/control jobs, and submit all 40–59 only if collect writes
`extension_decision.json` with `extend: true` and exits for missing extension
products. No mock is generated in collect.

4. Independent dense seeds: **16 CPUs, 48 GB RAM each** (28 pair catalogues of ~1e11 pixel pairs; ~5 h at 8 threads).
   `campaign4.sh` in this directory submits steps 2-6 in this order with these requests and waits between them.

```bash
for seed in $(seq 200 209); do
  "$LYALENSER_PYTHON" run_mock_validation.py --phase seed --seed "$seed" --variant dense --scale 1 --mock-root "$LYALENSER_DATA/mocks/iteration4"
done
```

Each job performs both g-on and g-off with 100 sightlines/deg², zero forest
pixel noise and A_true={0,0.5,1,2}. Dense ranges never overlap development or
sparse seeds. The map-level audit is persisted for dense g-on as well.

5. Independent controls after sparse seed 0: **8 CPUs, 48 GB RAM each**.

```bash
for name in numerical injection flags random benchmark; do
  "$LYALENSER_PYTHON" run_mock_validation.py --phase control --name "$name" --scale 1 --mock-root "$LYALENSER_DATA/mocks/iteration4"
done
```

`flags` generates the three flag-off mocks within its separately schedulable job.
`numerical` computes analytic derivative convergence and runs pytest. `random`
collects the 100 masked Gaussian-template fits saved by sparse seed 0. `injection`
and `benchmark` read seed 0. All controls are idempotent and provenance checked.

6. After required seeds and controls finish, collect: **8 CPUs, 48 GB RAM**.

```bash
"$LYALENSER_PYTHON" run_mock_validation.py --phase collect --scale 1 --mock-root "$LYALENSER_DATA/mocks/iteration4" --output "$LYALENSER_REPO/report"
```

This verifies hashes, applies the stopping rule once, aggregates absolute
statistics, blocks deprojection interpretation on failed map prerequisites, and
renders JSON/Markdown/figures. `--phase rebuild` repeats the same read-only
assembly. Neither phase generates a mock. Missing or changed products are errors.
A completion JSON with matching provenance and artifact hashes makes rerunning
any dev/seed/control phase a no-op.

The implementation session uses the same sequence at `--scale 0.25` and
`--mock-root "$LYALENSER_DATA/mocks/iteration4_smoke"`, four threads, sparse
seeds 0 and 1 and dense seeds 200 and 201. The result is smoke, never acceptance.
