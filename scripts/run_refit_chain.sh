#!/bin/bash
# Iteration 15 (derivative-map source-distance treatment): refit the fiducial measurement from the saved
# catalogues of the previous runs with the templates that carry derivative maps, and redo every product
# that depends on the fit (docs/pipeline.md stages 5-7). Sequential (memory); about 4 h on the workstation.
#   ./run_refit_chain.sh            (defaults below; the robustness variants run on NERSC, slurm/refit.sbatch)
set -e
PY=${PY:-/data/LyaLenser/envs/lyalenser/bin/python}
D=${LYALENSER_DATA:-/data/LyaLenser}
R=$(cd "$(dirname "$0")/.." && pwd)
BANDS=${BANDS:-"40 200 400 600 800 1000"}
LOWZ=${LOWZ:-lowz_v5}; TAG=${TAG:-v5}
OLD_AUTO=${OLD_AUTO:-dr1_lowz_v7d}; OLD_CROSS=${OLD_CROSS:-dr1_qso_v1d}
AUTO=${AUTO:-dr1_lowz_v8}; CROSS=${CROSS:-dr1_qso_v2}
ZEFF=2.3476
cd "$R/scripts"
echo "== forest auto refit $(date)";  NUMBA_NUM_THREADS=22 OMP_NUM_THREADS=4 $PY run_dr1_lowz.py --out $D/stageb/$AUTO --refit-from $D/stageb/$OLD_AUTO --lowz $D/$LOWZ --basis $D/stageb/basis_dr1_ab_z196.h5 \
    --regions lya lyb --xi-correction spline --randoms 20 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff $ZEFF --bands $BANDS --nside-alpha 2048
cp $D/stageb/$AUTO/dr1_lowz.json $R/results/$AUTO.json; cp $D/stageb/$AUTO/dr1_lowz.md $R/results/$AUTO.md
echo "== quasar x forest refit $(date)"; NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=4 $PY run_dr1_qso_lowz.py --out $D/stageb/$CROSS --refit-from $D/stageb/$OLD_CROSS --auto-run $D/stageb/$AUTO --lowz $D/$LOWZ \
    --basis $D/stageb/basis_cross_dr1_z196.h5 --randoms 20 --sub-slabs 1.96 2.25 2.55 3.0 --spline-fixed-base --bands $BANDS --nside-alpha 2048
cp $D/stageb/$CROSS/dr1_qso.json $R/results/$CROSS.json; cp $D/stageb/$CROSS/dr1_qso.md $R/results/$CROSS.md
echo "== joint response fit $(date)"; NUMBA_NUM_THREADS=12 $PY joint_response_fit.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --bands $BANDS --tag $TAG
echo "== joint fit, no derivative $(date)"; NUMBA_NUM_THREADS=12 $PY joint_response_fit.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --bands $BANDS --tag ${TAG}_noderiv --no-derivative
echo "== combination $(date)";      $PY combine_auto_cross.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --out $R/results/auto_cross_combination_$TAG.json
echo "== auto redshift split $(date)"; NUMBA_NUM_THREADS=16 $PY auto_subslabs.py --run $D/stageb/$AUTO --lowz $D/$LOWZ --bands $BANDS --nside-alpha 2048
cp $D/stageb/$AUTO/auto_subslabs.json $R/results/auto_subslabs_${AUTO#dr1_lowz_}.json
echo "== size of the source-distance term $(date)"; $PY source_distance_expansion.py --run $D/stageb/$AUTO --lowz $D/$LOWZ --summary $R/results/$AUTO.json --previous $R/results/$OLD_AUTO.json
echo "== scale sensitivity $(date)"; NUMBA_NUM_THREADS=16 $PY scale_sensitivity.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --bands $BANDS --tag $TAG
echo "== random-template nulls, 100 draws $(date)"; mkdir -p $D/audit_$TAG
NUMBA_NUM_THREADS=16 $PY random_template_null.py --statistic auto --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --out $D/audit_$TAG/null_auto_0.json
NUMBA_NUM_THREADS=16 $PY random_template_null.py --statistic cross --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --out $D/audit_$TAG/null_cross_0.json
$PY summarise_audit.py --directory $D/audit_$TAG --out $R/results/random_template_null_100_$TAG.json --auto-run $AUTO --cross-run $CROSS
echo "== injection $(date)";        NUMBA_NUM_THREADS=16 $PY injection_dr1.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --bands $BANDS --tag $TAG
echo "== CHAIN DONE $(date)"
