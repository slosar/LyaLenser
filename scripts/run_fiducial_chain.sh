#!/bin/bash
# The fiducial DR1 measurement, end to end (docs/pipeline.md stages 2-7), sequential because of memory
# (peaks 34 GB auto, 26 GB cross; about 2 h on a 24-thread workstation). The stage-1 basis tables
# ($D/stageb/basis_dr1_ab_z196.h5, basis_cross_dr1_z196.h5) are built once with build_basis_dr1.py /
# build_basis_cross_dr1.py. Override the variables below in the environment, e.g.
#   LMAX=1300 BANDS="40 200 400 600 800 1000 1300" LOWZ=lowz_v4_l1300 AUTO=dr1_lowz_l1300 CROSS=dr1_qso_l1300 TAG=l1300 ./run_fiducial_chain.sh
set -e
PY=${PY:-/data/LyaLenser/envs/lyalenser/bin/python}
D=${LYALENSER_DATA:-/data/LyaLenser}
R=$(cd "$(dirname "$0")/.." && pwd)
BANDS=${BANDS:-"40 200 400 600 800 1000"}
LMAX=${LMAX:-1000}
LOWZ=${LOWZ:-lowz_v4}; AUTO=${AUTO:-dr1_lowz_v7d}; CROSS=${CROSS:-dr1_qso_v1d}; TAG=${TAG:-v4}
ZEFF=2.3476            # weight-averaged pixel redshift of the 1.96 <= z <= 3.0 sample: the templates' source plane
cd "$R/scripts"
echo "== templates $(date)";        OMP_NUM_THREADS=12 $PY lowz_catalogues.py --out $D/$LOWZ --lmax $LMAX --zref $ZEFF --tracer-zmax 1.6 --wiener measured
echo "== template x CMB $(date)";   OMP_NUM_THREADS=8 $PY deflection_cmb_check.py --lowz $D/$LOWZ --lmax $LMAX --bands $BANDS --out $R/results/deflection_cmb_check_$TAG.json
echo "== forest auto $(date)";      NUMBA_NUM_THREADS=22 OMP_NUM_THREADS=4 $PY run_dr1_lowz.py --out $D/stageb/$AUTO --lowz $D/$LOWZ --basis $D/stageb/basis_dr1_ab_z196.h5 \
    --regions lya lyb --xi-correction spline --randoms 20 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff $ZEFF --bands $BANDS --nside-alpha 2048
cp $D/stageb/$AUTO/dr1_lowz.json $R/results/$AUTO.json; cp $D/stageb/$AUTO/dr1_lowz.md $R/results/$AUTO.md
echo "== quasar x forest $(date)";  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=4 $PY run_dr1_qso_lowz.py --out $D/stageb/$CROSS --auto-run $D/stageb/$AUTO --lowz $D/$LOWZ \
    --basis $D/stageb/basis_cross_dr1_z196.h5 --randoms 20 --sub-slabs 1.96 2.25 2.55 3.0 --spline-fixed-base --bands $BANDS --nside-alpha 2048
cp $D/stageb/$CROSS/dr1_qso.json $R/results/$CROSS.json; cp $D/stageb/$CROSS/dr1_qso.md $R/results/$CROSS.md
echo "== auto redshift split $(date)"; NUMBA_NUM_THREADS=16 $PY auto_subslabs.py --run $D/stageb/$AUTO --lowz $D/$LOWZ --bands $BANDS --nside-alpha 2048
cp $D/stageb/$AUTO/auto_subslabs.json $R/results/auto_subslabs_${AUTO#dr1_lowz_}.json
echo "== joint response fit $(date)"; NUMBA_NUM_THREADS=12 $PY joint_response_fit.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --bands $BANDS --tag $TAG
echo "== combination $(date)";      $PY combine_auto_cross.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --out $R/results/auto_cross_combination_$TAG.json
echo "== injection $(date)";        NUMBA_NUM_THREADS=16 $PY injection_dr1.py --auto $D/stageb/$AUTO --cross $D/stageb/$CROSS --lowz $D/$LOWZ --bands $BANDS --tag $TAG
echo "== CHAIN DONE $(date)"
