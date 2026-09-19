#!/bin/bash
# Iteration 13b: the lmax = 1300 convergence test. Templates rebuilt to lmax 1300, deflection validation, the
# cross and auto runs with bands 40-200 ... 1000-1300 (deflection at nside 2048), the auto redshift split and the
# auto x cross combination. Sequential (memory). Products: lowz_v3_l1300, dr1_qso_v1c, dr1_lowz_v7c.
set -e
PY=/data/LyaLenser/envs/lyalenser/bin/python; D=/data/LyaLenser; R=/home/anze/Dropbox/work/LyaLenser; S=$R/code/stageb
B="40 200 400 600 800 1000 1300"
cd $S
echo "== templates lmax 1300 $(date)"; OMP_NUM_THREADS=12 $PY lowz_catalogues.py --out $D/lowz_v3_l1300 --lmax 1300 --zref 2.3476 --tracer-zmax 1.6
echo "== deflection check $(date)"; OMP_NUM_THREADS=8 $PY deflection_cmb_check.py --lowz $D/lowz_v3_l1300 --lmax 1300 --bands $B --out $R/report/stageb/deflection_cmb_check_l1300.json
echo "== cross $(date)"; NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=4 $PY run_dr1_qso_lowz.py --out $D/stageb/dr1_qso_v1c --lowz $D/lowz_v3_l1300 --randoms 20 --sub-slabs 1.96 2.25 2.55 3.0 --spline-fixed-base --bands $B --nside-alpha 2048
cp $D/stageb/dr1_qso_v1c/dr1_qso.json $R/report/stageb/dr1_qso_v1c.json; cp $D/stageb/dr1_qso_v1c/dr1_qso.md $R/report/stageb/dr1_qso_v1c.md
echo "== auto $(date)"; NUMBA_NUM_THREADS=22 OMP_NUM_THREADS=4 $PY run_dr1_lowz.py --out $D/stageb/dr1_lowz_v7c --lowz $D/lowz_v3_l1300 --basis $D/stageb/basis_dr1_ab_z196.h5 --regions lya lyb --xi-correction spline --randoms 20 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff 2.3476 --bands $B --nside-alpha 2048
cp $D/stageb/dr1_lowz_v7c/dr1_lowz.json $R/report/stageb/dr1_lowz_v7c.json; cp $D/stageb/dr1_lowz_v7c/dr1_lowz.md $R/report/stageb/dr1_lowz_v7c.md
echo "== auto split $(date)"; NUMBA_NUM_THREADS=16 $PY auto_subslabs.py --run $D/stageb/dr1_lowz_v7c --lowz $D/lowz_v3_l1300 --bands $B --nside-alpha 2048
cp $D/stageb/dr1_lowz_v7c/auto_subslabs.json $R/report/stageb/auto_subslabs_v7c.json
echo "== combination $(date)"; $PY combine_auto_cross.py --auto $D/stageb/dr1_lowz_v7c --cross $D/stageb/dr1_qso_v1c --out $R/report/stageb/auto_cross_combination_l1300.json
echo "== CHAIN DONE $(date)"
