#!/bin/bash
# Iteration 14: the fiducial (lmax 1000) rebuilt with the measured-covariance Wiener combination (lowz_v4): deflection
# validation, cross and auto runs (fiducial bands, nside 2048), auto redshift split, combination. Products dr1_qso_v1d, dr1_lowz_v7d.
# auto x cross combination. Sequential (memory). Products: lowz_v4, dr1_qso_v1d, dr1_lowz_v7d.
set -e
PY=/data/LyaLenser/envs/lyalenser/bin/python; D=/data/LyaLenser; R=/home/anze/Dropbox/work/LyaLenser; S=$R/code/stageb
B="40 200 400 600 800 1000"
cd $S
echo "== templates lmax 1300 $(date)"; OMP_NUM_THREADS=12 $PY lowz_catalogues.py --out $D/lowz_v4 --zref 2.3476 --tracer-zmax 1.6 --wiener measured
echo "== deflection check $(date)"; OMP_NUM_THREADS=8 $PY deflection_cmb_check.py --lowz $D/lowz_v4 --bands $B --out $R/report/stageb/deflection_cmb_check_v4.json
echo "== cross $(date)"; NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=4 $PY run_dr1_qso_lowz.py --out $D/stageb/dr1_qso_v1d --lowz $D/lowz_v4 --randoms 20 --sub-slabs 1.96 2.25 2.55 3.0 --spline-fixed-base --bands $B --nside-alpha 2048
cp $D/stageb/dr1_qso_v1d/dr1_qso.json $R/report/stageb/dr1_qso_v1d.json; cp $D/stageb/dr1_qso_v1d/dr1_qso.md $R/report/stageb/dr1_qso_v1d.md
echo "== auto $(date)"; NUMBA_NUM_THREADS=22 OMP_NUM_THREADS=4 $PY run_dr1_lowz.py --out $D/stageb/dr1_lowz_v7d --lowz $D/lowz_v4 --basis $D/stageb/basis_dr1_ab_z196.h5 --regions lya lyb --xi-correction spline --randoms 20 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff 2.3476 --bands $B --nside-alpha 2048
cp $D/stageb/dr1_lowz_v7d/dr1_lowz.json $R/report/stageb/dr1_lowz_v7d.json; cp $D/stageb/dr1_lowz_v7d/dr1_lowz.md $R/report/stageb/dr1_lowz_v7d.md
echo "== auto split $(date)"; NUMBA_NUM_THREADS=16 $PY auto_subslabs.py --run $D/stageb/dr1_lowz_v7d --lowz $D/lowz_v4 --bands $B --nside-alpha 2048
cp $D/stageb/dr1_lowz_v7d/auto_subslabs.json $R/report/stageb/auto_subslabs_v7d.json
echo "== combination $(date)"; $PY combine_auto_cross.py --auto $D/stageb/dr1_lowz_v7d --cross $D/stageb/dr1_qso_v1d --out $R/report/stageb/auto_cross_combination_v4.json
echo "== CHAIN DONE $(date)"
