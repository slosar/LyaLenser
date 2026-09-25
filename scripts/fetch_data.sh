#!/bin/bash
# Reproducible download of the public data used by LyaLenser (2026-09-10). Resumable (wget -c).
# Usage: bash fetch_data.sh [desi|act|act_sims|planck|all]   Target: /data/LyaLenser/raw
set -e
RAW=/data/LyaLenser/raw
DESI=https://data.desi.lbl.gov/public/dr1
ACT=https://portal.nersc.gov/project/act/dr6_lensing_v1
what=${1:-all}

desi() {
  mkdir -p $RAW/desi/lya-deltas/delta-lya-0-0/{Delta,Log} $RAW/desi/lss_v1.5 $RAW/desi/qso_iron
  B=$DESI/vac/dr1/lya-deltas/v1.0/delta-lya-0-0
  # Lya forest deltas (picca, 1028 HEALPix files, 5.8 GB) + metadata
  ( cd $RAW/desi/lya-deltas && wget -q -c -O README.md $DESI/vac/dr1/lya-deltas/v1.0/README.md )
  ( cd $RAW/desi/lya-deltas/delta-lya-0-0 && wget -q -c $B/continuum_fitting_mask.txt $B/picca_delta.ini $B/dr1_vac_dr1_lya-deltas_v1.0_delta-lya-0-0.sha256sum )
  ( cd $RAW/desi/lya-deltas/delta-lya-0-0/Log && wget -q -c $B/Log/delta_attributes.fits.gz $B/Log/rejection_log.fits.gz $B/Log/dr1_vac_dr1_lya-deltas_v1.0_delta-lya-0-0_Log.sha256sum )
  curl -sL $B/Delta/ | grep -oE 'delta-[0-9]+\.fits\.gz' | sort -u > $RAW/desi/delta_lya_files.txt
  ( cd $RAW/desi/lya-deltas/delta-lya-0-0/Delta && xargs -P 8 -I{} wget -q -c $B/Delta/{} < $RAW/desi/delta_lya_files.txt )
  # Quasar redshift catalogue used for the Lya analysis (816 MB)
  ( cd $RAW/desi/qso_iron && wget -q -c $DESI/survey/catalogs/dr1/QSO/iron/QSO_cat_iron_cumulative_v0.fits $DESI/survey/catalogs/dr1/QSO/iron/dr1_survey_catalogs_dr1_QSO_iron.sha256sum )
  # LSS clustering catalogues and randoms (v1.5): data + n(z) + 4 randoms per cap (~8 GB)
  L=$DESI/survey/catalogs/dr1/LSS/iron/LSScats/v1.5
  ( cd $RAW/desi/lss_v1.5 && wget -q -c $L/QSO_NGC_clustering.dat.fits $L/QSO_SGC_clustering.dat.fits $L/QSO_NGC_nz.txt $L/QSO_SGC_nz.txt $L/QSO_frac_tlobs.fits $L/dr1_survey_catalogs_dr1_LSS_iron_LSScats_v1.5.sha256sum
    for i in 0 1 2 3; do wget -q -c $L/QSO_NGC_${i}_clustering.ran.fits $L/QSO_SGC_${i}_clustering.ran.fits; done )
}
act() {
  mkdir -p $RAW/act/baseline
  ( cd $RAW/act && wget -q -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/dr6_lensing_release.tar.gz )   # all 13 variants, 1.37 GB
  ( cd $RAW/act/baseline && wget -q -c $ACT/README $ACT/maps/baseline/N_L_kk_act_dr6_lensing_v1_baseline.txt $ACT/maps/baseline/kappa_alm_data_act_dr6_lensing_v1_baseline.fits $ACT/maps/baseline/kappa_filter_act_dr6_lensing_v1_baseline.txt $ACT/maps/baseline/mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits )
}
act_sims() {   # 400 x 153 MB = 61 GB
  mkdir -p $RAW/act/baseline/simulations
  ( cd $RAW/act/baseline/simulations && for i in $(seq -f "%04g" 1 400); do echo $ACT/maps/baseline/simulations/kappa_alm_sim_act_dr6_lensing_v1_baseline_$i.fits; done | xargs -P 2 -n 1 wget -q -c )
}
planck() {   # Carron, Mirmelstein & Lewis 2022 "2018-like" PR4 lensing maps (klm, mask, N_L); simulations are NERSC-only
  mkdir -p $RAW/planck
  ( cd $RAW/planck && wget -q -c https://github.com/carronj/planck_PR4_lensing/releases/download/Data/PR42018like_maps.tar https://github.com/carronj/planck_PR4_lensing/releases/download/Data/NPIPE.zip && tar xf PR42018like_maps.tar && unzip -o -q NPIPE.zip -d NPIPE )
}
case $what in desi) desi;; act) act;; act_sims) act_sims;; planck) planck;; all) desi; act; planck; act_sims;; esac
echo "done: $what"
