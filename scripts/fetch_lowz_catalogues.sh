#!/bin/bash
# DESI DR1 LSS v1.5 clustering catalogues (LRG, ELG_LOPnotqso, BGS_BRIGHT-21.5) with two randoms per cap,
# and the BOSS DR12v5 CMASSLOWZTOT combined sample with one random per cap.
set -u
base=https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5
cd /data/LyaLenser/raw/desi/lss_v1.5
for t in LRG ELG_LOPnotqso BGS_BRIGHT-21.5; do
  for cap in NGC SGC; do
    for f in ${t}_${cap}_clustering.dat.fits ${t}_${cap}_0_clustering.ran.fits ${t}_${cap}_1_clustering.ran.fits ${t}_${cap}_nz.txt; do
      [ -s $f ] && continue
      wget -q -c --retry-on-http-error=503 --tries=20 --waitretry=30 $base/$f -O $f.part && mv $f.part $f && echo "done $f"
    done
  done
done
cd /data/LyaLenser/raw/boss
for f in galaxy_DR12v5_CMASSLOWZTOT_North.fits.gz galaxy_DR12v5_CMASSLOWZTOT_South.fits.gz random0_DR12v5_CMASSLOWZTOT_North.fits.gz random0_DR12v5_CMASSLOWZTOT_South.fits.gz; do
  [ -s $f ] && continue
  wget -q -c --tries=20 --waitretry=30 https://data.sdss.org/sas/dr12/boss/lss/$f -O $f.part && mv $f.part $f && echo "done $f"
done
echo ALL DONE
