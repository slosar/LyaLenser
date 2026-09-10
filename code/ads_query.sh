#!/bin/bash
Q="$1"; ROWS="${2:-50}"
curl -s -G "https://api.adsabs.harvard.edu/v1/search/query" \
  -H "Authorization: Bearer $ADS_API_TOKEN" \
  --data-urlencode "q=$Q" \
  --data-urlencode "fl=bibcode,title,author,year,pub,volume,page,identifier,citation_count" \
  --data-urlencode "rows=$ROWS" --data-urlencode "sort=date asc" \
| python3 /tmp/claude-1000/-home-anze-Dropbox-work-LyaLenser/6ea34145-6049-4f0e-a9e5-b61457afca9a/scratchpad/parse.py
