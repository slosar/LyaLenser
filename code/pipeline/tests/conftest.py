import sys
from pathlib import Path

PIPE=Path(__file__).resolve().parents[1]
CODE=PIPE.parent
for p in (str(PIPE),str(CODE)):
    if p not in sys.path: sys.path.insert(0,p)

