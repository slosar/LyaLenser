"""Write report/main_short.bib from main.bib with author lists truncated to N names + 'and others'."""
import re, sys
N = 6
src = open("../report/main.bib").read()
def trunc(m):
    authors = [a.strip() for a in re.split(r"\s+and\s+", m.group(2))]
    if len(authors) > N:
        authors = authors[:N] + ["others"]
    return m.group(1) + " and ".join(authors) + m.group(3)
out = re.sub(r"(author\s*=\s*\{)(.*?)(\},?\n)", trunc, src, flags=re.S)
open("../report/main_short.bib", "w").write(out)
print("wrote main_short.bib")
