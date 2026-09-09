import openpyxl
import collections

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
M = CH3 + "data/251027_16S_metadata_final.xlsx"
S = CH3 + "results/rurik_16s/matched_samples.txt"
OUT = CH3 + "results/rurik_16s/host_summary.txt"

matched = set(l.strip() for l in open(S) if l.strip())
print("matched samples:", len(matched))

wb = openpyxl.load_workbook(M, read_only=True, data_only=True)
ws = wb.worksheets[0]
rows = ws.iter_rows(values_only=True)
hdr = ['' if c is None else str(c) for c in next(rows)]
print("metadata columns:", len(hdr))
print(hdr[:20])
print("")

def find(*keys):
    for i, h in enumerate(hdr):
        hl = h.lower()
        if any(k in hl for k in keys):
            return i
    return None

ci = find("sample", "sampleid", "id")
hi = find("host_species", "species")
oi = find("order", "host_order")
mi = find("management", "captiv", "wild")
print("using columns: sample=%s host=%s order=%s status=%s" % (
    hdr[ci] if ci is not None else None,
    hdr[hi] if hi is not None else None,
    hdr[oi] if oi is not None else None,
    hdr[mi] if mi is not None else None))
print("")

hit_h, hit_o, hit_m = collections.Counter(), collections.Counter(), collections.Counter()
all_o = collections.Counter()
found = 0
for row in rows:
    c = ['' if x is None else str(x).strip() for x in row]
    if ci is None or ci >= len(c):
        continue
    sid = c[ci]
    if oi is not None and oi < len(c):
        all_o[c[oi]] += 1
    if sid in matched:
        found += 1
        if hi is not None and hi < len(c): hit_h[c[hi]] += 1
        if oi is not None and oi < len(c): hit_o[c[oi]] += 1
        if mi is not None and mi < len(c): hit_m[c[mi]] += 1
wb.close()

L = []
def say(s):
    print(s)
    L.append(s)

say("matched samples found in metadata: %d of %d" % (found, len(matched)))
say("")
say("=== HOST ORDER OF MATCHED SAMPLES ===")
for k, v in hit_o.most_common():
    tot = all_o.get(k, 0)
    say("   %-22s %4d of %4d in the dataset" % (k or "(blank)", v, tot))
say("")
say("=== MANAGEMENT STATUS ===")
for k, v in hit_m.most_common():
    say("   %-22s %4d" % (k or "(blank)", v))
say("")
say("=== HOST SPECIES, top 30 ===")
for k, v in hit_h.most_common(30):
    say("   %-40s %4d" % (k or "(blank)", v))
say("")
say("distinct host species with a match: %d" % len(hit_h))

open(OUT, "w").write("\n".join(L) + "\n")
print("")
print("WROTE:", OUT)
