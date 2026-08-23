#!/usr/bin/env python3
# Oxidative pentose phosphate genes scored across collections
# Source: ch3-chitin-evolution/scripts/build_ppp_presence.py
# Output: results/pangenome/ppp_presence_verified.tsv
import os, collections
P="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome"
OUT=f"{P}/ppp_presence_verified.tsv"
PID, COVLEN = 30.0, 0    # DIAMOND outfmt here has no qlen; use pident>=30 as the call, evalue already <=1e-5

trio={"UHM979.41089_R.bin.103_CDS_0654":"zwf",
      "UHM1210.23070_R.bin.101_CDS_0293":"gnd",
      "EHM058980_CDS_1432":"g6pd_sub"}

# denominators: genomes actually in each Akkermansia database
def faa_ids(d):
    dd=f"{P}/{d}"
    return set(f[:-4] for f in os.listdir(dd)) if os.path.isdir(dd) else set()
akkfam=faa_ids("akkfam")      # 172
gtdbakk=faa_ids("gtdb_akk")   # 60

# amphibian Akkermansia genomes from genome_groups.tsv (AKK_AMPH = 94)
amph=set()
for l in open(f"{P}/genome_groups.tsv"):
    x=l.rstrip("\n").split("\t")
    if len(x)>1 and x[1]=="AKK_AMPH": amph.add(x[0])

def subj_genome(s):
    # subject id forms: "GB_GCA_x.1|contig_n"  or "Genus::GB_..|.."  or "UHM..._bin.x|contig"
    s=s.split("::")[-1]
    s=s.split("|")[0]
    # strip trailing _contig index if present on MAG names
    return s

def parse(fn):
    """return {query_cds: set(subject_genome)} for pident>=PID"""
    hits=collections.defaultdict(set)
    p=f"{P}/{fn}"
    if not os.path.isfile(p): return hits
    for l in open(p,errors="ignore"):
        x=l.rstrip("\n").split("\t")
        if len(x)<5: continue
        q=x[0]; s=x[1]
        try: pid=float(x[2])
        except: continue
        if pid<PID: continue
        hits[q].add(subj_genome(s))
    return hits

# Which file covers which Akkermansia set:
#  all199    -> NOVEL(105)+AKK_AMPH(94) selfsearch; we take AKK_AMPH subjects only
#  akkother  -> akkfam (172)
#  gtdbakk / gtdbakkdb -> gtdb_akk (60)
files={"all199":"n22_all199.tsv","akkother":"n22_akkother.tsv",
       "gtdbakk":"n22_gtdbakk.tsv","gtdbakkdb":"n22_gtdbakkdb.tsv"}
H={k:parse(v) for k,v in files.items()}

def count_in(query, genomeset, hitmap):
    present={g for g in hitmap.get(query,set()) if g in genomeset}
    return len(present), len(genomeset)

print("PER-GENE PRESENCE (pident>=30, evalue<=1e-5), REAL denominators")
print("="*72)
rows=[]
for cds,gene in trio.items():
    # amphibian Akk: subjects in all199 that are AKK_AMPH genomes
    na,da=count_in(cds, amph, H["all199"])
    # akkfam 172
    nf,df=count_in(cds, akkfam, H["akkother"])
    # gtdb 60 (try both possible files)
    hg = H["gtdbakk"] if H["gtdbakk"] else H["gtdbakkdb"]
    ng,dg=count_in(cds, gtdbakk, hg)
    tot_n=na+nf+ng; tot_d=da+df+dg
    print(f"\n{gene} ({cds})")
    print(f"  AKK_AMPH   {na}/{da}")
    print(f"  akkfam172  {nf}/{df}")
    print(f"  gtdb_akk60 {ng}/{dg}")
    print(f"  TOTAL Akkermansia carrying {gene}: {tot_n}/{tot_d}")
    rows.append((gene,na,da,nf,df,ng,dg,tot_n,tot_d))

with open(OUT,"w") as fh:
    fh.write("gene\tn_akkamph\td_akkamph\tn_akkfam\td_akkfam\tn_gtdb\td_gtdb\tn_total\td_total\n")
    for r in rows: fh.write("\t".join(map(str,r))+"\n")
print(f"\nwrote {OUT}")

print("\n"+"="*72)
print("ALL 22 nohomolog proteins: how many hit ANY Akkermansia genome? (pident>=30)")
print("="*72)
allq=set()
for l in open(f"{P}/nohomolog22.faa"):
    if l.startswith(">"): allq.add(l[1:].strip())
akk_all = akkfam|gtdbakk|amph
merged=collections.defaultdict(set)
for k in ["all199","akkother","gtdbakk","gtdbakkdb"]:
    for q,ss in H[k].items():
        merged[q] |= {g for g in ss if g in akk_all}
clean=[q for q in allq if len(merged.get(q,set()))==0]
print(f"  queries with ZERO Akkermansia hit: {len(clean)}/{len(allq)}")
for q in sorted(clean): print(f"    {q}")
print("\nDONE. Paste all output.")
