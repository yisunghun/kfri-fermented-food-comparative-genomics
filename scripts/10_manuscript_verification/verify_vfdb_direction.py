"""
VFDB burden의 실제 방향(어느 쪽이 순위기반으로 유의하게 더 높은지)을
1) dereplicated 데이터(3.9절), 2) genus-stratified 데이터(3.10절) 각각에서
one-sided Mann-Whitney로 명확히 확인한다.
"""
import pandas as pd
import re, os
from scipy import stats

def derive_sample_id(raw):
    base = os.path.basename(str(raw))
    return re.sub(r'\.(tab|fna|fa|fasta|txt|csv|tsv)$', '', base)

# ===== 1) Dereplicated (3.9절) =====
vfdb_derep = pd.read_csv(
    "/mnt/f/WGS_Consolidated/dereplication/filtered/vfdb_summary_dereplicated.tsv", sep="\t"
)
vfdb_derep["sample_id"] = vfdb_derep[vfdb_derep.columns[0]].apply(derive_sample_id)
vfdb_derep = vfdb_derep.set_index("sample_id")

master_derep = pd.read_csv(
    "/mnt/f/WGS_Consolidated/dereplication/filtered/master_table_qc_dereplicated.tsv", sep="\t"
)
master_derep.columns = [re.sub(r"^\ufeff", "", c) for c in master_derep.columns]
master_derep = master_derep.set_index("sample_id")

import glob
fg_map = {}
fg_root = "/mnt/f/WGS_Consolidated/grouped/by_functional_group"
for fg_name in os.listdir(fg_root):
    fg_dir = os.path.join(fg_root, fg_name)
    gff_dir = os.path.join(fg_dir, "genomes_gff")
    search_dir = gff_dir if os.path.isdir(gff_dir) else fg_dir
    for fpath in glob.glob(os.path.join(search_dir, "*.gff")):
        sid = os.path.splitext(os.path.basename(fpath))[0]
        fg_map[sid] = fg_name

vfdb_derep["functional_group"] = vfdb_derep.index.map(lambda s: fg_map.get(s))
vfdb_derep = vfdb_derep[vfdb_derep["functional_group"].isin(["LAB", "Bacillus_group"])]

lab = vfdb_derep.loc[vfdb_derep["functional_group"] == "LAB", "NUM_FOUND"]
bac = vfdb_derep.loc[vfdb_derep["functional_group"] == "Bacillus_group", "NUM_FOUND"]

print("=== 1) Dereplicated (3.9절 대상 데이터) ===")
print(f"LAB n={len(lab)}, 평균={lab.mean():.3f}, 양성비율={100*(lab>0).mean():.1f}%")
print(f"Bacillus_group n={len(bac)}, 평균={bac.mean():.3f}, 양성비율={100*(bac>0).mean():.1f}%")
u, p2 = stats.mannwhitneyu(lab, bac, alternative="two-sided")
_, p_greater = stats.mannwhitneyu(lab, bac, alternative="greater")
_, p_less = stats.mannwhitneyu(lab, bac, alternative="less")
print(f"양측 p={p2:.4g} | LAB>Bacillus_group 단측 p={p_greater:.4g} | LAB<Bacillus_group 단측 p={p_less:.4g}")
print()

# ===== 2) Genus-stratified (3.10절) =====
print("=== 2) Genus-stratified (3.10절 대상 데이터) ===")
vfdb_full = pd.read_csv("/mnt/f/WGS_Consolidated/abricate_out/vfdb_summary.tsv", sep="\t")
vfdb_full["sample_id"] = vfdb_full[vfdb_full.columns[0]].apply(derive_sample_id)
vfdb_full = vfdb_full.set_index("sample_id")
vfdb_full["functional_group"] = vfdb_full.index.map(lambda s: fg_map.get(s))

master_full = pd.read_csv("/mnt/f/WGS_Consolidated/master_table_qc.tsv", sep="\t")
master_full.columns = [re.sub(r"^\ufeff", "", c) for c in master_full.columns]
master_full = master_full.set_index("sample_id")
vfdb_full["genus_final"] = vfdb_full.index.map(lambda s: master_full.loc[s, "genus_final"] if s in master_full.index else None)

df = vfdb_full[vfdb_full["functional_group"].isin(["LAB", "Bacillus_group"])].dropna(subset=["genus_final"])
genus_n = df.groupby("genus_final").size()
valid_genera = genus_n[genus_n >= 3].index
df = df[df["genus_final"].isin(valid_genera)]
genus_means = df.groupby(["genus_final", "functional_group"])["NUM_FOUND"].mean().reset_index()

lab_g = genus_means.loc[genus_means["functional_group"] == "LAB", "NUM_FOUND"]
bac_g = genus_means.loc[genus_means["functional_group"] == "Bacillus_group", "NUM_FOUND"]
print(f"LAB {len(lab_g)}개 속: {sorted(lab_g.round(3).tolist())}")
print(f"Bacillus_group {len(bac_g)}개 속: {sorted(bac_g.round(3).tolist())}")
u2, p2b = stats.mannwhitneyu(lab_g, bac_g, alternative="two-sided")
_, pg_greater = stats.mannwhitneyu(lab_g, bac_g, alternative="greater")
_, pg_less = stats.mannwhitneyu(lab_g, bac_g, alternative="less")
print(f"양측 p={p2b:.4g} | LAB>Bacillus_group 단측 p={pg_greater:.4g} | LAB<Bacillus_group 단측 p={pg_less:.4g}")
