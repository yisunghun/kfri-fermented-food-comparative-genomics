import pandas as pd
import re, os, glob

def derive_sample_id(raw):
    base = os.path.basename(str(raw))
    return re.sub(r'\.(tab|fna|fa|fasta|txt|csv|tsv)$', '', base)

vfdb = pd.read_csv('/mnt/f/WGS_Consolidated/abricate_out/vfdb_summary.tsv', sep='\t')
vfdb['sample_id'] = vfdb[vfdb.columns[0]].apply(derive_sample_id)
vfdb = vfdb.set_index('sample_id')

fg_map = {}
fg_root = '/mnt/f/WGS_Consolidated/grouped/by_functional_group'
for fg_name in os.listdir(fg_root):
    fg_dir = os.path.join(fg_root, fg_name)
    gff_dir = os.path.join(fg_dir, 'genomes_gff')
    search_dir = gff_dir if os.path.isdir(gff_dir) else fg_dir
    for fpath in glob.glob(os.path.join(search_dir, '*.gff')):
        sid = os.path.splitext(os.path.basename(fpath))[0]
        fg_map[sid] = fg_name

vfdb['functional_group'] = vfdb.index.map(lambda s: fg_map.get(s))
vfdb = vfdb[vfdb['functional_group'].isin(['LAB', 'Bacillus_group'])]

for fg in ['LAB', 'Bacillus_group']:
    sub = vfdb[vfdb['functional_group'] == fg]['NUM_FOUND']
    pct_positive = (sub > 0).mean() * 100
    print(f'{fg}: n={len(sub)}, {pct_positive:.1f}% 시료가 1개 이상 보유, 평균={sub.mean():.3f}, 최댓값={sub.max()}')
    print(f'  분포: {sub.value_counts().sort_index().to_dict()}')
