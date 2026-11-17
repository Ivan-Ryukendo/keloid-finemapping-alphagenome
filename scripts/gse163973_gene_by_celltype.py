"""Per-cell-type mean log1p(CP10k) and %cells>0 for GWAS-mapped genes in GSE163973. Run with -I."""
import glob, gzip, re
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
R = 'data/raw/GSE163973/'
f = pd.read_csv('data/processed/finemap_candidates.csv'); g = pd.read_csv('data/raw/keloid_gwas_variants.csv')
mg = g[g.rs_id.isin(f.lead_rsid.unique())].mapped_genes.dropna()
genes = sorted({x.strip() for s in mg for x in s.split(';')})
controls = ['POSTN', 'COL1A1', 'PTPRC', 'PECAM1', 'KRT14']  # sanity markers, not GWAS genes
want = genes + controls
meta = pd.read_csv(R + 'GSE163973_integrate.all.NS.all.KL_cell.meta.data.csv.gz', index_col=0)
rows, found = [], set()
for d in sorted(glob.glob(R + 'raw_tar/*/*_matrix')):
    s = re.search(r'/(\w\w\d)_matrix', d.replace(chr(92), '/')).group(1)
    X = sio.mmread(d + '/matrix.mtx.gz').tocsr()  # genes x cells
    feat = pd.read_csv(d + '/features.tsv.gz', sep='\t', header=None)[1].values
    bc = pd.read_csv(d + '/barcodes.tsv.gz', header=None)[0].str.replace('-1$', '', regex=True).values
    lib = np.asarray(X.sum(0)).ravel(); lib[lib == 0] = 1
    ct = meta.loc[meta['orig.ident'] == s, 'cellType']; ct.index = ct.index.str.replace(r'_\d+$', '', regex=True)
    keep = pd.Series(range(len(bc)), index=bc).reindex(ct.index)
    assert keep.notna().mean() > 0.99, (s, keep.notna().mean())  # barcode join sanity
    idx = keep.dropna().astype(int).values; ctk = ct[keep.notna()].values
    for gn in want:
        w = np.where(feat == gn)[0]
        if len(w) == 0: continue
        found.add(gn)
        v = X[w[0], :].toarray().ravel()[idx]; v = np.log1p(v / lib[idx] * 1e4)
        rows.append(pd.DataFrame({'gene': gn, 'sample': s, 'cell_type': ctk, 'x': v}))
d = pd.concat(rows); d['condition'] = np.where(d['sample'].str.startswith('KF'), 'keloid', 'normal_scar')
out = []
for cond, dd in [('all', d), ('keloid', d[d.condition == 'keloid']), ('normal_scar', d[d.condition == 'normal_scar'])]:
    a = dd.groupby(['gene', 'cell_type']).x.agg(n_cells='size', mean_log1p_cp10k='mean', pct_expressing=lambda x: 100 * (x > 0).mean()).reset_index()
    a.insert(1, 'condition', cond); out.append(a)
o = pd.concat(out); o['role'] = np.where(o.gene.isin(controls), 'marker_control', 'gwas_mapped')
o.to_csv('results/alphagenome/gse163973_gene_by_celltype.csv', index=False)
print('cell types:', sorted(meta.cellType.unique()))
print('found', len(found & set(genes)), 'of', len(genes), '; missing:', sorted(set(genes) - found))
print(meta.groupby('cellType').size().to_string())
