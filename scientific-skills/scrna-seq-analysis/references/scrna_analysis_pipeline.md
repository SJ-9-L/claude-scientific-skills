# scRNA-seq Analysis Pipeline Reference

## Complete Analysis Workflow

### Step 1: Data Loading

```python
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

# Load 10X Genomics data
adata = sc.read_10x_mtx(
    'path/to/filtered_feature_bc_matrix/',
    var_names='gene_symbols',
    cache=True
)

# Load H5AD file
adata = sc.read_h5ad('data.h5ad')

# Load from H5 file (10X)
adata = sc.read_10x_h5('data.h5')

# Load from CSV/TSV
adata = sc.read_csv('expression_matrix.csv')
adata = adata.T  # Transpose if genes are rows

# Load multiple samples
samples = ['sample1', 'sample2', 'sample3']
adatas = []
for sample in samples:
    a = sc.read_10x_h5(f'{sample}/outs/filtered_feature_bc_matrix.h5')
    a.obs['sample'] = sample
    adatas.append(a)
adata = ad.concat(adatas, label='sample')
```

### Step 2: Quality Control

```python
# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Human
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # Mouse
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt', 'ribo'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# QC Visualizations
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_qc.pdf')

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
              color='pct_counts_mt', save='_qc_scatter.pdf')

# Apply filters
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=5000)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 10, :]

# Doublet detection (optional)
import scrublet as scr
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets
adata = adata[~adata.obs['predicted_doublet'], :]
```

### Step 3: Normalization

```python
# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Library size normalization
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transformation
sc.pp.log1p(adata)

# Store normalized data for visualization
adata.raw = adata

# Alternative: scran normalization (for heterogeneous populations)
# import scran
# from scipy import sparse
# adata_pp = adata.copy()
# sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
# sc.pp.log1p(adata_pp)
# sc.pp.pca(adata_pp, n_comps=15)
# sc.pp.neighbors(adata_pp)
# sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
# input_groups = adata_pp.obs['groups']
# data_mat = adata.X.T
# size_factors = scran.computeSumFactors(data_mat, clusters=input_groups)
# adata.obs['size_factors'] = size_factors
# adata.X /= adata.obs['size_factors'].values[:, None]
# sc.pp.log1p(adata)
```

### Step 4: Feature Selection

```python
# Highly Variable Genes (HVGs)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor='seurat_v3',  # or 'seurat', 'cell_ranger'
    batch_key='sample'   # for batch-aware HVG selection
)

# Visualize HVGs
sc.pl.highly_variable_genes(adata, save='_hvg.pdf')

# Subset to HVGs
adata = adata[:, adata.var.highly_variable]

# Regress out confounding factors
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Scale data
sc.pp.scale(adata, max_value=10)
```

### Step 5: Dimensionality Reduction

```python
# PCA
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Determine optimal number of PCs
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='_variance.pdf')

# Neighborhood graph
sc.pp.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=30,  # Adjust based on variance ratio plot
    metric='cosine'
)

# UMAP
sc.tl.umap(
    adata,
    min_dist=0.3,
    spread=1.0,
    n_components=2
)

# t-SNE (optional)
sc.tl.tsne(adata, perplexity=30, n_pcs=30)

# Diffusion map (for trajectory analysis)
sc.tl.diffmap(adata, n_comps=15)
```

### Step 6: Clustering

```python
# Leiden clustering (recommended)
resolutions = [0.3, 0.5, 0.8, 1.0, 1.2]
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')

# Visualize different resolutions
sc.pl.umap(adata, color=['leiden_0.3', 'leiden_0.5', 'leiden_0.8', 'leiden_1.0'],
           legend_loc='on data', save='_clustering_resolutions.pdf')

# Select optimal resolution
adata.obs['leiden'] = adata.obs['leiden_0.8']  # Adjust based on biological knowledge

# Hierarchical clustering of clusters
sc.tl.dendrogram(adata, groupby='leiden')
```

### Step 7: Differential Expression

```python
# Find marker genes for each cluster
sc.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    method='wilcoxon',  # or 't-test', 'logreg'
    pts=True  # Calculate percentage of cells expressing
)

# Visualize markers
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_markers.pdf')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save='_markers_dotplot.pdf')
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, save='_markers_heatmap.pdf')

# Get marker table
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv('marker_genes.csv', index=False)

# Filter significant markers
sig_markers = markers[
    (markers['pvals_adj'] < 0.05) &
    (markers['logfoldchanges'] > 1)
]

# Pairwise comparison
sc.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    groups=['0'],  # Test cluster
    reference='1',  # Reference cluster
    method='wilcoxon'
)
```

### Step 8: Cell Type Annotation

```python
# Manual annotation with known markers
marker_genes = {
    'T cells': ['CD3D', 'CD3E', 'CD4', 'CD8A'],
    'B cells': ['MS4A1', 'CD79A', 'CD19'],
    'NK cells': ['NKG7', 'GNLY', 'NCAM1'],
    'Monocytes': ['CD14', 'LYZ', 'S100A8'],
    'Dendritic': ['FCER1A', 'CD1C'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
    'Fibroblasts': ['COL1A1', 'DCN', 'LUM']
}

# Score cells for each signature
for cell_type, genes in marker_genes.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if genes_present:
        sc.tl.score_genes(adata, genes_present, score_name=f'{cell_type}_score')

# Visualize scores
sc.pl.umap(adata, color=[f'{ct}_score' for ct in marker_genes.keys()],
           save='_celltype_scores.pdf')

# Dot plot of markers
sc.pl.dotplot(adata, var_names=marker_genes, groupby='leiden',
              save='_marker_dotplot.pdf')

# Assign cell types
cluster_to_celltype = {
    '0': 'T cells',
    '1': 'B cells',
    '2': 'Monocytes',
    # ... add more
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)

# Automated annotation with CellTypist
# import celltypist
# predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl')
# adata.obs['celltypist'] = predictions.predicted_labels
```

### Step 9: Trajectory Analysis

```python
# PAGA for trajectory inference
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color='leiden', threshold=0.03, save='_paga.pdf')

# Initialize UMAP with PAGA
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color='leiden', save='_paga_graph.pdf')

# Diffusion pseudotime
root_cluster = '0'  # Set root cluster
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == root_cluster)[0]
sc.tl.dpt(adata)

# Visualize pseudotime
sc.pl.umap(adata, color=['dpt_pseudotime', 'leiden'],
           save='_pseudotime.pdf')

# Gene expression along pseudotime
sc.pl.paga_path(
    adata,
    nodes=['0', '1', '2', '3'],  # Path through clusters
    keys=['gene1', 'gene2', 'gene3'],
    save='_path.pdf'
)
```

### Step 10: Integration & Batch Correction

```python
# Harmony integration
import scanpy.external as sce
sce.pp.harmony_integrate(adata, 'batch', basis='X_pca', adjusted_basis='X_pca_harmony')
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# scVI integration
# import scvi
# adata_scvi = adata.copy()
# scvi.model.SCVI.setup_anndata(adata_scvi, batch_key='batch')
# model = scvi.model.SCVI(adata_scvi)
# model.train()
# adata.obsm['X_scVI'] = model.get_latent_representation()
# sc.pp.neighbors(adata, use_rep='X_scVI')
# sc.tl.umap(adata)

# Combat batch correction
sc.pp.combat(adata, key='batch')
```

## Novel Gene Discovery Pipeline

```python
def discover_novel_genes(adata, reference_markers, cluster_key='leiden',
                         n_top=50, novelty_threshold=0.8):
    """
    Identify potentially novel marker genes not in existing literature.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    reference_markers : dict
        Known markers per cell type from literature
    cluster_key : str
        Cluster column in obs
    n_top : int
        Top N markers to consider per cluster
    novelty_threshold : float
        Minimum fraction of expression in target vs other clusters
    """

    # Get marker genes
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')

    novel_candidates = {}

    for cluster in adata.obs[cluster_key].unique():
        # Get top markers
        markers_df = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        top_markers = markers_df.head(n_top)['names'].tolist()

        # Get known markers for this cluster
        known = set()
        for cell_type, genes in reference_markers.items():
            known.update(genes)

        # Find novel candidates
        novel = []
        for gene in top_markers:
            if gene not in known:
                # Check expression specificity
                expr_in_cluster = adata[adata.obs[cluster_key] == cluster, gene].X.mean()
                expr_in_others = adata[adata.obs[cluster_key] != cluster, gene].X.mean()

                if expr_in_cluster > 0 and expr_in_others > 0:
                    specificity = expr_in_cluster / (expr_in_cluster + expr_in_others)
                    if specificity >= novelty_threshold:
                        novel.append({
                            'gene': gene,
                            'specificity': specificity,
                            'expr_in_cluster': expr_in_cluster,
                            'expr_in_others': expr_in_others
                        })

        novel_candidates[cluster] = novel

    return novel_candidates
```

## Parameter Guidelines

### QC Thresholds by Tissue Type

| Tissue | min_genes | max_genes | max_mt% |
|--------|-----------|-----------|---------|
| PBMC | 200 | 4000 | 10 |
| Tumor | 200 | 6000 | 20 |
| Brain | 500 | 8000 | 5 |
| Liver | 200 | 5000 | 15 |
| Kidney | 200 | 5000 | 15 |
| Embryo | 500 | 6000 | 10 |

### Clustering Resolution Guidelines

| Cell Number | Low Resolution | Medium | High |
|-------------|----------------|--------|------|
| < 5,000 | 0.2-0.4 | 0.5-0.8 | 1.0+ |
| 5,000-20,000 | 0.3-0.5 | 0.6-1.0 | 1.2+ |
| 20,000-100,000 | 0.4-0.6 | 0.8-1.2 | 1.5+ |
| > 100,000 | 0.5-0.8 | 1.0-1.5 | 2.0+ |

### HVG Selection Guidelines

| Analysis Type | n_top_genes |
|--------------|-------------|
| Quick exploration | 1000-2000 |
| Standard analysis | 2000-3000 |
| Rare population detection | 3000-5000 |
| Trajectory analysis | 2000-4000 |
