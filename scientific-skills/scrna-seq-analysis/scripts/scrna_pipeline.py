#!/usr/bin/env python3
"""
scRNA-seq Analysis Pipeline

Comprehensive pipeline for single-cell RNA sequencing data analysis
including QC, normalization, clustering, DE analysis, and novel gene discovery.
"""

import warnings
warnings.filterwarnings('ignore')

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd


class ScRNAAnalyzer:
    """
    Comprehensive scRNA-seq analysis pipeline.

    Features:
    - Data loading (10X, h5ad, CSV)
    - Quality control
    - Normalization and preprocessing
    - Dimensionality reduction
    - Clustering
    - Differential expression
    - Novel gene discovery
    """

    def __init__(self, verbosity: int = 1, figure_dir: str = './figures'):
        """
        Initialize analyzer.

        Parameters
        ----------
        verbosity : int
            Scanpy verbosity level (0-3)
        figure_dir : str
            Directory for saving figures
        """
        import scanpy as sc

        self.sc = sc
        sc.settings.verbosity = verbosity
        sc.settings.figdir = figure_dir
        sc.settings.set_figure_params(dpi=100, facecolor='white')

        self.adata = None
        self.qc_metrics = {}
        self.clustering_results = {}
        self.marker_genes = {}

    def load_data(
        self,
        path: str,
        format: str = 'auto'
    ) -> 'anndata.AnnData':
        """
        Load scRNA-seq data from various formats.

        Parameters
        ----------
        path : str
            Path to data file or directory
        format : str
            Data format: 'auto', '10x_mtx', '10x_h5', 'h5ad', 'csv', 'loom'

        Returns
        -------
        AnnData
            Loaded data object
        """
        path = Path(path)

        if format == 'auto':
            if path.is_dir():
                format = '10x_mtx'
            elif path.suffix == '.h5':
                format = '10x_h5'
            elif path.suffix == '.h5ad':
                format = 'h5ad'
            elif path.suffix in ['.csv', '.tsv']:
                format = 'csv'
            elif path.suffix == '.loom':
                format = 'loom'

        if format == '10x_mtx':
            self.adata = self.sc.read_10x_mtx(path, var_names='gene_symbols')
        elif format == '10x_h5':
            self.adata = self.sc.read_10x_h5(str(path))
        elif format == 'h5ad':
            self.adata = self.sc.read_h5ad(str(path))
        elif format == 'csv':
            self.adata = self.sc.read_csv(str(path))
        elif format == 'loom':
            self.adata = self.sc.read_loom(str(path))

        print(f"Loaded data: {self.adata.n_obs} cells x {self.adata.n_vars} genes")

        return self.adata

    def quality_control(
        self,
        adata: Optional['anndata.AnnData'] = None,
        min_genes: int = 200,
        max_genes: int = 5000,
        min_cells: int = 3,
        max_mt_pct: float = 10.0,
        max_ribo_pct: Optional[float] = None,
        detect_doublets: bool = False,
        plot: bool = True
    ) -> 'anndata.AnnData':
        """
        Perform quality control filtering.

        Parameters
        ----------
        adata : AnnData, optional
            Data object (uses self.adata if not provided)
        min_genes : int
            Minimum genes per cell
        max_genes : int
            Maximum genes per cell
        min_cells : int
            Minimum cells per gene
        max_mt_pct : float
            Maximum mitochondrial percentage
        max_ribo_pct : float, optional
            Maximum ribosomal percentage
        detect_doublets : bool
            Run doublet detection
        plot : bool
            Generate QC plots

        Returns
        -------
        AnnData
            Filtered data
        """
        if adata is None:
            adata = self.adata

        n_cells_before = adata.n_obs
        n_genes_before = adata.n_vars

        # Identify mitochondrial and ribosomal genes
        adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
        adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL', 'Rps', 'Rpl'))

        # Calculate QC metrics
        self.sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=['mt', 'ribo'],
            percent_top=None,
            log1p=False,
            inplace=True
        )

        # Store pre-filter metrics
        self.qc_metrics['pre_filter'] = {
            'n_cells': n_cells_before,
            'n_genes': n_genes_before,
            'median_genes': np.median(adata.obs['n_genes_by_counts']),
            'median_counts': np.median(adata.obs['total_counts']),
            'median_mt_pct': np.median(adata.obs['pct_counts_mt'])
        }

        # Generate QC plots
        if plot:
            self.sc.pl.violin(
                adata,
                ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                jitter=0.4,
                multi_panel=True,
                save='_qc_violin.pdf'
            )

            self.sc.pl.scatter(
                adata,
                x='total_counts',
                y='n_genes_by_counts',
                color='pct_counts_mt',
                save='_qc_scatter.pdf'
            )

        # Apply filters
        self.sc.pp.filter_cells(adata, min_genes=min_genes)
        self.sc.pp.filter_genes(adata, min_cells=min_cells)

        adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
        adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :]

        if max_ribo_pct is not None:
            adata = adata[adata.obs.pct_counts_ribo < max_ribo_pct, :]

        # Doublet detection
        if detect_doublets:
            try:
                import scrublet as scr
                scrub = scr.Scrublet(adata.X)
                doublet_scores, predicted_doublets = scrub.scrub_doublets()
                adata.obs['doublet_score'] = doublet_scores
                adata.obs['predicted_doublet'] = predicted_doublets
                adata = adata[~adata.obs['predicted_doublet'], :]
            except ImportError:
                print("Scrublet not installed. Skipping doublet detection.")

        # Store post-filter metrics
        self.qc_metrics['post_filter'] = {
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'cells_removed': n_cells_before - adata.n_obs,
            'genes_removed': n_genes_before - adata.n_vars
        }

        print(f"QC complete: {adata.n_obs} cells, {adata.n_vars} genes remaining")
        print(f"  Removed {n_cells_before - adata.n_obs} cells, {n_genes_before - adata.n_vars} genes")

        self.adata = adata
        return adata

    def preprocess(
        self,
        adata: Optional['anndata.AnnData'] = None,
        normalize_total: bool = True,
        target_sum: float = 1e4,
        log_transform: bool = True,
        n_top_genes: int = 2000,
        regress_out: Optional[List[str]] = None,
        scale: bool = True,
        max_value: float = 10
    ) -> 'anndata.AnnData':
        """
        Preprocess data (normalization, HVG selection, scaling).

        Parameters
        ----------
        adata : AnnData, optional
            Data object
        normalize_total : bool
            Apply library size normalization
        target_sum : float
            Target sum for normalization
        log_transform : bool
            Apply log1p transformation
        n_top_genes : int
            Number of highly variable genes
        regress_out : list, optional
            Variables to regress out (e.g., ['total_counts', 'pct_counts_mt'])
        scale : bool
            Scale data to unit variance
        max_value : float
            Clip values after scaling

        Returns
        -------
        AnnData
            Preprocessed data
        """
        if adata is None:
            adata = self.adata

        # Store raw counts
        adata.layers['counts'] = adata.X.copy()

        # Normalize
        if normalize_total:
            self.sc.pp.normalize_total(adata, target_sum=target_sum)

        # Log transform
        if log_transform:
            self.sc.pp.log1p(adata)

        # Save for visualization
        adata.raw = adata

        # Highly variable genes
        self.sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor='seurat_v3' if 'counts' in adata.layers else 'seurat'
        )

        # Filter to HVGs
        adata = adata[:, adata.var.highly_variable]

        # Regress out confounders
        if regress_out:
            self.sc.pp.regress_out(adata, regress_out)

        # Scale
        if scale:
            self.sc.pp.scale(adata, max_value=max_value)

        print(f"Preprocessing complete: {adata.n_vars} HVGs selected")

        self.adata = adata
        return adata

    def reduce_dimensions(
        self,
        adata: Optional['anndata.AnnData'] = None,
        n_pcs: int = 50,
        n_neighbors: int = 15,
        use_pcs: Optional[int] = None,
        compute_umap: bool = True,
        compute_tsne: bool = False,
        min_dist: float = 0.3
    ) -> 'anndata.AnnData':
        """
        Perform dimensionality reduction (PCA, UMAP, t-SNE).

        Parameters
        ----------
        adata : AnnData, optional
            Data object
        n_pcs : int
            Number of principal components
        n_neighbors : int
            Number of neighbors for graph construction
        use_pcs : int, optional
            Number of PCs to use (if None, determined from elbow)
        compute_umap : bool
            Compute UMAP embedding
        compute_tsne : bool
            Compute t-SNE embedding
        min_dist : float
            UMAP min_dist parameter

        Returns
        -------
        AnnData
            Data with embeddings
        """
        if adata is None:
            adata = self.adata

        # PCA
        self.sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')

        # Determine optimal number of PCs
        if use_pcs is None:
            use_pcs = self._find_elbow(adata.uns['pca']['variance_ratio'])
            print(f"Using {use_pcs} PCs (elbow method)")

        # Build neighbor graph
        self.sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=use_pcs)

        # UMAP
        if compute_umap:
            self.sc.tl.umap(adata, min_dist=min_dist)

        # t-SNE
        if compute_tsne:
            self.sc.tl.tsne(adata, n_pcs=use_pcs)

        self.adata = adata
        return adata

    def _find_elbow(self, variance_ratio: np.ndarray, threshold: float = 0.9) -> int:
        """Find elbow point in variance ratio curve."""
        cumsum = np.cumsum(variance_ratio)
        n_pcs = np.argmax(cumsum >= threshold) + 1
        return max(n_pcs, 10)  # Minimum 10 PCs

    def cluster(
        self,
        adata: Optional['anndata.AnnData'] = None,
        resolution: float = 0.8,
        method: str = 'leiden',
        resolutions: Optional[List[float]] = None
    ) -> 'anndata.AnnData':
        """
        Perform clustering.

        Parameters
        ----------
        adata : AnnData, optional
            Data object
        resolution : float
            Clustering resolution
        method : str
            'leiden' or 'louvain'
        resolutions : list, optional
            Multiple resolutions to try

        Returns
        -------
        AnnData
            Data with cluster assignments
        """
        if adata is None:
            adata = self.adata

        if resolutions is None:
            resolutions = [resolution]

        for res in resolutions:
            key = f'{method}_{res}'
            if method == 'leiden':
                self.sc.tl.leiden(adata, resolution=res, key_added=key)
            else:
                self.sc.tl.louvain(adata, resolution=res, key_added=key)

            self.clustering_results[key] = {
                'n_clusters': len(adata.obs[key].unique()),
                'resolution': res
            }

        # Set default cluster column
        adata.obs[method] = adata.obs[f'{method}_{resolution}']

        print(f"Clustering complete: {len(adata.obs[method].unique())} clusters at resolution {resolution}")

        self.adata = adata
        return adata

    def find_markers(
        self,
        adata: Optional['anndata.AnnData'] = None,
        groupby: str = 'leiden',
        method: str = 'wilcoxon',
        n_genes: int = 100,
        min_logfc: float = 0.5,
        pval_cutoff: float = 0.05,
        reference: str = 'rest'
    ) -> pd.DataFrame:
        """
        Find marker genes for clusters.

        Parameters
        ----------
        adata : AnnData, optional
            Data object
        groupby : str
            Column to group by
        method : str
            Test method ('wilcoxon', 't-test', 'logreg')
        n_genes : int
            Number of genes per cluster
        min_logfc : float
            Minimum log fold change
        pval_cutoff : float
            P-value cutoff
        reference : str
            Reference group ('rest' or specific group)

        Returns
        -------
        DataFrame
            Marker genes
        """
        if adata is None:
            adata = self.adata

        # Run DE analysis
        self.sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method=method,
            pts=True,
            reference=reference
        )

        # Get results as DataFrame
        markers = self.sc.get.rank_genes_groups_df(adata, group=None)

        # Filter
        markers = markers[
            (markers['logfoldchanges'].abs() >= min_logfc) &
            (markers['pvals_adj'] < pval_cutoff)
        ]

        self.marker_genes[groupby] = markers
        self.adata = adata

        return markers

    def discover_novel_genes(
        self,
        adata: Optional['anndata.AnnData'] = None,
        cluster_key: str = 'leiden',
        reference_markers: Optional[Dict[str, List[str]]] = None,
        novelty_threshold: float = 0.8,
        min_pct: float = 0.25,
        n_top: int = 50
    ) -> Dict[str, List[Dict]]:
        """
        Discover potentially novel marker genes.

        Parameters
        ----------
        adata : AnnData, optional
            Data object
        cluster_key : str
            Cluster column
        reference_markers : dict, optional
            Known markers {cell_type: [genes]}
        novelty_threshold : float
            Specificity threshold for novelty
        min_pct : float
            Minimum percentage expressing
        n_top : int
            Top N genes to consider per cluster

        Returns
        -------
        dict
            Novel candidates per cluster
        """
        if adata is None:
            adata = self.adata

        if reference_markers is None:
            reference_markers = {}

        # Get all known markers
        known_genes = set()
        for genes in reference_markers.values():
            known_genes.update(genes)

        # Find markers if not already done
        if cluster_key not in self.marker_genes:
            self.find_markers(adata, groupby=cluster_key)

        markers_df = self.marker_genes[cluster_key]

        novel_candidates = {}

        for cluster in adata.obs[cluster_key].unique():
            cluster_markers = markers_df[markers_df['group'] == str(cluster)]
            top_genes = cluster_markers.head(n_top)['names'].tolist()

            novel = []
            for gene in top_genes:
                if gene in known_genes:
                    continue

                # Calculate specificity
                expr_in_cluster = adata[
                    adata.obs[cluster_key] == cluster, gene
                ].X.mean() if hasattr(adata[:, gene].X, 'mean') else np.mean(adata[adata.obs[cluster_key] == cluster, gene].X)

                expr_in_others = adata[
                    adata.obs[cluster_key] != cluster, gene
                ].X.mean() if hasattr(adata[:, gene].X, 'mean') else np.mean(adata[adata.obs[cluster_key] != cluster, gene].X)

                if expr_in_cluster + expr_in_others > 0:
                    specificity = expr_in_cluster / (expr_in_cluster + expr_in_others)

                    # Check percentage expressing
                    pct_expr = np.mean(
                        adata[adata.obs[cluster_key] == cluster, gene].X > 0
                    )

                    if specificity >= novelty_threshold and pct_expr >= min_pct:
                        # Get stats from marker table
                        gene_stats = cluster_markers[
                            cluster_markers['names'] == gene
                        ].iloc[0] if len(cluster_markers[cluster_markers['names'] == gene]) > 0 else None

                        novel.append({
                            'gene': gene,
                            'cluster': cluster,
                            'specificity': float(specificity),
                            'pct_expressing': float(pct_expr),
                            'expr_in_cluster': float(expr_in_cluster),
                            'expr_in_others': float(expr_in_others),
                            'logfc': float(gene_stats['logfoldchanges']) if gene_stats is not None else None,
                            'pval_adj': float(gene_stats['pvals_adj']) if gene_stats is not None else None
                        })

            novel_candidates[str(cluster)] = sorted(
                novel,
                key=lambda x: x['specificity'],
                reverse=True
            )

        return novel_candidates

    def full_pipeline(
        self,
        path: str,
        min_genes: int = 200,
        max_genes: int = 5000,
        max_mt_pct: float = 10,
        n_hvgs: int = 2000,
        resolution: float = 0.8,
        output_dir: str = './results'
    ) -> Tuple['anndata.AnnData', pd.DataFrame, Dict]:
        """
        Run complete analysis pipeline.

        Parameters
        ----------
        path : str
            Path to input data
        min_genes : int
            Minimum genes per cell
        max_genes : int
            Maximum genes per cell
        max_mt_pct : float
            Maximum mitochondrial percentage
        n_hvgs : int
            Number of highly variable genes
        resolution : float
            Clustering resolution
        output_dir : str
            Output directory

        Returns
        -------
        tuple
            (adata, markers, novel_candidates)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

        # Load data
        print("Loading data...")
        self.load_data(path)

        # QC
        print("Running QC...")
        self.quality_control(
            min_genes=min_genes,
            max_genes=max_genes,
            max_mt_pct=max_mt_pct
        )

        # Preprocess
        print("Preprocessing...")
        self.preprocess(n_top_genes=n_hvgs)

        # Dimensionality reduction
        print("Dimensionality reduction...")
        self.reduce_dimensions()

        # Clustering
        print("Clustering...")
        self.cluster(resolution=resolution)

        # Find markers
        print("Finding markers...")
        markers = self.find_markers()

        # Discover novel genes
        print("Discovering novel genes...")
        novel = self.discover_novel_genes()

        # Save results
        self.adata.write(output_dir / 'processed.h5ad')
        markers.to_csv(output_dir / 'markers.csv', index=False)

        # Save novel candidates
        import json
        with open(output_dir / 'novel_candidates.json', 'w') as f:
            json.dump(novel, f, indent=2)

        print(f"Results saved to {output_dir}")

        return self.adata, markers, novel

    def annotate_cells(
        self,
        adata: Optional['anndata.AnnData'] = None,
        marker_dict: Dict[str, List[str]] = None,
        cluster_key: str = 'leiden'
    ) -> 'anndata.AnnData':
        """
        Annotate cell types based on marker genes.

        Parameters
        ----------
        adata : AnnData, optional
            Data object
        marker_dict : dict
            {cell_type: [marker_genes]}
        cluster_key : str
            Cluster column

        Returns
        -------
        AnnData
            Annotated data
        """
        if adata is None:
            adata = self.adata

        if marker_dict is None:
            # Default markers for common cell types
            marker_dict = {
                'T cells': ['CD3D', 'CD3E', 'CD4', 'CD8A'],
                'B cells': ['MS4A1', 'CD79A', 'CD19'],
                'NK cells': ['NKG7', 'GNLY', 'NCAM1'],
                'Monocytes': ['CD14', 'LYZ', 'S100A8'],
                'Dendritic': ['FCER1A', 'CD1C'],
                'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
                'Fibroblasts': ['COL1A1', 'DCN', 'LUM']
            }

        # Score each cell type
        for cell_type, markers in marker_dict.items():
            present_markers = [m for m in markers if m in adata.var_names]
            if present_markers:
                self.sc.tl.score_genes(
                    adata,
                    present_markers,
                    score_name=f'{cell_type}_score'
                )

        self.adata = adata
        return adata


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='scRNA-seq Analysis Pipeline')
    parser.add_argument('input', help='Input file path')
    parser.add_argument('--output', '-o', default='./results', help='Output directory')
    parser.add_argument('--min-genes', type=int, default=200)
    parser.add_argument('--max-genes', type=int, default=5000)
    parser.add_argument('--max-mt', type=float, default=10)
    parser.add_argument('--n-hvgs', type=int, default=2000)
    parser.add_argument('--resolution', type=float, default=0.8)

    args = parser.parse_args()

    analyzer = ScRNAAnalyzer()
    analyzer.full_pipeline(
        args.input,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_mt_pct=args.max_mt,
        n_hvgs=args.n_hvgs,
        resolution=args.resolution,
        output_dir=args.output
    )
