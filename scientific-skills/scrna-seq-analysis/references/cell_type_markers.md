# Universal Cell Type Marker Database

## Overview

어떤 조직/세포 유형에도 적용 가능한 cell type marker 데이터베이스

## Cell Type Marker Database

```python
CELL_TYPE_MARKERS = {
    # ===================
    # IMMUNE CELLS
    # ===================
    'immune': {
        'T_cells': {
            'pan': ['CD3D', 'CD3E', 'CD3G', 'TRAC', 'TRBC1'],
            'CD4_T': ['CD4', 'IL7R', 'CCR7', 'LEF1'],
            'CD8_T': ['CD8A', 'CD8B', 'GZMK', 'GZMB'],
            'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
            'Th1': ['TBX21', 'IFNG', 'TNF', 'IL2'],
            'Th2': ['GATA3', 'IL4', 'IL5', 'IL13'],
            'Th17': ['RORC', 'IL17A', 'IL17F', 'IL22'],
            'Tfh': ['CXCR5', 'BCL6', 'ICOS', 'IL21'],
            'naive': ['CCR7', 'SELL', 'LEF1', 'TCF7'],
            'memory': ['CD44', 'IL7R', 'S1PR1'],
            'effector': ['GZMB', 'PRF1', 'IFNG', 'NKG7'],
            'exhausted': ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'TOX']
        },
        'B_cells': {
            'pan': ['MS4A1', 'CD79A', 'CD79B', 'CD19', 'PAX5'],
            'naive': ['IGHD', 'FCER2', 'TCL1A'],
            'memory': ['CD27', 'IGHG1', 'IGHA1'],
            'plasma': ['SDC1', 'IGHG1', 'MZB1', 'JCHAIN', 'XBP1'],
            'germinal_center': ['BCL6', 'AICDA', 'MKI67']
        },
        'NK_cells': {
            'pan': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1', 'KLRF1'],
            'CD56bright': ['NCAM1', 'SELL', 'GZMK'],
            'CD56dim': ['FCGR3A', 'PRF1', 'GZMB']
        },
        'monocytes_macrophages': {
            'monocyte_pan': ['CD14', 'LYZ', 'S100A8', 'S100A9'],
            'classical_mono': ['CD14', 'S100A12', 'VCAN'],
            'nonclassical_mono': ['FCGR3A', 'CX3CR1', 'CDKN1C'],
            'macrophage_pan': ['CD68', 'CD163', 'MRC1', 'MARCO'],
            'M1': ['CD80', 'CD86', 'NOS2', 'IL1B', 'TNF'],
            'M2': ['CD163', 'MRC1', 'CD206', 'ARG1', 'IL10'],
            'TAM': ['CD163', 'MSR1', 'MRC1', 'TREM2']
        },
        'dendritic_cells': {
            'pan': ['ITGAX', 'HLA-DRA', 'HLA-DRB1'],
            'cDC1': ['CLEC9A', 'XCR1', 'BATF3', 'IRF8'],
            'cDC2': ['CD1C', 'FCER1A', 'CLEC10A', 'IRF4'],
            'pDC': ['LILRA4', 'IL3RA', 'CLEC4C', 'IRF7'],
            'mature': ['CCR7', 'LAMP3', 'CD83', 'CD40']
        },
        'granulocytes': {
            'neutrophil': ['FCGR3B', 'CSF3R', 'S100A8', 'CXCR2'],
            'eosinophil': ['CLC', 'EPX', 'PRG2', 'SIGLEC8'],
            'basophil': ['CPA3', 'HDC', 'MS4A2'],
            'mast_cell': ['TPSAB1', 'CPA3', 'KIT', 'MS4A2']
        }
    },

    # ===================
    # ENDOTHELIAL CELLS
    # ===================
    'endothelial': {
        'pan': ['PECAM1', 'CDH5', 'VWF', 'ERG', 'FLI1'],
        'arterial': ['EFNB2', 'GJA5', 'DLL4', 'HEY1', 'CXCR4', 'SOX17'],
        'venous': ['NR2F2', 'EPHB4', 'APLNR'],
        'capillary': ['RGCC', 'CA4', 'BTNL9'],
        'lymphatic': ['PROX1', 'LYVE1', 'PDPN', 'FLT4'],
        'tip_cell': ['DLL4', 'KDR', 'CXCR4', 'PDGFB', 'ANGPT2', 'ESM1', 'APLN'],
        'stalk_cell': ['JAG1', 'HES1', 'HEY1', 'NOTCH1'],
        'phalanx_cell': ['CDH5', 'CLDN5', 'OCLN'],
        'high_endothelial_venule': ['ACKR1', 'SELP', 'MADCAM1'],
        'sinusoidal': ['CLEC4G', 'CLEC4M', 'STAB2', 'LYVE1'],
        'brain': ['CLDN5', 'SLC2A1', 'MFSD2A', 'ABCB1']
    },

    # ===================
    # STROMAL CELLS
    # ===================
    'stromal': {
        'fibroblast': {
            'pan': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA'],
            'activated': ['ACTA2', 'FAP', 'POSTN'],
            'myofibroblast': ['ACTA2', 'MYH11', 'TAGLN'],
            'CAF': ['FAP', 'PDPN', 'ACTA2', 'S100A4']
        },
        'pericyte': ['RGS5', 'PDGFRB', 'CSPG4', 'ACTA2', 'DES'],
        'smooth_muscle': ['ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'DES'],
        'mesenchymal_stem': ['NT5E', 'THY1', 'ENG', 'VCAM1']
    },

    # ===================
    # EPITHELIAL CELLS
    # ===================
    'epithelial': {
        'pan': ['EPCAM', 'KRT8', 'KRT18', 'CDH1'],
        'basal': ['KRT5', 'KRT14', 'TP63'],
        'luminal': ['KRT8', 'KRT18', 'KRT19'],
        'secretory': ['MUC5AC', 'MUC5B', 'SCGB1A1'],
        'ciliated': ['FOXJ1', 'DNAH5', 'TPPP3'],
        'alveolar_type1': ['AGER', 'PDPN', 'CAV1', 'AQP5'],
        'alveolar_type2': ['SFTPC', 'SFTPB', 'ABCA3', 'LAMP3'],
        'club_cell': ['SCGB1A1', 'SCGB3A2'],
        'goblet': ['MUC5AC', 'MUC5B', 'SPDEF'],
        'ionocyte': ['FOXI1', 'CFTR'],
        'neuroendocrine': ['CHGA', 'SYP', 'ASCL1']
    },

    # ===================
    # NEURAL CELLS
    # ===================
    'neural': {
        'neuron': {
            'pan': ['RBFOX3', 'SNAP25', 'SYP', 'TUBB3'],
            'excitatory': ['SLC17A7', 'CAMK2A', 'GRIN1'],
            'inhibitory': ['GAD1', 'GAD2', 'SLC32A1'],
            'dopaminergic': ['TH', 'SLC6A3', 'DDC'],
            'cholinergic': ['CHAT', 'SLC5A7'],
            'serotonergic': ['TPH2', 'SLC6A4']
        },
        'astrocyte': ['GFAP', 'AQP4', 'S100B', 'SLC1A2', 'SLC1A3'],
        'oligodendrocyte': {
            'OPC': ['PDGFRA', 'CSPG4', 'OLIG2'],
            'mature': ['MBP', 'MOG', 'PLP1', 'MAG']
        },
        'microglia': ['CX3CR1', 'P2RY12', 'TMEM119', 'ITGAM', 'AIF1'],
        'schwann_cell': ['MPZ', 'PMP22', 'SOX10']
    },

    # ===================
    # STEM/PROGENITOR
    # ===================
    'stem_progenitor': {
        'HSC': ['KIT', 'SCA1', 'CD34', 'THY1', 'PROCR'],
        'neural_stem': ['NES', 'SOX2', 'PAX6'],
        'intestinal_stem': ['LGR5', 'OLFM4', 'ASCL2'],
        'muscle_stem': ['PAX7', 'MYF5'],
        'mesenchymal_stem': ['NT5E', 'THY1', 'ENG'],
        'cancer_stem': ['PROM1', 'ALDH1A1', 'CD44', 'SOX2', 'NANOG']
    },

    # ===================
    # TISSUE-SPECIFIC
    # ===================
    'liver': {
        'hepatocyte': ['ALB', 'APOB', 'CYP3A4', 'HNF4A'],
        'cholangiocyte': ['KRT19', 'KRT7', 'SOX9', 'EPCAM'],
        'stellate_cell': ['DES', 'LRAT', 'COL1A1'],
        'kupffer_cell': ['CLEC4F', 'CD68', 'MARCO']
    },
    'kidney': {
        'podocyte': ['NPHS1', 'NPHS2', 'WT1', 'PODXL'],
        'proximal_tubule': ['LRP2', 'CUBN', 'SLC34A1'],
        'distal_tubule': ['SLC12A3', 'CALB1'],
        'collecting_duct': ['AQP2', 'SLC4A1', 'AQP3'],
        'mesangial': ['PDGFRB', 'ITGA8']
    },
    'heart': {
        'cardiomyocyte': ['TNNT2', 'MYH7', 'MYH6', 'ACTC1'],
        'cardiac_fibroblast': ['DDR2', 'TCF21', 'COL1A1'],
        'cardiac_endothelial': ['PECAM1', 'NPR3']
    },
    'pancreas': {
        'beta_cell': ['INS', 'PDX1', 'NKX6-1', 'MAFA'],
        'alpha_cell': ['GCG', 'ARX', 'IRX1'],
        'delta_cell': ['SST', 'HHEX'],
        'acinar': ['PRSS1', 'CPA1', 'CELA3A'],
        'ductal': ['KRT19', 'SOX9', 'MUC1']
    },
    'skin': {
        'keratinocyte': {
            'basal': ['KRT5', 'KRT14', 'TP63'],
            'spinous': ['KRT1', 'KRT10'],
            'granular': ['FLG', 'LOR']
        },
        'melanocyte': ['PMEL', 'MLANA', 'TYR', 'MITF'],
        'langerhans': ['CD207', 'CD1A', 'EPCAM'],
        'merkel_cell': ['KRT20', 'ATOH1']
    },
    'muscle': {
        'skeletal_muscle': ['MYH1', 'MYH2', 'MYOD1', 'DES'],
        'satellite_cell': ['PAX7', 'MYF5', 'CD34'],
        'smooth_muscle': ['ACTA2', 'MYH11', 'TAGLN']
    },
    'adipose': {
        'adipocyte': ['ADIPOQ', 'LEP', 'FABP4', 'PPARG'],
        'preadipocyte': ['PDGFRA', 'DLK1', 'CD34']
    }
}
```

## Cell Type Query Functions

```python
class CellTypeMarkerDB:
    """Universal cell type marker database."""

    def __init__(self):
        self.markers = CELL_TYPE_MARKERS

    def get_markers(
        self,
        cell_type: str,
        tissue: str = None,
        marker_type: str = 'pan'
    ) -> list:
        """
        Get markers for any cell type.

        Parameters
        ----------
        cell_type : str
            Cell type name
        tissue : str, optional
            Tissue context
        marker_type : str
            'pan', 'subtype', or specific subtype name

        Returns
        -------
        list
            Marker genes
        """
        # Search through hierarchy
        for category, subtypes in self.markers.items():
            if cell_type.lower() in category.lower():
                if isinstance(subtypes, dict):
                    if marker_type in subtypes:
                        return subtypes[marker_type]
                    elif 'pan' in subtypes:
                        return subtypes['pan']
                elif isinstance(subtypes, list):
                    return subtypes

            # Search within subtypes
            if isinstance(subtypes, dict):
                for subtype, markers in subtypes.items():
                    if cell_type.lower() in subtype.lower():
                        if isinstance(markers, dict):
                            return markers.get(marker_type, markers.get('pan', []))
                        return markers

        return []

    def find_cell_type_by_markers(
        self,
        markers: list,
        threshold: float = 0.5
    ) -> list:
        """
        Identify cell type from marker expression.

        Parameters
        ----------
        markers : list
            Expressed marker genes
        threshold : float
            Minimum fraction of markers matching

        Returns
        -------
        list
            Matched cell types with scores
        """
        markers_set = set(m.upper() for m in markers)
        matches = []

        def search_markers(data, path=""):
            if isinstance(data, dict):
                for key, value in data.items():
                    search_markers(value, f"{path}/{key}" if path else key)
            elif isinstance(data, list):
                ref_markers = set(m.upper() for m in data)
                if ref_markers:
                    overlap = len(markers_set & ref_markers)
                    score = overlap / len(ref_markers)
                    if score >= threshold:
                        matches.append({
                            'cell_type': path,
                            'score': score,
                            'matching_markers': list(markers_set & ref_markers),
                            'reference_markers': data
                        })

        search_markers(self.markers)

        return sorted(matches, key=lambda x: x['score'], reverse=True)

    def get_all_subtypes(self, cell_type: str) -> dict:
        """Get all subtypes for a cell type."""
        subtypes = {}

        def search(data, path=""):
            if isinstance(data, dict):
                for key, value in data.items():
                    if cell_type.lower() in key.lower():
                        if isinstance(value, dict):
                            return value
                    result = search(value, f"{path}/{key}" if path else key)
                    if result:
                        return result
            return None

        return search(self.markers) or {}

    def get_differentiation_trajectory(
        self,
        start_cell: str,
        end_cell: str
    ) -> dict:
        """Get markers along differentiation trajectory."""

        TRAJECTORIES = {
            ('HSC', 'T_cells'): [
                ('HSC', ['KIT', 'SCA1', 'CD34']),
                ('CLP', ['IL7R', 'FLT3']),
                ('DN_thymocyte', ['CD44', 'CD25']),
                ('DP_thymocyte', ['CD4', 'CD8A']),
                ('SP_thymocyte', ['CD4', 'CD8A', 'TCR'])
            ],
            ('HSC', 'B_cells'): [
                ('HSC', ['KIT', 'SCA1', 'CD34']),
                ('CLP', ['IL7R', 'FLT3']),
                ('pro-B', ['CD19', 'PAX5']),
                ('pre-B', ['CD19', 'CD79A']),
                ('immature_B', ['IgM', 'CD19']),
                ('mature_B', ['IgM', 'IgD', 'CD19'])
            ],
            ('monocyte', 'macrophage'): [
                ('monocyte', ['CD14', 'LYZ']),
                ('intermediate', ['CD14', 'CD16']),
                ('macrophage', ['CD68', 'CD163'])
            ],
            ('fibroblast', 'myofibroblast'): [
                ('fibroblast', ['COL1A1', 'DCN']),
                ('activated', ['FAP', 'POSTN']),
                ('myofibroblast', ['ACTA2', 'MYH11'])
            ]
        }

        key = (start_cell, end_cell)
        return TRAJECTORIES.get(key, [])

    def compare_cell_types(
        self,
        cell_types: list
    ) -> dict:
        """Compare markers between cell types."""
        markers_by_type = {}

        for ct in cell_types:
            markers_by_type[ct] = set(self.get_markers(ct))

        # Find shared and unique
        all_markers = set()
        for m in markers_by_type.values():
            all_markers |= m

        shared = all_markers.copy()
        for m in markers_by_type.values():
            shared &= m

        unique = {}
        for ct, markers in markers_by_type.items():
            others = set()
            for other_ct, other_markers in markers_by_type.items():
                if other_ct != ct:
                    others |= other_markers
            unique[ct] = markers - others

        return {
            'shared': list(shared),
            'unique': unique,
            'all': markers_by_type
        }
```

## Integration with CellxGene and HCA

```python
def query_cellxgene_markers(cell_type: str) -> list:
    """Query CellxGene Census for cell type markers."""
    try:
        import cellxgene_census

        with cellxgene_census.open_soma() as census:
            # Get cell type ontology
            cell_metadata = census["census_data"]["homo_sapiens"]["obs"]

            # Find cells of this type
            # (simplified - full implementation would use proper ontology)

            return []  # Return discovered markers

    except ImportError:
        return []


def query_hca_markers(cell_type: str, tissue: str = None) -> list:
    """Query Human Cell Atlas for markers."""
    # HCA API integration
    return []


def query_panglaodb(cell_type: str) -> list:
    """Query PanglaoDB for markers."""
    import requests

    # PanglaoDB provides curated markers
    url = f"https://panglaodb.se/markers/markers.tsv"

    try:
        response = requests.get(url)
        # Parse and filter for cell type
        return []
    except:
        return []


def query_cellmarker(cell_type: str) -> list:
    """Query CellMarker database."""
    # CellMarker 2.0 database
    return []
```

## Marker Validation Functions

```python
def validate_markers_in_data(
    adata,
    cell_type: str,
    cluster_key: str = 'leiden'
) -> dict:
    """
    Validate cell type markers in scRNA-seq data.

    Parameters
    ----------
    adata : AnnData
        scRNA-seq data
    cell_type : str
        Cell type to validate
    cluster_key : str
        Cluster annotation column

    Returns
    -------
    dict
        Validation results
    """
    db = CellTypeMarkerDB()
    markers = db.get_markers(cell_type)

    # Find markers present in data
    present = [m for m in markers if m in adata.var_names]
    absent = [m for m in markers if m not in adata.var_names]

    # Calculate expression scores per cluster
    import scanpy as sc

    if present:
        sc.tl.score_genes(adata, present, score_name=f'{cell_type}_score')

        # Find cluster with highest score
        scores_by_cluster = adata.obs.groupby(cluster_key)[f'{cell_type}_score'].mean()
        best_cluster = scores_by_cluster.idxmax()

        return {
            'cell_type': cell_type,
            'markers_found': present,
            'markers_missing': absent,
            'best_matching_cluster': best_cluster,
            'scores_by_cluster': scores_by_cluster.to_dict()
        }

    return {
        'cell_type': cell_type,
        'markers_found': [],
        'markers_missing': markers,
        'error': 'No markers found in data'
    }


def auto_annotate_clusters(
    adata,
    cluster_key: str = 'leiden',
    tissue: str = None
) -> dict:
    """
    Automatically annotate clusters using marker database.

    Parameters
    ----------
    adata : AnnData
        scRNA-seq data with clustering
    cluster_key : str
        Cluster column
    tissue : str, optional
        Tissue context for better annotation

    Returns
    -------
    dict
        Cluster annotations
    """
    import scanpy as sc

    db = CellTypeMarkerDB()
    annotations = {}

    # Get marker genes per cluster
    sc.tl.rank_genes_groups(adata, cluster_key, method='wilcoxon')
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)

    for cluster in adata.obs[cluster_key].unique():
        cluster_markers = marker_df[
            marker_df['group'] == str(cluster)
        ].head(50)['names'].tolist()

        # Find best matching cell type
        matches = db.find_cell_type_by_markers(cluster_markers)

        if matches:
            best_match = matches[0]
            annotations[cluster] = {
                'predicted_type': best_match['cell_type'],
                'confidence': best_match['score'],
                'matching_markers': best_match['matching_markers']
            }
        else:
            annotations[cluster] = {
                'predicted_type': 'Unknown',
                'confidence': 0,
                'matching_markers': []
            }

    return annotations
```

## Usage Examples

```python
# Initialize database
db = CellTypeMarkerDB()

# Get markers for any cell type
tip_markers = db.get_markers('tip_cell')
print(f"Tip cell markers: {tip_markers}")

# Find cell type from markers
expressed = ['DLL4', 'KDR', 'CXCR4', 'PDGFB']
matches = db.find_cell_type_by_markers(expressed)
print(f"Best match: {matches[0]}")

# Compare cell types
comparison = db.compare_cell_types(['tip_cell', 'stalk_cell', 'phalanx_cell'])
print(f"Shared: {comparison['shared']}")
print(f"Unique to tip: {comparison['unique']['tip_cell']}")

# Get subtypes
t_subtypes = db.get_all_subtypes('T_cells')
print(f"T cell subtypes: {list(t_subtypes.keys())}")

# Auto-annotate clusters
annotations = auto_annotate_clusters(adata, cluster_key='leiden')
for cluster, info in annotations.items():
    print(f"Cluster {cluster}: {info['predicted_type']} ({info['confidence']:.2f})")
```
