# Universal Pathway Figure Templates

## Overview

어떤 pathway에도 적용 가능한 Nature Reviews 스타일 figure 템플릿 시스템

## Pathway Template Library

```python
PATHWAY_TEMPLATES = {
    # ===================
    # RECEPTOR TYROSINE KINASE PATHWAYS
    # ===================
    'RTK_signaling': {
        'name': 'Receptor Tyrosine Kinase Signaling',
        'compartments': ['extracellular', 'membrane', 'cytoplasm', 'nucleus'],
        'components': {
            'ligand': {'position': 'extracellular', 'symbol': 'circle'},
            'receptor': {'position': 'membrane', 'symbol': 'y_shape'},
            'adaptor': {'position': 'cytoplasm', 'symbol': 'rectangle'},
            'kinase': {'position': 'cytoplasm', 'symbol': 'rounded_rect'},
            'transcription_factor': {'position': 'nucleus', 'symbol': 'pentagon'}
        },
        'canonical_flow': [
            'ligand → receptor (binding)',
            'receptor → receptor (dimerization)',
            'receptor → adaptor (recruitment)',
            'adaptor → kinase (activation)',
            'kinase → transcription_factor (phosphorylation)',
            'transcription_factor → gene (transcription)'
        ],
        'examples': ['VEGF', 'EGF', 'FGF', 'PDGF', 'insulin']
    },

    'GPCR_signaling': {
        'name': 'G Protein-Coupled Receptor Signaling',
        'compartments': ['extracellular', 'membrane', 'cytoplasm'],
        'components': {
            'ligand': {'position': 'extracellular', 'symbol': 'circle'},
            'receptor': {'position': 'membrane', 'symbol': 'serpentine'},
            'g_protein': {'position': 'membrane', 'symbol': 'trimer'},
            'effector': {'position': 'cytoplasm', 'symbol': 'rectangle'},
            'second_messenger': {'position': 'cytoplasm', 'symbol': 'small_circle'}
        },
        'canonical_flow': [
            'ligand → receptor (binding)',
            'receptor → g_protein (activation)',
            'g_protein → effector (modulation)',
            'effector → second_messenger (production)'
        ]
    },

    'notch_signaling': {
        'name': 'Notch Signaling Pathway',
        'compartments': ['signal_cell', 'receiving_cell', 'nucleus'],
        'components': {
            'ligand': {'position': 'signal_cell', 'symbol': 'dll_jagged'},
            'receptor': {'position': 'receiving_cell', 'symbol': 'notch'},
            'nicd': {'position': 'cytoplasm', 'symbol': 'cleaved'},
            'transcription_complex': {'position': 'nucleus', 'symbol': 'complex'}
        },
        'canonical_flow': [
            'ligand → receptor (trans-activation)',
            'receptor → ADAM (S2 cleavage)',
            'receptor → gamma-secretase (S3 cleavage)',
            'nicd → nucleus (translocation)',
            'nicd + RBPJ + MAML → gene (transcription)'
        ],
        'special_features': ['lateral_inhibition', 'cell_cell_contact']
    },

    'wnt_signaling': {
        'name': 'Wnt/β-catenin Signaling',
        'compartments': ['extracellular', 'membrane', 'cytoplasm', 'nucleus'],
        'components': {
            'wnt': {'position': 'extracellular', 'symbol': 'circle'},
            'frizzled': {'position': 'membrane', 'symbol': 'serpentine'},
            'lrp': {'position': 'membrane', 'symbol': 'coreceptor'},
            'destruction_complex': {'position': 'cytoplasm', 'symbol': 'complex'},
            'beta_catenin': {'position': 'cytoplasm', 'symbol': 'oval'},
            'tcf_lef': {'position': 'nucleus', 'symbol': 'pentagon'}
        },
        'states': {
            'off': 'destruction_complex active → β-catenin degraded',
            'on': 'Wnt binding → destruction_complex inhibited → β-catenin stabilized'
        }
    },

    'tgfb_signaling': {
        'name': 'TGF-β/SMAD Signaling',
        'compartments': ['extracellular', 'membrane', 'cytoplasm', 'nucleus'],
        'components': {
            'ligand': {'position': 'extracellular', 'symbol': 'dimer'},
            'type2_receptor': {'position': 'membrane', 'symbol': 'receptor'},
            'type1_receptor': {'position': 'membrane', 'symbol': 'receptor'},
            'rsmad': {'position': 'cytoplasm', 'symbol': 'oval'},
            'smad4': {'position': 'cytoplasm', 'symbol': 'oval'},
            'smad_complex': {'position': 'nucleus', 'symbol': 'complex'}
        },
        'canonical_flow': [
            'ligand → type2_receptor (binding)',
            'type2_receptor → type1_receptor (phosphorylation)',
            'type1_receptor → rsmad (phosphorylation)',
            'rsmad + smad4 → complex (formation)',
            'complex → nucleus (translocation)'
        ]
    },

    'nfkb_signaling': {
        'name': 'NF-κB Signaling',
        'compartments': ['membrane', 'cytoplasm', 'nucleus'],
        'components': {
            'receptor': {'position': 'membrane', 'symbol': 'receptor'},
            'ikk': {'position': 'cytoplasm', 'symbol': 'kinase'},
            'ikb': {'position': 'cytoplasm', 'symbol': 'inhibitor'},
            'nfkb': {'position': 'cytoplasm', 'symbol': 'dimer'},
            'proteasome': {'position': 'cytoplasm', 'symbol': 'barrel'}
        },
        'canonical_flow': [
            'stimulus → receptor → IKK (activation)',
            'IKK → IκB (phosphorylation)',
            'IκB → proteasome (degradation)',
            'NF-κB → nucleus (translocation)'
        ]
    },

    'jak_stat_signaling': {
        'name': 'JAK-STAT Signaling',
        'compartments': ['extracellular', 'membrane', 'cytoplasm', 'nucleus'],
        'components': {
            'cytokine': {'position': 'extracellular', 'symbol': 'ligand'},
            'receptor': {'position': 'membrane', 'symbol': 'dimer_receptor'},
            'jak': {'position': 'membrane', 'symbol': 'kinase'},
            'stat': {'position': 'cytoplasm', 'symbol': 'monomer'},
            'stat_dimer': {'position': 'nucleus', 'symbol': 'dimer'}
        }
    },

    'hippo_signaling': {
        'name': 'Hippo Signaling Pathway',
        'compartments': ['membrane', 'cytoplasm', 'nucleus'],
        'components': {
            'mst': {'position': 'cytoplasm', 'symbol': 'kinase'},
            'lats': {'position': 'cytoplasm', 'symbol': 'kinase'},
            'yap_taz': {'position': 'cytoplasm', 'symbol': 'effector'},
            'tead': {'position': 'nucleus', 'symbol': 'tf'}
        }
    },

    'mapk_cascade': {
        'name': 'MAPK Signaling Cascade',
        'compartments': ['membrane', 'cytoplasm', 'nucleus'],
        'components': {
            'mapkkk': {'position': 'cytoplasm', 'symbol': 'kinase', 'examples': ['RAF', 'MEKK']},
            'mapkk': {'position': 'cytoplasm', 'symbol': 'kinase', 'examples': ['MEK', 'MKK']},
            'mapk': {'position': 'cytoplasm', 'symbol': 'kinase', 'examples': ['ERK', 'JNK', 'p38']}
        },
        'cascade_structure': 'Linear three-tier kinase cascade',
        'variants': {
            'ERK': ['RAF', 'MEK1/2', 'ERK1/2'],
            'JNK': ['MLK/TAK1', 'MKK4/7', 'JNK'],
            'p38': ['MLK/TAK1', 'MKK3/6', 'p38']
        }
    },

    'pi3k_akt_mtor': {
        'name': 'PI3K-AKT-mTOR Pathway',
        'compartments': ['membrane', 'cytoplasm'],
        'components': {
            'rtk': {'position': 'membrane', 'symbol': 'receptor'},
            'pi3k': {'position': 'membrane', 'symbol': 'kinase'},
            'pip3': {'position': 'membrane', 'symbol': 'lipid'},
            'pdk1': {'position': 'membrane', 'symbol': 'kinase'},
            'akt': {'position': 'cytoplasm', 'symbol': 'kinase'},
            'mtorc1': {'position': 'cytoplasm', 'symbol': 'complex'},
            'mtorc2': {'position': 'membrane', 'symbol': 'complex'}
        }
    }
}
```

## Universal Figure Generator

```python
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple


class UniversalPathwayFigure:
    """
    Generate pathway figures for any signaling pathway.

    Supports automatic layout based on pathway templates.
    """

    def __init__(
        self,
        pathway_type: str = None,
        width: float = 180,
        height: float = 140,
        style: str = 'nature_reviews'
    ):
        """
        Initialize figure generator.

        Parameters
        ----------
        pathway_type : str, optional
            Pathway template type (e.g., 'RTK_signaling', 'notch_signaling')
        width : float
            Figure width in mm
        height : float
            Figure height in mm
        style : str
            Figure style
        """
        self.pathway_type = pathway_type
        self.width = width
        self.height = height
        self.style = style

        self.template = PATHWAY_TEMPLATES.get(pathway_type, {})
        self.components = {}
        self.interactions = []

        self._init_svg()

    def _init_svg(self):
        """Initialize SVG structure."""
        self.svg = ET.Element('svg', {
            'xmlns': 'http://www.w3.org/2000/svg',
            'width': f'{self.width}mm',
            'height': f'{self.height}mm',
            'viewBox': f'0 0 {self.width * 3.78} {self.height * 3.78}'
        })

        self._add_defs()
        self._add_styles()

        # Create layers
        self.layers = {
            'background': ET.SubElement(self.svg, 'g', {'id': 'background'}),
            'compartments': ET.SubElement(self.svg, 'g', {'id': 'compartments'}),
            'interactions': ET.SubElement(self.svg, 'g', {'id': 'interactions'}),
            'molecules': ET.SubElement(self.svg, 'g', {'id': 'molecules'}),
            'labels': ET.SubElement(self.svg, 'g', {'id': 'labels'}),
            'annotations': ET.SubElement(self.svg, 'g', {'id': 'annotations'}),
            'legend': ET.SubElement(self.svg, 'g', {'id': 'legend'})
        }

    def _add_defs(self):
        """Add SVG definitions (markers, gradients)."""
        defs = ET.SubElement(self.svg, 'defs')

        # Arrow markers
        arrow = ET.SubElement(defs, 'marker', {
            'id': 'arrowhead',
            'markerWidth': '10',
            'markerHeight': '7',
            'refX': '9',
            'refY': '3.5',
            'orient': 'auto'
        })
        ET.SubElement(arrow, 'polygon', {
            'points': '0 0, 10 3.5, 0 7',
            'fill': '#333'
        })

        # T-bar for inhibition
        tbar = ET.SubElement(defs, 'marker', {
            'id': 'tbar',
            'markerWidth': '10',
            'markerHeight': '10',
            'refX': '10',
            'refY': '5',
            'orient': 'auto'
        })
        ET.SubElement(tbar, 'line', {
            'x1': '0', 'y1': '0',
            'x2': '0', 'y2': '10',
            'stroke': '#E74C3C',
            'stroke-width': '3'
        })

        # Phosphorylation marker
        phos = ET.SubElement(defs, 'marker', {
            'id': 'phospho',
            'markerWidth': '12',
            'markerHeight': '12',
            'refX': '6',
            'refY': '6',
            'orient': 'auto'
        })
        ET.SubElement(phos, 'circle', {
            'cx': '6', 'cy': '6', 'r': '5',
            'fill': '#3498DB'
        })
        ET.SubElement(phos, 'text', {
            'x': '6', 'y': '9',
            'text-anchor': 'middle',
            'fill': 'white',
            'font-size': '8',
            'font-weight': 'bold'
        }).text = 'P'

    def _add_styles(self):
        """Add CSS styles."""
        style = ET.SubElement(self.svg, 'style')
        style.text = '''
            .molecule { stroke: #333; stroke-width: 1.5; }
            .molecule-label { font-family: Arial, sans-serif; font-size: 9px; font-weight: bold; }
            .compartment-label { font-family: Arial, sans-serif; font-size: 11px; font-style: italic; fill: #666; }
            .annotation { font-family: Arial, sans-serif; font-size: 8px; fill: #E74C3C; }
            .interaction { stroke-width: 2; fill: none; }
            .activation { stroke: #27AE60; }
            .inhibition { stroke: #E74C3C; }
            .phosphorylation { stroke: #3498DB; stroke-dasharray: 5,3; }
            .translocation { stroke: #9B59B6; stroke-dasharray: 2,2; }
        '''

    def auto_layout(self, molecules: Dict[str, dict]) -> Dict[str, Tuple[float, float]]:
        """
        Automatically layout molecules based on pathway template.

        Parameters
        ----------
        molecules : dict
            {name: {'type': 'receptor', 'compartment': 'membrane', ...}}

        Returns
        -------
        dict
            {name: (x, y)}
        """
        positions = {}

        # Define compartment positions
        compartment_y = {
            'extracellular': 50,
            'membrane': 150,
            'cytoplasm': 300,
            'nucleus': 450
        }

        # Group molecules by compartment
        by_compartment = {}
        for name, info in molecules.items():
            comp = info.get('compartment', 'cytoplasm')
            if comp not in by_compartment:
                by_compartment[comp] = []
            by_compartment[comp].append(name)

        # Layout within each compartment
        for comp, mol_list in by_compartment.items():
            y = compartment_y.get(comp, 300)
            n = len(mol_list)
            spacing = (self.width * 3.78 - 100) / (n + 1)

            for i, mol in enumerate(mol_list):
                x = 50 + (i + 1) * spacing
                positions[mol] = (x, y)

        return positions

    def add_compartment(
        self,
        name: str,
        y_start: float,
        height: float,
        color: str = '#f5f5f5'
    ):
        """Add cellular compartment."""
        comp = ET.SubElement(self.layers['compartments'], 'g', {'class': 'compartment'})

        # Background
        ET.SubElement(comp, 'rect', {
            'x': '10',
            'y': str(y_start),
            'width': str(self.width * 3.78 - 20),
            'height': str(height),
            'fill': color,
            'stroke': '#ddd',
            'stroke-width': '1',
            'rx': '5'
        })

        # Label
        ET.SubElement(comp, 'text', {
            'x': '20',
            'y': str(y_start + 20),
            'class': 'compartment-label'
        }).text = name

    def add_molecule(
        self,
        name: str,
        mol_type: str,
        x: float,
        y: float,
        color: str = None,
        state: str = None
    ):
        """
        Add molecule to figure.

        Parameters
        ----------
        name : str
            Molecule name
        mol_type : str
            Type: 'receptor', 'kinase', 'transcription_factor', 'ligand', etc.
        x, y : float
            Position
        color : str, optional
            Override color
        state : str, optional
            State annotation (e.g., 'phosphorylated', 'active')
        """
        if color is None:
            color = self._get_default_color(mol_type)

        mol_group = ET.SubElement(self.layers['molecules'], 'g', {
            'class': f'molecule {mol_type}',
            'transform': f'translate({x}, {y})'
        })

        # Draw shape based on type
        if mol_type == 'receptor':
            self._draw_receptor(mol_group, color)
        elif mol_type == 'kinase':
            self._draw_kinase(mol_group, color)
        elif mol_type == 'transcription_factor':
            self._draw_tf(mol_group, color)
        elif mol_type == 'ligand':
            self._draw_ligand(mol_group, color)
        elif mol_type == 'complex':
            self._draw_complex(mol_group, color)
        elif mol_type == 'adaptor':
            self._draw_adaptor(mol_group, color)
        else:
            self._draw_generic(mol_group, color)

        # Add label
        ET.SubElement(mol_group, 'text', {
            'y': '30',
            'text-anchor': 'middle',
            'class': 'molecule-label'
        }).text = name

        # Add state annotation
        if state:
            ET.SubElement(mol_group, 'text', {
                'y': '42',
                'text-anchor': 'middle',
                'class': 'annotation',
                'font-size': '7'
            }).text = f'({state})'

        self.components[name] = {'x': x, 'y': y, 'type': mol_type}

    def _get_default_color(self, mol_type: str) -> str:
        """Get default color for molecule type."""
        colors = {
            'receptor': '#3498DB',
            'kinase': '#9B59B6',
            'transcription_factor': '#E67E22',
            'ligand': '#27AE60',
            'adaptor': '#1ABC9C',
            'complex': '#34495E',
            'inhibitor': '#E74C3C'
        }
        return colors.get(mol_type, '#666666')

    def _draw_receptor(self, parent, color):
        """Draw receptor (Y-shape)."""
        # Extracellular domain
        ET.SubElement(parent, 'path', {
            'd': 'M-6,-15 L-6,-5 M6,-15 L6,-5',
            'stroke': color, 'stroke-width': '3', 'fill': 'none'
        })
        ET.SubElement(parent, 'circle', {'cx': '-6', 'cy': '-15', 'r': '5', 'fill': color})
        ET.SubElement(parent, 'circle', {'cx': '6', 'cy': '-15', 'r': '5', 'fill': color})

        # Transmembrane
        ET.SubElement(parent, 'rect', {
            'x': '-8', 'y': '-5', 'width': '16', 'height': '8',
            'fill': color, 'rx': '2'
        })

        # Intracellular
        ET.SubElement(parent, 'path', {
            'd': 'M0,3 L0,15',
            'stroke': color, 'stroke-width': '4', 'fill': 'none'
        })

    def _draw_kinase(self, parent, color):
        """Draw kinase (rounded rectangle with P)."""
        ET.SubElement(parent, 'rect', {
            'x': '-15', 'y': '-10', 'width': '30', 'height': '20',
            'fill': color, 'rx': '5'
        })
        ET.SubElement(parent, 'text', {
            'x': '0', 'y': '5', 'text-anchor': 'middle',
            'fill': 'white', 'font-size': '12', 'font-weight': 'bold'
        }).text = 'P'

    def _draw_tf(self, parent, color):
        """Draw transcription factor (pentagon)."""
        ET.SubElement(parent, 'polygon', {
            'points': '-12,8 -12,-4 0,-12 12,-4 12,8',
            'fill': color
        })

    def _draw_ligand(self, parent, color):
        """Draw ligand (circle)."""
        ET.SubElement(parent, 'circle', {'r': '10', 'fill': color})

    def _draw_complex(self, parent, color):
        """Draw protein complex (overlapping circles)."""
        ET.SubElement(parent, 'circle', {'cx': '-6', 'r': '10', 'fill': color, 'opacity': '0.8'})
        ET.SubElement(parent, 'circle', {'cx': '6', 'r': '10', 'fill': '#2C3E50', 'opacity': '0.8'})

    def _draw_adaptor(self, parent, color):
        """Draw adaptor protein (rectangle)."""
        ET.SubElement(parent, 'rect', {
            'x': '-12', 'y': '-8', 'width': '24', 'height': '16',
            'fill': color
        })

    def _draw_generic(self, parent, color):
        """Draw generic molecule (oval)."""
        ET.SubElement(parent, 'ellipse', {
            'rx': '15', 'ry': '10', 'fill': color
        })

    def add_interaction(
        self,
        from_mol: str,
        to_mol: str,
        interaction_type: str = 'activation',
        label: str = None
    ):
        """
        Add interaction between molecules.

        Parameters
        ----------
        from_mol, to_mol : str
            Molecule names
        interaction_type : str
            'activation', 'inhibition', 'phosphorylation', 'translocation', 'binding'
        label : str, optional
            Interaction label
        """
        if from_mol not in self.components or to_mol not in self.components:
            return

        from_pos = self.components[from_mol]
        to_pos = self.components[to_mol]

        # Calculate line endpoints
        x1, y1 = from_pos['x'], from_pos['y'] + 15
        x2, y2 = to_pos['x'], to_pos['y'] - 15

        # Create interaction line
        line_attrs = {
            'x1': str(x1), 'y1': str(y1),
            'x2': str(x2), 'y2': str(y2),
            'class': f'interaction {interaction_type}'
        }

        if interaction_type == 'activation':
            line_attrs['marker-end'] = 'url(#arrowhead)'
            line_attrs['stroke'] = '#27AE60'
        elif interaction_type == 'inhibition':
            line_attrs['marker-end'] = 'url(#tbar)'
            line_attrs['stroke'] = '#E74C3C'
        elif interaction_type == 'phosphorylation':
            line_attrs['marker-end'] = 'url(#phospho)'
            line_attrs['stroke'] = '#3498DB'
            line_attrs['stroke-dasharray'] = '5,3'

        ET.SubElement(self.layers['interactions'], 'line', line_attrs)

        # Add label if provided
        if label:
            mid_x = (x1 + x2) / 2
            mid_y = (y1 + y2) / 2
            ET.SubElement(self.layers['labels'], 'text', {
                'x': str(mid_x + 10),
                'y': str(mid_y),
                'class': 'annotation',
                'font-size': '7'
            }).text = label

        self.interactions.append({
            'from': from_mol,
            'to': to_mol,
            'type': interaction_type
        })

    def add_annotation(
        self,
        text: str,
        x: float,
        y: float,
        annotation_type: str = 'note'
    ):
        """Add annotation to figure."""
        style_map = {
            'note': {'fill': '#666', 'font-style': 'italic'},
            'ko_effect': {'fill': '#E74C3C', 'font-weight': 'bold'},
            'drug_target': {'fill': '#9B59B6', 'font-weight': 'bold'},
            'phenotype': {'fill': '#27AE60'}
        }

        attrs = {
            'x': str(x),
            'y': str(y),
            'class': 'annotation',
            'font-size': '8'
        }
        attrs.update(style_map.get(annotation_type, {}))

        ET.SubElement(self.layers['annotations'], 'text', attrs).text = text

    def add_legend(self, x: float = None, y: float = None):
        """Add figure legend."""
        if x is None:
            x = self.width * 3.78 - 100
        if y is None:
            y = 20

        legend = ET.SubElement(self.layers['legend'], 'g', {
            'transform': f'translate({x}, {y})'
        })

        ET.SubElement(legend, 'text', {
            'y': '0', 'font-weight': 'bold', 'font-size': '10'
        }).text = 'Legend'

        items = [
            ('Activation', '#27AE60', 'solid'),
            ('Inhibition', '#E74C3C', 'tbar'),
            ('Phosphorylation', '#3498DB', 'dashed')
        ]

        for i, (label, color, style) in enumerate(items):
            y_pos = 20 + i * 18

            line_attrs = {
                'x1': '0', 'y1': str(y_pos),
                'x2': '25', 'y2': str(y_pos),
                'stroke': color, 'stroke-width': '2'
            }
            if style == 'dashed':
                line_attrs['stroke-dasharray'] = '5,3'

            ET.SubElement(legend, 'line', line_attrs)
            ET.SubElement(legend, 'text', {
                'x': '30', 'y': str(y_pos + 4),
                'font-size': '8'
            }).text = label

    def save(self, filename: str):
        """Save figure to file."""
        tree = ET.ElementTree(self.svg)
        with open(filename, 'wb') as f:
            f.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
            tree.write(f, encoding='unicode')
        print(f"Figure saved to {filename}")

    def to_string(self) -> str:
        """Return SVG as string."""
        return ET.tostring(self.svg, encoding='unicode')


def quick_pathway_figure(
    pathway: str,
    molecules: dict,
    interactions: list,
    output: str = 'pathway.svg'
) -> str:
    """
    Quick function to generate pathway figure.

    Parameters
    ----------
    pathway : str
        Pathway name
    molecules : dict
        {name: {'type': 'receptor', 'compartment': 'membrane'}}
    interactions : list
        [('A', 'B', 'activation'), ...]
    output : str
        Output filename

    Returns
    -------
    str
        SVG content
    """
    fig = UniversalPathwayFigure(pathway_type='RTK_signaling')

    # Auto-layout molecules
    positions = fig.auto_layout(molecules)

    # Add compartments
    fig.add_compartment('Extracellular', 30, 80, '#E8F6F3')
    fig.add_compartment('Membrane', 110, 60, '#FCF3CF')
    fig.add_compartment('Cytoplasm', 170, 150, '#F5F5F5')
    fig.add_compartment('Nucleus', 320, 100, '#EBF5FB')

    # Add molecules
    for name, info in molecules.items():
        x, y = positions[name]
        fig.add_molecule(name, info['type'], x, y)

    # Add interactions
    for from_mol, to_mol, int_type in interactions:
        fig.add_interaction(from_mol, to_mol, int_type)

    fig.add_legend()
    fig.save(output)

    return fig.to_string()
```

## Usage Examples

```python
# Example 1: VEGF signaling
molecules = {
    'VEGF-A': {'type': 'ligand', 'compartment': 'extracellular'},
    'VEGFR2': {'type': 'receptor', 'compartment': 'membrane'},
    'PLCγ': {'type': 'kinase', 'compartment': 'cytoplasm'},
    'ERK': {'type': 'kinase', 'compartment': 'cytoplasm'},
    'ETS': {'type': 'transcription_factor', 'compartment': 'nucleus'}
}

interactions = [
    ('VEGF-A', 'VEGFR2', 'activation'),
    ('VEGFR2', 'PLCγ', 'phosphorylation'),
    ('PLCγ', 'ERK', 'activation'),
    ('ERK', 'ETS', 'phosphorylation')
]

svg = quick_pathway_figure('VEGF', molecules, interactions, 'vegf_pathway.svg')


# Example 2: Any RTK pathway
def create_rtk_figure(ligand, receptor, downstream):
    """Create figure for any RTK pathway."""
    molecules = {
        ligand: {'type': 'ligand', 'compartment': 'extracellular'},
        receptor: {'type': 'receptor', 'compartment': 'membrane'}
    }

    for name, info in downstream.items():
        molecules[name] = info

    # Auto-generate interactions based on order
    # ...

    return quick_pathway_figure('RTK', molecules, [], f'{receptor}_pathway.svg')
```
