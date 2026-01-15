#!/usr/bin/env python3
"""
Nature Reviews Figure Generator

Generate publication-quality pathway schematics and figures
following Nature Reviews style guidelines.
"""

import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from pathlib import Path
import math


# Color palettes
COLOR_PALETTES = {
    'signaling': {
        'receptor': '#2E86AB',
        'kinase': '#A23B72',
        'transcription_factor': '#F18F01',
        'effector': '#C73E1D',
        'adaptor': '#3B1F2B',
        'ligand': '#2ECC71'
    },
    'cell_types': {
        'tip_cell': '#E63946',
        'stalk_cell': '#457B9D',
        'phalanx_cell': '#2A9D8F',
        'pericyte': '#E9C46A',
        'smooth_muscle': '#264653',
        'endothelial': '#1D3557'
    },
    'interactions': {
        'activation': '#2ECC71',
        'inhibition': '#E74C3C',
        'phosphorylation': '#3498DB',
        'transcription': '#9B59B6',
        'translocation': '#F39C12',
        'binding': '#34495E'
    },
    'colorblind_safe': {
        'blue': '#0072B2',
        'orange': '#E69F00',
        'green': '#009E73',
        'yellow': '#F0E442',
        'sky_blue': '#56B4E9',
        'vermillion': '#D55E00',
        'purple': '#CC79A7',
        'black': '#000000'
    }
}


@dataclass
class Molecule:
    """Represents a molecule in the pathway."""
    name: str
    mol_type: str  # receptor, kinase, tf, ligand, gene, complex
    x: float
    y: float
    color: Optional[str] = None
    label: Optional[str] = None
    size: float = 1.0


@dataclass
class Interaction:
    """Represents an interaction between molecules."""
    source: str
    target: str
    interaction_type: str  # activation, inhibition, phosphorylation, etc.
    label: Optional[str] = None
    curved: bool = False
    bidirectional: bool = False


@dataclass
class Annotation:
    """Represents an annotation on the figure."""
    text: str
    x: float
    y: float
    style: str = 'note'  # note, ko_effect, inhibitor, highlight


class NatureReviewsFigure:
    """
    Generate Nature Reviews-style pathway figures.

    Features:
    - Standard molecular symbols
    - Consistent color schemes
    - Publication-quality output
    - Multiple export formats
    """

    def __init__(
        self,
        width: float = 180,
        height: float = 120,
        style: str = 'nature_reviews',
        title: str = '',
        resolution: int = 300
    ):
        """
        Initialize figure.

        Parameters
        ----------
        width : float
            Figure width in mm
        height : float
            Figure height in mm
        style : str
            Style template
        title : str
            Figure title
        resolution : int
            Output resolution in dpi
        """
        self.width = width
        self.height = height
        self.style = style
        self.title = title
        self.resolution = resolution

        # Convert mm to pixels (assuming 96 dpi for SVG)
        self.px_width = width * 3.78
        self.px_height = height * 3.78

        # Storage
        self.molecules: Dict[str, Molecule] = {}
        self.interactions: List[Interaction] = []
        self.annotations: List[Annotation] = []
        self.compartments: List[Dict] = []

        # Initialize SVG
        self._init_svg()

    def _init_svg(self):
        """Initialize SVG document."""
        self.svg = ET.Element('svg', {
            'xmlns': 'http://www.w3.org/2000/svg',
            'xmlns:xlink': 'http://www.w3.org/1999/xlink',
            'width': f'{self.width}mm',
            'height': f'{self.height}mm',
            'viewBox': f'0 0 {self.px_width} {self.px_height}'
        })

        # Add definitions (markers, gradients, etc.)
        self._add_defs()

        # Add styles
        self._add_styles()

        # Create layer groups
        self.layers = {
            'background': ET.SubElement(self.svg, 'g', {'id': 'layer-background'}),
            'compartments': ET.SubElement(self.svg, 'g', {'id': 'layer-compartments'}),
            'interactions': ET.SubElement(self.svg, 'g', {'id': 'layer-interactions'}),
            'molecules': ET.SubElement(self.svg, 'g', {'id': 'layer-molecules'}),
            'labels': ET.SubElement(self.svg, 'g', {'id': 'layer-labels'}),
            'annotations': ET.SubElement(self.svg, 'g', {'id': 'layer-annotations'}),
            'legend': ET.SubElement(self.svg, 'g', {'id': 'layer-legend'})
        }

    def _add_defs(self):
        """Add SVG definitions (markers, gradients)."""
        defs = ET.SubElement(self.svg, 'defs')

        # Arrowhead marker (activation)
        marker = ET.SubElement(defs, 'marker', {
            'id': 'arrowhead',
            'markerWidth': '10',
            'markerHeight': '7',
            'refX': '9',
            'refY': '3.5',
            'orient': 'auto',
            'markerUnits': 'strokeWidth'
        })
        ET.SubElement(marker, 'polygon', {
            'points': '0 0, 10 3.5, 0 7',
            'fill': COLOR_PALETTES['interactions']['activation']
        })

        # T-bar marker (inhibition)
        tbar = ET.SubElement(defs, 'marker', {
            'id': 'tbar',
            'markerWidth': '10',
            'markerHeight': '10',
            'refX': '10',
            'refY': '5',
            'orient': 'auto'
        })
        ET.SubElement(tbar, 'line', {
            'x1': '0', 'y1': '0', 'x2': '0', 'y2': '10',
            'stroke': COLOR_PALETTES['interactions']['inhibition'],
            'stroke-width': '3'
        })

        # Phosphorylation circle
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
            'fill': COLOR_PALETTES['interactions']['phosphorylation']
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
            .molecule-label {
                font-family: Arial, Helvetica, sans-serif;
                font-size: 8pt;
                fill: #333;
            }
            .section-label {
                font-family: Arial, Helvetica, sans-serif;
                font-size: 10pt;
                font-weight: bold;
                fill: #333;
            }
            .annotation-note {
                font-family: Arial, Helvetica, sans-serif;
                font-size: 7pt;
                fill: #666;
                font-style: italic;
            }
            .annotation-ko {
                font-family: Arial, Helvetica, sans-serif;
                font-size: 7pt;
                fill: #E74C3C;
                font-weight: bold;
            }
            .annotation-inhibitor {
                font-family: Arial, Helvetica, sans-serif;
                font-size: 7pt;
                fill: #3498DB;
                font-weight: bold;
            }
            .legend-label {
                font-family: Arial, Helvetica, sans-serif;
                font-size: 7pt;
                fill: #333;
            }
        '''

    def add_compartment(
        self,
        name: str,
        x: float,
        y: float,
        width: float,
        height: float,
        color: str = '#f5f5f5',
        border_color: str = '#999',
        label_position: str = 'top-left'
    ):
        """
        Add a cellular compartment.

        Parameters
        ----------
        name : str
            Compartment name
        x, y : float
            Position
        width, height : float
            Dimensions
        color : str
            Fill color
        border_color : str
            Border color
        label_position : str
            Label position
        """
        comp = ET.SubElement(self.layers['compartments'], 'g', {
            'class': 'compartment',
            'id': f'comp-{name.replace(" ", "_")}'
        })

        # Background
        ET.SubElement(comp, 'rect', {
            'x': str(x),
            'y': str(y),
            'width': str(width),
            'height': str(height),
            'fill': color,
            'stroke': border_color,
            'stroke-width': '1.5',
            'rx': '5'
        })

        # Label
        label_x = x + 5
        label_y = y + 15
        if label_position == 'top-right':
            label_x = x + width - 5
        elif label_position == 'bottom-left':
            label_y = y + height - 5

        ET.SubElement(comp, 'text', {
            'x': str(label_x),
            'y': str(label_y),
            'class': 'section-label'
        }).text = name

        self.compartments.append({
            'name': name,
            'x': x,
            'y': y,
            'width': width,
            'height': height
        })

    def add_molecule(
        self,
        name: str,
        mol_type: str,
        x: float,
        y: float,
        color: Optional[str] = None,
        label: Optional[str] = None,
        size: float = 1.0
    ):
        """
        Add a molecule to the figure.

        Parameters
        ----------
        name : str
            Molecule name (used as ID)
        mol_type : str
            Type: receptor, kinase, tf, ligand, gene, complex
        x, y : float
            Position
        color : str, optional
            Override color
        label : str, optional
            Display label (if different from name)
        size : float
            Scale factor
        """
        if color is None:
            color = COLOR_PALETTES['signaling'].get(mol_type, '#666666')

        mol = Molecule(
            name=name,
            mol_type=mol_type,
            x=x,
            y=y,
            color=color,
            label=label or name,
            size=size
        )
        self.molecules[name] = mol

        # Create molecule group
        mol_group = ET.SubElement(self.layers['molecules'], 'g', {
            'class': f'molecule {mol_type}',
            'id': f'mol-{name}',
            'transform': f'translate({x}, {y}) scale({size})'
        })

        # Draw based on type
        if mol_type == 'receptor':
            self._draw_receptor(mol_group, color)
        elif mol_type == 'kinase':
            self._draw_kinase(mol_group, color)
        elif mol_type == 'tf':
            self._draw_tf(mol_group, color)
        elif mol_type == 'ligand':
            self._draw_ligand(mol_group, color)
        elif mol_type == 'gene':
            self._draw_gene(mol_group, color)
        elif mol_type == 'complex':
            self._draw_complex(mol_group, color)
        else:
            self._draw_generic(mol_group, color)

        # Add label
        ET.SubElement(mol_group, 'text', {
            'y': '25',
            'text-anchor': 'middle',
            'class': 'molecule-label'
        }).text = label or name

    def _draw_receptor(self, parent: ET.Element, color: str):
        """Draw receptor symbol."""
        # Transmembrane domain
        ET.SubElement(parent, 'rect', {
            'x': '-10', 'y': '-4', 'width': '20', 'height': '8',
            'fill': color, 'rx': '4'
        })
        # Extracellular domains
        ET.SubElement(parent, 'path', {
            'd': 'M-5,-4 L-5,-18 M5,-4 L5,-18',
            'stroke': color, 'stroke-width': '3', 'fill': 'none'
        })
        ET.SubElement(parent, 'circle', {
            'cx': '-5', 'cy': '-18', 'r': '5', 'fill': color
        })
        ET.SubElement(parent, 'circle', {
            'cx': '5', 'cy': '-18', 'r': '5', 'fill': color
        })
        # Intracellular domain
        ET.SubElement(parent, 'path', {
            'd': 'M0,4 L0,15',
            'stroke': color, 'stroke-width': '4', 'fill': 'none'
        })

    def _draw_kinase(self, parent: ET.Element, color: str):
        """Draw kinase symbol."""
        ET.SubElement(parent, 'rect', {
            'x': '-15', 'y': '-10', 'width': '30', 'height': '20',
            'fill': color, 'rx': '5'
        })
        ET.SubElement(parent, 'text', {
            'x': '0', 'y': '5',
            'text-anchor': 'middle',
            'fill': 'white',
            'font-size': '12',
            'font-weight': 'bold'
        }).text = 'P'

    def _draw_tf(self, parent: ET.Element, color: str):
        """Draw transcription factor symbol."""
        ET.SubElement(parent, 'polygon', {
            'points': '-12,10 -12,-5 0,-15 12,-5 12,10',
            'fill': color
        })

    def _draw_ligand(self, parent: ET.Element, color: str):
        """Draw ligand symbol."""
        ET.SubElement(parent, 'circle', {
            'r': '10', 'fill': color
        })

    def _draw_gene(self, parent: ET.Element, color: str):
        """Draw gene symbol."""
        # DNA line
        ET.SubElement(parent, 'line', {
            'x1': '-30', 'y1': '0', 'x2': '30', 'y2': '0',
            'stroke': '#333', 'stroke-width': '2'
        })
        # Arrow
        ET.SubElement(parent, 'polygon', {
            'points': '-5,-8 18,0 -5,8',
            'fill': color
        })

    def _draw_complex(self, parent: ET.Element, color: str):
        """Draw protein complex symbol."""
        ET.SubElement(parent, 'circle', {
            'cx': '-8', 'cy': '0', 'r': '14',
            'fill': color, 'opacity': '0.7'
        })
        ET.SubElement(parent, 'circle', {
            'cx': '8', 'cy': '0', 'r': '14',
            'fill': color, 'opacity': '0.9'
        })

    def _draw_generic(self, parent: ET.Element, color: str):
        """Draw generic molecule symbol."""
        ET.SubElement(parent, 'ellipse', {
            'rx': '15', 'ry': '10', 'fill': color
        })

    def add_interaction(
        self,
        source: str,
        target: str,
        interaction_type: str = 'activation',
        label: Optional[str] = None,
        curved: bool = False,
        control_offset: float = 30
    ):
        """
        Add an interaction between molecules.

        Parameters
        ----------
        source : str
            Source molecule name
        target : str
            Target molecule name
        interaction_type : str
            Type: activation, inhibition, phosphorylation, transcription, binding
        label : str, optional
            Interaction label
        curved : bool
            Use curved path
        control_offset : float
            Curve control point offset
        """
        if source not in self.molecules or target not in self.molecules:
            print(f"Warning: Molecule not found: {source} or {target}")
            return

        interaction = Interaction(
            source=source,
            target=target,
            interaction_type=interaction_type,
            label=label,
            curved=curved
        )
        self.interactions.append(interaction)

        # Get molecule positions
        src = self.molecules[source]
        tgt = self.molecules[target]

        # Create interaction group
        int_group = ET.SubElement(self.layers['interactions'], 'g', {
            'class': f'interaction {interaction_type}',
            'id': f'int-{source}-{target}'
        })

        # Determine color and marker
        color = COLOR_PALETTES['interactions'].get(interaction_type, '#666')

        marker_end = 'url(#arrowhead)'
        if interaction_type == 'inhibition':
            marker_end = 'url(#tbar)'
        elif interaction_type == 'phosphorylation':
            marker_end = 'url(#phospho)'

        # Draw line/path
        line_attrs = {
            'stroke': color,
            'stroke-width': '2',
            'fill': 'none',
            'marker-end': marker_end
        }

        if interaction_type == 'phosphorylation':
            line_attrs['stroke-dasharray'] = '5,3'

        if curved:
            # Bezier curve
            mx = (src.x + tgt.x) / 2
            my = (src.y + tgt.y) / 2 - control_offset
            path_d = f'M{src.x},{src.y + 15} Q{mx},{my} {tgt.x},{tgt.y - 20}'
            ET.SubElement(int_group, 'path', {**line_attrs, 'd': path_d})
        else:
            # Straight line
            ET.SubElement(int_group, 'line', {
                **line_attrs,
                'x1': str(src.x),
                'y1': str(src.y + 15),
                'x2': str(tgt.x),
                'y2': str(tgt.y - 20)
            })

        # Add label if provided
        if label:
            mx = (src.x + tgt.x) / 2
            my = (src.y + tgt.y) / 2
            ET.SubElement(int_group, 'text', {
                'x': str(mx),
                'y': str(my - 5),
                'text-anchor': 'middle',
                'class': 'annotation-note',
                'fill': color
            }).text = label

    def add_annotation(
        self,
        text: str,
        x: float,
        y: float,
        style: str = 'note',
        target: Optional[str] = None
    ):
        """
        Add annotation to figure.

        Parameters
        ----------
        text : str
            Annotation text
        x, y : float
            Position
        style : str
            Style: note, ko_effect, inhibitor, highlight
        target : str, optional
            Target molecule to annotate
        """
        annotation = Annotation(text=text, x=x, y=y, style=style)
        self.annotations.append(annotation)

        ann_group = ET.SubElement(self.layers['annotations'], 'g', {
            'class': f'annotation annotation-{style}',
            'transform': f'translate({x}, {y})'
        })

        css_class = f'annotation-{style.replace("_", "-")}'
        if style == 'note':
            css_class = 'annotation-note'
        elif style == 'ko_effect':
            css_class = 'annotation-ko'

        ET.SubElement(ann_group, 'text', {
            'class': css_class
        }).text = text

    def add_legend(
        self,
        x: float,
        y: float,
        include_interactions: bool = True,
        include_molecules: bool = True,
        include_cell_types: bool = False
    ):
        """
        Add figure legend.

        Parameters
        ----------
        x, y : float
            Legend position
        include_interactions : bool
            Include interaction types
        include_molecules : bool
            Include molecule types
        include_cell_types : bool
            Include cell type colors
        """
        legend = ET.SubElement(self.layers['legend'], 'g', {
            'transform': f'translate({x}, {y})'
        })

        # Title
        ET.SubElement(legend, 'text', {
            'class': 'section-label',
            'y': '0'
        }).text = 'Legend'

        y_offset = 25

        # Interaction types
        if include_interactions:
            interaction_items = [
                ('Activation', 'activation', '#2ECC71'),
                ('Inhibition', 'inhibition', '#E74C3C'),
                ('Phosphorylation', 'phosphorylation', '#3498DB')
            ]

            for label, itype, color in interaction_items:
                # Line sample
                ET.SubElement(legend, 'line', {
                    'x1': '0', 'y1': str(y_offset),
                    'x2': '25', 'y2': str(y_offset),
                    'stroke': color,
                    'stroke-width': '2',
                    'stroke-dasharray': '5,3' if itype == 'phosphorylation' else 'none'
                })

                # Add marker
                if itype == 'activation':
                    ET.SubElement(legend, 'polygon', {
                        'points': '22,-3 30,0 22,3',
                        'fill': color,
                        'transform': f'translate(0, {y_offset})'
                    })
                elif itype == 'inhibition':
                    ET.SubElement(legend, 'line', {
                        'x1': '25', 'y1': str(y_offset - 4),
                        'x2': '25', 'y2': str(y_offset + 4),
                        'stroke': color,
                        'stroke-width': '2'
                    })

                # Label
                ET.SubElement(legend, 'text', {
                    'x': '35',
                    'y': str(y_offset + 4),
                    'class': 'legend-label'
                }).text = label

                y_offset += 18

        # Molecule types
        if include_molecules:
            y_offset += 10
            ET.SubElement(legend, 'text', {
                'x': '0',
                'y': str(y_offset),
                'class': 'legend-label',
                'font-weight': 'bold'
            }).text = 'Molecules'
            y_offset += 15

            mol_items = [
                ('Receptor', 'receptor', '#2E86AB'),
                ('Kinase', 'kinase', '#A23B72'),
                ('TF', 'tf', '#F18F01'),
                ('Ligand', 'ligand', '#2ECC71')
            ]

            for label, mtype, color in mol_items:
                # Symbol
                if mtype == 'receptor':
                    ET.SubElement(legend, 'circle', {
                        'cx': '10', 'cy': str(y_offset - 3),
                        'r': '5', 'fill': color
                    })
                elif mtype == 'kinase':
                    ET.SubElement(legend, 'rect', {
                        'x': '3', 'y': str(y_offset - 8),
                        'width': '14', 'height': '10',
                        'fill': color, 'rx': '2'
                    })
                elif mtype == 'tf':
                    ET.SubElement(legend, 'polygon', {
                        'points': '2,3 10,-5 18,3',
                        'fill': color,
                        'transform': f'translate(0, {y_offset - 3})'
                    })
                else:
                    ET.SubElement(legend, 'circle', {
                        'cx': '10', 'cy': str(y_offset - 3),
                        'r': '6', 'fill': color
                    })

                # Label
                ET.SubElement(legend, 'text', {
                    'x': '25',
                    'y': str(y_offset),
                    'class': 'legend-label'
                }).text = label

                y_offset += 18

    def add_pathway(
        self,
        pathway_data: Any,
        layout: str = 'hierarchical',
        show_interactions: bool = True,
        show_regulation: bool = True,
        highlight_novel: Optional[List[str]] = None
    ):
        """
        Add pathway from analysis data.

        Parameters
        ----------
        pathway_data : PathwayResult
            Pathway analysis result
        layout : str
            Layout algorithm: hierarchical, circular, force
        show_interactions : bool
            Show molecular interactions
        show_regulation : bool
            Show regulatory relationships
        highlight_novel : list, optional
            Genes to highlight as novel
        """
        # Auto-layout molecules
        genes = pathway_data.genes if hasattr(pathway_data, 'genes') else []
        n_genes = len(genes)

        if layout == 'hierarchical':
            # Arrange in rows
            cols = 5
            for i, gene in enumerate(genes[:20]):  # Limit to 20
                row = i // cols
                col = i % cols
                x = 80 + col * 80
                y = 80 + row * 80

                color = '#E74C3C' if highlight_novel and gene in highlight_novel else None
                self.add_molecule(gene, 'kinase', x, y, color=color)

        # Add interactions
        if show_interactions and hasattr(pathway_data, 'interactions'):
            for interaction in pathway_data.interactions[:30]:  # Limit
                src = interaction.get('source', '')
                tgt = interaction.get('target', '')
                if src in self.molecules and tgt in self.molecules:
                    self.add_interaction(src, tgt, 'activation')

    def export(
        self,
        filename: str,
        format: Optional[str] = None
    ) -> str:
        """
        Export figure to file.

        Parameters
        ----------
        filename : str
            Output filename
        format : str, optional
            Format override (svg, pdf, png)

        Returns
        -------
        str
            Output filepath
        """
        if format is None:
            format = Path(filename).suffix.lstrip('.') or 'svg'

        filepath = Path(filename)
        if filepath.suffix == '':
            filepath = filepath.with_suffix(f'.{format}')

        if format == 'svg':
            tree = ET.ElementTree(self.svg)
            with open(filepath, 'wb') as f:
                f.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
                tree.write(f, encoding='unicode')

        elif format in ['pdf', 'png']:
            # First save SVG, then convert
            svg_path = filepath.with_suffix('.svg')
            self.export(str(svg_path), 'svg')

            # Use inkscape or rsvg-convert for conversion
            import subprocess
            try:
                if format == 'pdf':
                    subprocess.run([
                        'inkscape', str(svg_path),
                        '--export-filename', str(filepath),
                        '--export-type', 'pdf'
                    ], check=True, capture_output=True)
                elif format == 'png':
                    subprocess.run([
                        'inkscape', str(svg_path),
                        '--export-filename', str(filepath),
                        '--export-type', 'png',
                        '--export-dpi', str(self.resolution)
                    ], check=True, capture_output=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                print(f"Inkscape not available. SVG saved to {svg_path}")
                return str(svg_path)

        print(f"Figure exported to {filepath}")
        return str(filepath)

    def to_string(self) -> str:
        """Return SVG as string."""
        return ET.tostring(self.svg, encoding='unicode')


class PathwaySchematic(NatureReviewsFigure):
    """Specialized pathway schematic generator."""

    def create_signaling_cascade(
        self,
        cascade: List[Tuple[str, str]],
        title: str = 'Signaling Cascade',
        vertical: bool = True
    ):
        """
        Create a signaling cascade diagram.

        Parameters
        ----------
        cascade : list
            List of (molecule, type) tuples in order
        title : str
            Cascade title
        vertical : bool
            Vertical or horizontal layout
        """
        n = len(cascade)

        for i, (mol_name, mol_type) in enumerate(cascade):
            if vertical:
                x = self.px_width / 2
                y = 60 + i * 60
            else:
                x = 60 + i * 80
                y = self.px_height / 2

            self.add_molecule(mol_name, mol_type, x, y)

            # Add arrows between consecutive molecules
            if i > 0:
                prev_mol = cascade[i - 1][0]
                self.add_interaction(prev_mol, mol_name, 'phosphorylation')


if __name__ == '__main__':
    # Example usage
    fig = NatureReviewsFigure(
        width=180,
        height=120,
        title='VEGF Signaling Pathway'
    )

    # Add compartments
    fig.add_compartment('Extracellular', 10, 10, 160, 30, '#e8f4e8')
    fig.add_compartment('Cytoplasm', 10, 50, 160, 80, '#f5f5f5')

    # Add molecules
    fig.add_molecule('VEGF-A', 'ligand', 85, 25)
    fig.add_molecule('VEGFR2', 'receptor', 85, 60)
    fig.add_molecule('PLCγ', 'kinase', 50, 100)
    fig.add_molecule('ERK', 'kinase', 120, 100)

    # Add interactions
    fig.add_interaction('VEGF-A', 'VEGFR2', 'activation')
    fig.add_interaction('VEGFR2', 'PLCγ', 'phosphorylation')
    fig.add_interaction('VEGFR2', 'ERK', 'phosphorylation')

    # Add legend
    fig.add_legend(180, 20)

    # Export
    fig.export('vegf_pathway.svg')
