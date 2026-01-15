# Figure Generation Reference

## Nature Reviews Style Guide

### General Principles
1. **Clarity**: Every element should be immediately understandable
2. **Consistency**: Use consistent symbols, colors, and styles throughout
3. **Scientific accuracy**: All representations must be biologically accurate
4. **Accessibility**: Color-blind friendly palettes, clear labels

### Dimensions and Resolution
```python
figure_specs = {
    'single_column': {
        'width': 85,  # mm
        'max_height': 225  # mm
    },
    'double_column': {
        'width': 180,  # mm
        'max_height': 225  # mm
    },
    'resolution': {
        'print': 300,  # dpi minimum
        'high_quality': 600  # dpi for line art
    },
    'fonts': {
        'family': 'Arial, Helvetica, sans-serif',
        'size_min': 6,  # pt
        'size_labels': 8,  # pt
        'size_titles': 10  # pt
    }
}
```

### Color Palettes
```python
# Nature Reviews color schemes
color_palettes = {
    'signaling': {
        'receptor': '#2E86AB',      # Blue
        'kinase': '#A23B72',        # Magenta
        'transcription_factor': '#F18F01',  # Orange
        'effector': '#C73E1D',      # Red
        'adaptor': '#3B1F2B'        # Dark purple
    },

    'cell_types': {
        'tip_cell': '#E63946',      # Red
        'stalk_cell': '#457B9D',    # Blue
        'phalanx_cell': '#2A9D8F',  # Teal
        'pericyte': '#E9C46A',      # Yellow
        'smooth_muscle': '#264653'  # Dark blue
    },

    'interactions': {
        'activation': '#2ECC71',    # Green
        'inhibition': '#E74C3C',    # Red
        'phosphorylation': '#3498DB',  # Blue
        'transcription': '#9B59B6', # Purple
        'translocation': '#F39C12'  # Orange
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
```

## Pathway Schematic Elements

### Standard Molecular Symbols
```python
molecular_symbols = {
    'receptor': {
        'shape': 'Y-shape or mushroom',
        'description': 'Transmembrane receptor',
        'svg_template': '''
            <g class="receptor">
                <rect x="-10" y="-5" width="20" height="10" rx="5" fill="{color}"/>
                <path d="M-5,-5 L-5,-20 M5,-5 L5,-20" stroke="{color}" stroke-width="2"/>
                <circle cx="-5" cy="-20" r="5" fill="{color}"/>
                <circle cx="5" cy="-20" r="5" fill="{color}"/>
            </g>
        '''
    },

    'kinase': {
        'shape': 'rounded rectangle with P',
        'description': 'Protein kinase',
        'svg_template': '''
            <g class="kinase">
                <rect x="-15" y="-10" width="30" height="20" rx="5" fill="{color}"/>
                <text x="0" y="5" text-anchor="middle" fill="white" font-size="10">P</text>
            </g>
        '''
    },

    'transcription_factor': {
        'shape': 'pentagon or house shape',
        'description': 'DNA-binding protein',
        'svg_template': '''
            <g class="tf">
                <polygon points="-12,10 -12,-5 0,-15 12,-5 12,10" fill="{color}"/>
            </g>
        '''
    },

    'ligand': {
        'shape': 'circle or diamond',
        'description': 'Secreted signaling molecule',
        'svg_template': '''
            <g class="ligand">
                <circle r="10" fill="{color}"/>
            </g>
        '''
    },

    'gene': {
        'shape': 'arrow on DNA line',
        'description': 'Gene locus',
        'svg_template': '''
            <g class="gene">
                <line x1="-30" y1="0" x2="30" y2="0" stroke="#333" stroke-width="2"/>
                <polygon points="0,-8 20,0 0,8" fill="{color}"/>
            </g>
        '''
    },

    'complex': {
        'shape': 'interlocking shapes',
        'description': 'Protein complex',
        'svg_template': '''
            <g class="complex">
                <circle cx="-8" cy="0" r="12" fill="{color1}" opacity="0.8"/>
                <circle cx="8" cy="0" r="12" fill="{color2}" opacity="0.8"/>
            </g>
        '''
    }
}
```

### Interaction Arrows
```python
arrow_types = {
    'activation': {
        'style': 'solid line with arrowhead',
        'color': '#2ECC71',
        'svg': '''
            <defs>
                <marker id="arrowhead" markerWidth="10" markerHeight="7"
                    refX="9" refY="3.5" orient="auto">
                    <polygon points="0 0, 10 3.5, 0 7" fill="{color}"/>
                </marker>
            </defs>
            <line marker-end="url(#arrowhead)" stroke="{color}" stroke-width="2"/>
        '''
    },

    'inhibition': {
        'style': 'solid line with flat head (T-bar)',
        'color': '#E74C3C',
        'svg': '''
            <defs>
                <marker id="tbar" markerWidth="10" markerHeight="10"
                    refX="10" refY="5" orient="auto">
                    <line x1="0" y1="0" x2="0" y2="10" stroke="{color}" stroke-width="3"/>
                </marker>
            </defs>
            <line marker-end="url(#tbar)" stroke="{color}" stroke-width="2"/>
        '''
    },

    'phosphorylation': {
        'style': 'dashed line with P in circle',
        'color': '#3498DB',
        'svg': '''
            <line stroke="{color}" stroke-width="2" stroke-dasharray="5,3"/>
            <circle cx="{end_x}" cy="{end_y}" r="8" fill="{color}"/>
            <text x="{end_x}" y="{end_y}" dy="4" text-anchor="middle"
                fill="white" font-size="10">P</text>
        '''
    },

    'transcription': {
        'style': 'wavy line with arrowhead',
        'color': '#9B59B6',
        'description': 'Transcriptional regulation'
    },

    'translocation': {
        'style': 'dotted line with arrowhead',
        'color': '#F39C12',
        'description': 'Cellular translocation'
    },

    'binding': {
        'style': 'double-headed arrow',
        'color': '#34495E',
        'description': 'Physical binding/complex formation'
    }
}
```

## SVG Generation Code

### Base Pathway Figure Class
```python
import xml.etree.ElementTree as ET

class PathwayFigure:
    """Generate publication-quality pathway schematics."""

    def __init__(self, width=180, height=120, title=""):
        """
        Initialize pathway figure.

        Parameters:
        -----------
        width : float
            Figure width in mm
        height : float
            Figure height in mm
        title : str
            Figure title
        """
        self.width = width
        self.height = height
        self.title = title

        # Create SVG root
        self.svg = ET.Element('svg', {
            'xmlns': 'http://www.w3.org/2000/svg',
            'width': f'{width}mm',
            'height': f'{height}mm',
            'viewBox': f'0 0 {width * 3.78} {height * 3.78}'  # Convert to pixels
        })

        # Add style definitions
        self._add_styles()

        # Layer groups
        self.layers = {
            'background': ET.SubElement(self.svg, 'g', {'id': 'background'}),
            'compartments': ET.SubElement(self.svg, 'g', {'id': 'compartments'}),
            'molecules': ET.SubElement(self.svg, 'g', {'id': 'molecules'}),
            'interactions': ET.SubElement(self.svg, 'g', {'id': 'interactions'}),
            'labels': ET.SubElement(self.svg, 'g', {'id': 'labels'}),
            'annotations': ET.SubElement(self.svg, 'g', {'id': 'annotations'}),
            'legend': ET.SubElement(self.svg, 'g', {'id': 'legend'})
        }

    def _add_styles(self):
        """Add CSS styles."""
        style = ET.SubElement(self.svg, 'style')
        style.text = '''
            .molecule-label { font-family: Arial, sans-serif; font-size: 8pt; }
            .section-label { font-family: Arial, sans-serif; font-size: 10pt; font-weight: bold; }
            .annotation { font-family: Arial, sans-serif; font-size: 7pt; font-style: italic; }
            .receptor { stroke: #333; stroke-width: 1; }
            .kinase { stroke: #333; stroke-width: 1; }
        '''

    def add_compartment(self, name, x, y, width, height, color='#f0f0f0'):
        """Add cellular compartment (membrane, cytoplasm, nucleus, etc.)."""
        compartment = ET.SubElement(self.layers['compartments'], 'g', {
            'class': 'compartment',
            'id': f'compartment-{name}'
        })

        # Background
        ET.SubElement(compartment, 'rect', {
            'x': str(x),
            'y': str(y),
            'width': str(width),
            'height': str(height),
            'fill': color,
            'stroke': '#999',
            'stroke-width': '1',
            'rx': '5'
        })

        # Label
        ET.SubElement(compartment, 'text', {
            'x': str(x + 5),
            'y': str(y + 15),
            'class': 'section-label',
            'fill': '#666'
        }).text = name

        return compartment

    def add_molecule(self, name, mol_type, x, y, color=None):
        """
        Add molecule to figure.

        Parameters:
        -----------
        name : str
            Molecule name
        mol_type : str
            Type: 'receptor', 'kinase', 'tf', 'ligand', 'gene', 'complex'
        x, y : float
            Position
        color : str
            Hex color code
        """
        if color is None:
            color = color_palettes['signaling'].get(mol_type, '#666')

        mol_group = ET.SubElement(self.layers['molecules'], 'g', {
            'class': f'molecule {mol_type}',
            'id': f'mol-{name}',
            'transform': f'translate({x}, {y})'
        })

        # Add shape based on type
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

        # Add label
        ET.SubElement(mol_group, 'text', {
            'y': '25',
            'text-anchor': 'middle',
            'class': 'molecule-label'
        }).text = name

        return mol_group

    def _draw_receptor(self, parent, color):
        """Draw receptor symbol."""
        # Transmembrane domain
        ET.SubElement(parent, 'rect', {
            'x': '-8', 'y': '-3', 'width': '16', 'height': '6',
            'fill': color, 'rx': '3'
        })
        # Extracellular domains
        ET.SubElement(parent, 'path', {
            'd': 'M-4,-3 L-4,-15 M4,-3 L4,-15',
            'stroke': color, 'stroke-width': '2', 'fill': 'none'
        })
        ET.SubElement(parent, 'circle', {
            'cx': '-4', 'cy': '-15', 'r': '4', 'fill': color
        })
        ET.SubElement(parent, 'circle', {
            'cx': '4', 'cy': '-15', 'r': '4', 'fill': color
        })
        # Intracellular domain
        ET.SubElement(parent, 'path', {
            'd': 'M0,3 L0,12',
            'stroke': color, 'stroke-width': '3', 'fill': 'none'
        })

    def _draw_kinase(self, parent, color):
        """Draw kinase symbol."""
        ET.SubElement(parent, 'rect', {
            'x': '-12', 'y': '-8', 'width': '24', 'height': '16',
            'fill': color, 'rx': '4'
        })
        ET.SubElement(parent, 'text', {
            'x': '0', 'y': '4', 'text-anchor': 'middle',
            'fill': 'white', 'font-size': '10', 'font-weight': 'bold'
        }).text = 'P'

    def _draw_tf(self, parent, color):
        """Draw transcription factor symbol."""
        ET.SubElement(parent, 'polygon', {
            'points': '-10,8 -10,-4 0,-12 10,-4 10,8',
            'fill': color
        })

    def _draw_ligand(self, parent, color):
        """Draw ligand symbol."""
        ET.SubElement(parent, 'circle', {
            'r': '8', 'fill': color
        })

    def _draw_gene(self, parent, color):
        """Draw gene symbol."""
        # DNA line
        ET.SubElement(parent, 'line', {
            'x1': '-25', 'y1': '0', 'x2': '25', 'y2': '0',
            'stroke': '#333', 'stroke-width': '2'
        })
        # Arrow for transcription direction
        ET.SubElement(parent, 'polygon', {
            'points': '-5,-6 15,0 -5,6',
            'fill': color
        })

    def add_interaction(self, from_mol, to_mol, interaction_type, from_pos, to_pos):
        """
        Add interaction arrow between molecules.

        Parameters:
        -----------
        from_mol, to_mol : str
            Molecule names
        interaction_type : str
            'activation', 'inhibition', 'phosphorylation', etc.
        from_pos, to_pos : tuple
            (x, y) coordinates
        """
        color = arrow_types[interaction_type]['color']

        interaction_group = ET.SubElement(self.layers['interactions'], 'g', {
            'class': f'interaction {interaction_type}',
            'id': f'int-{from_mol}-{to_mol}'
        })

        # Draw line
        line_attrs = {
            'x1': str(from_pos[0]),
            'y1': str(from_pos[1]),
            'x2': str(to_pos[0]),
            'y2': str(to_pos[1]),
            'stroke': color,
            'stroke-width': '2'
        }

        if interaction_type == 'inhibition':
            line_attrs['marker-end'] = 'url(#tbar)'
        else:
            line_attrs['marker-end'] = 'url(#arrowhead)'

        if interaction_type == 'phosphorylation':
            line_attrs['stroke-dasharray'] = '5,3'

        ET.SubElement(interaction_group, 'line', line_attrs)

        return interaction_group

    def add_annotation(self, text, x, y, style='note'):
        """Add annotation text."""
        annotation = ET.SubElement(self.layers['annotations'], 'g', {
            'class': f'annotation {style}',
            'transform': f'translate({x}, {y})'
        })

        # Background
        text_elem = ET.SubElement(annotation, 'text', {
            'class': 'annotation',
            'fill': '#E74C3C' if style == 'ko_effect' else '#666'
        })
        text_elem.text = text

        return annotation

    def add_legend(self, x, y):
        """Add figure legend."""
        legend = ET.SubElement(self.layers['legend'], 'g', {
            'transform': f'translate({x}, {y})'
        })

        # Title
        ET.SubElement(legend, 'text', {
            'class': 'section-label',
            'y': '0'
        }).text = 'Legend'

        # Add legend items
        items = [
            ('Activation', 'activation', '#2ECC71'),
            ('Inhibition', 'inhibition', '#E74C3C'),
            ('Phosphorylation', 'phosphorylation', '#3498DB')
        ]

        for i, (label, itype, color) in enumerate(items):
            y_offset = 20 + i * 20

            # Line sample
            ET.SubElement(legend, 'line', {
                'x1': '0', 'y1': str(y_offset),
                'x2': '30', 'y2': str(y_offset),
                'stroke': color, 'stroke-width': '2'
            })

            # Label
            ET.SubElement(legend, 'text', {
                'x': '40', 'y': str(y_offset + 4),
                'class': 'molecule-label'
            }).text = label

        return legend

    def save(self, filename):
        """Save figure to file."""
        tree = ET.ElementTree(self.svg)

        # Add XML declaration and DOCTYPE
        with open(filename, 'wb') as f:
            f.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
            tree.write(f, encoding='unicode')

        print(f"Figure saved to {filename}")

    def to_string(self):
        """Return SVG as string."""
        return ET.tostring(self.svg, encoding='unicode')
```

## Complete Example: VEGF-Notch Pathway

```python
def create_vegf_notch_figure():
    """Create VEGF-Notch crosstalk figure."""

    fig = PathwayFigure(width=180, height=150, title="VEGF-Notch Crosstalk")

    # Add compartments
    fig.add_compartment("Tip Cell", 20, 20, 160, 60, '#FFE5E5')
    fig.add_compartment("Stalk Cell", 200, 20, 160, 60, '#E5E5FF')
    fig.add_compartment("Extracellular", 20, 0, 340, 20, '#E5FFE5')

    # Add molecules - Tip Cell
    fig.add_molecule("VEGF-A", "ligand", 100, 10, '#2ECC71')
    fig.add_molecule("VEGFR2", "receptor", 100, 40, '#3498DB')
    fig.add_molecule("DLL4", "receptor", 150, 40, '#E74C3C')
    fig.add_molecule("ERK", "kinase", 100, 80, '#9B59B6')

    # Add molecules - Stalk Cell
    fig.add_molecule("Notch1", "receptor", 220, 40, '#E74C3C')
    fig.add_molecule("HES1", "tf", 280, 80, '#F39C12')
    fig.add_molecule("Jagged1", "receptor", 280, 40, '#E74C3C')

    # Add interactions
    fig.add_interaction("VEGF-A", "VEGFR2", "activation", (100, 15), (100, 30))
    fig.add_interaction("VEGFR2", "ERK", "phosphorylation", (100, 55), (100, 70))
    fig.add_interaction("ERK", "DLL4", "activation", (110, 80), (145, 55))
    fig.add_interaction("DLL4", "Notch1", "activation", (160, 40), (210, 40))
    fig.add_interaction("Notch1", "HES1", "activation", (220, 55), (280, 70))
    fig.add_interaction("HES1", "VEGFR2", "inhibition", (290, 90), (110, 50))

    # Add annotations
    fig.add_annotation("DAPT blocks", 195, 35, style='inhibitor')
    fig.add_annotation("KO: No tip cells", 100, 100, style='ko_effect')

    # Add legend
    fig.add_legend(380, 20)

    # Save
    fig.save("vegf_notch_crosstalk.svg")

    return fig
```

## Export Options

### Multiple Formats
```python
import subprocess

def export_figure(svg_path, formats=['pdf', 'png', 'eps']):
    """Export SVG to multiple formats using Inkscape or rsvg-convert."""

    for fmt in formats:
        output_path = svg_path.replace('.svg', f'.{fmt}')

        if fmt == 'pdf':
            subprocess.run([
                'inkscape', svg_path,
                '--export-filename', output_path,
                '--export-type', 'pdf'
            ])
        elif fmt == 'png':
            subprocess.run([
                'inkscape', svg_path,
                '--export-filename', output_path,
                '--export-type', 'png',
                '--export-dpi', '600'
            ])
        elif fmt == 'eps':
            subprocess.run([
                'inkscape', svg_path,
                '--export-filename', output_path,
                '--export-type', 'eps'
            ])

        print(f"Exported to {output_path}")
```

### PowerPoint/Illustrator Compatible
```python
def make_editable(svg_content):
    """
    Modify SVG for better editability in vector programs.

    - Convert text to paths for font consistency
    - Ungroup complex elements
    - Add layer structure
    """
    # Use Inkscape for text-to-path conversion
    # inkscape input.svg --export-text-to-path --export-plain-svg output.svg
    pass
```
