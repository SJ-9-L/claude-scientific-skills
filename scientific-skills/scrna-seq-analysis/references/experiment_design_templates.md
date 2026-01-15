# Universal Experiment Design Templates

## Overview

어떤 유전자/pathway에도 적용 가능한 실험 설계 템플릿 시스템

## Knockout Comparison Experiment Template

### Template Structure

```python
class KOComparisonExperiment:
    """
    Universal knockout comparison experiment design.

    Supports any gene combination for systematic comparison.
    """

    def __init__(
        self,
        target_genes: list,
        species: str = 'mouse',
        context: str = None
    ):
        """
        Initialize KO comparison experiment.

        Parameters
        ----------
        target_genes : list
            Genes to knockout/compare (e.g., ['VEGFR2', 'VEGFR3', 'NOTCH1'])
        species : str
            'mouse', 'zebrafish', 'human_cell_line'
        context : str
            Biological context (e.g., 'angiogenesis', 'tumor', 'development')
        """
        self.target_genes = target_genes
        self.species = species
        self.context = context
        self.groups = self._generate_groups()

    def _generate_groups(self) -> dict:
        """Generate experimental groups automatically."""
        from itertools import combinations

        groups = {
            'wild_type': {
                'genotype': '+/+',
                'description': 'Wild-type control',
                'n_recommended': 10
            }
        }

        # Single KO groups
        for gene in self.target_genes:
            groups[f'{gene}_KO'] = {
                'genotype': f'{gene} -/-',
                'description': f'{gene} homozygous knockout',
                'n_recommended': 10
            }

            # Heterozygous if relevant
            groups[f'{gene}_het'] = {
                'genotype': f'{gene} +/-',
                'description': f'{gene} heterozygous',
                'n_recommended': 5
            }

        # Double KO groups (if 2+ genes)
        if len(self.target_genes) >= 2:
            for combo in combinations(self.target_genes, 2):
                key = f'{combo[0]}_{combo[1]}_DKO'
                groups[key] = {
                    'genotype': f'{combo[0]} -/-; {combo[1]} -/-',
                    'description': f'Double knockout of {combo[0]} and {combo[1]}',
                    'n_recommended': 10
                }

        return groups

    def generate_design(self) -> dict:
        """Generate complete experiment design."""
        return {
            'title': f'Knockout Comparison: {", ".join(self.target_genes)}',
            'species': self.species,
            'context': self.context,
            'groups': self.groups,
            'controls': self._get_controls(),
            'phenotypes_to_assess': self._get_phenotypes(),
            'timeline': self._get_timeline(),
            'statistical_analysis': self._get_statistics(),
            'expected_outcomes': self._predict_outcomes()
        }

    def _get_controls(self) -> dict:
        """Define control groups."""
        return {
            'negative_controls': ['wild_type', 'vehicle_treated'],
            'positive_controls': self._suggest_positive_controls(),
            'technical_controls': ['biological_replicates', 'batch_controls']
        }

    def _suggest_positive_controls(self) -> list:
        """Suggest positive controls based on context."""
        positive_controls = {
            'angiogenesis': ['anti-VEGF_treatment', 'hypoxia_induction'],
            'inflammation': ['LPS_treatment', 'TNFalpha_treatment'],
            'apoptosis': ['staurosporine_treatment', 'UV_irradiation'],
            'proliferation': ['serum_starvation', 'growth_factor_stimulation'],
            'migration': ['scratch_wound', 'chemotaxis_assay']
        }
        return positive_controls.get(self.context, ['standard_positive_control'])

    def _get_phenotypes(self) -> list:
        """Get phenotypes to assess based on context."""
        from .universal_phenotype_database import PATHWAY_PHENOTYPE_TEMPLATES

        if self.context in PATHWAY_PHENOTYPE_TEMPLATES:
            return PATHWAY_PHENOTYPE_TEMPLATES[self.context]['expected_phenotypes']

        # Default phenotypes
        return [
            'viability',
            'morphology',
            'cell_proliferation',
            'cell_migration',
            'gene_expression'
        ]

    def _get_timeline(self) -> dict:
        """Generate experiment timeline."""
        if self.species == 'mouse':
            return {
                'breeding': 'Week 1-8: Generate required genotypes',
                'validation': 'Week 9: Genotype confirmation',
                'baseline': 'Week 10: Baseline measurements',
                'intervention': 'Week 11-14: Experimental procedures',
                'endpoint': 'Week 15: Tissue collection',
                'analysis': 'Week 16-20: Data analysis'
            }
        elif self.species == 'zebrafish':
            return {
                'breeding': 'Day 1-5: Generate embryos',
                'injection': 'Day 1: Morpholino/CRISPR injection',
                'observation': 'Day 2-5: Live imaging',
                'endpoint': 'Day 5-7: Analysis'
            }
        else:
            return {
                'preparation': 'Day 1-3: Cell preparation',
                'transfection': 'Day 4: KO induction',
                'recovery': 'Day 5-7: Selection',
                'experiment': 'Day 8-14: Phenotype analysis'
            }

    def _get_statistics(self) -> dict:
        """Define statistical analysis plan."""
        n_groups = len(self.groups)

        return {
            'sample_size': {
                'per_group': 10,
                'total': n_groups * 10,
                'power': 0.8,
                'alpha': 0.05
            },
            'primary_analysis': 'One-way ANOVA' if n_groups <= 4 else 'Two-way ANOVA',
            'post_hoc': "Tukey's HSD for multiple comparisons",
            'effect_size': "Cohen's d for pairwise comparisons",
            'normality_test': 'Shapiro-Wilk test',
            'non_parametric_alternative': 'Kruskal-Wallis test'
        }

    def _predict_outcomes(self) -> dict:
        """Predict expected outcomes based on known data."""
        predictions = {}

        for gene in self.target_genes:
            predictions[f'{gene}_KO'] = {
                'predicted_severity': 'Query phenotype database',
                'expected_phenotypes': 'Based on pathway analysis',
                'compensatory_mechanisms': 'Check for related genes'
            }

        return predictions
```

### Usage Example

```python
# Design VEGFR2 vs VEGFR3 KO comparison
experiment = KOComparisonExperiment(
    target_genes=['KDR', 'FLT4'],  # VEGFR2, VEGFR3
    species='mouse',
    context='angiogenesis'
)

design = experiment.generate_design()
print(design)
```

## Treatment Comparison Experiment Template

```python
class TreatmentComparisonExperiment:
    """
    Universal treatment comparison experiment design.
    """

    def __init__(
        self,
        treatments: list,
        model_system: str,
        readouts: list = None
    ):
        """
        Initialize treatment comparison.

        Parameters
        ----------
        treatments : list
            Treatments to compare (e.g., ['DAPT', 'anti-VEGF', 'combination'])
        model_system : str
            Model system (e.g., 'mouse_tumor', 'cell_line', 'organoid')
        readouts : list, optional
            Phenotypes to measure
        """
        self.treatments = treatments
        self.model_system = model_system
        self.readouts = readouts or self._default_readouts()

    def _default_readouts(self) -> list:
        """Default readouts by model system."""
        readouts = {
            'mouse_tumor': [
                'tumor_volume',
                'tumor_weight',
                'vessel_density',
                'survival'
            ],
            'cell_line': [
                'viability',
                'proliferation',
                'migration',
                'invasion',
                'gene_expression'
            ],
            'organoid': [
                'organoid_size',
                'morphology',
                'differentiation_markers',
                'drug_response'
            ],
            'zebrafish': [
                'vessel_formation',
                'ISV_length',
                'DLAV_formation',
                'circulation'
            ]
        }
        return readouts.get(self.model_system, ['primary_endpoint'])

    def generate_design(self) -> dict:
        """Generate treatment comparison design."""
        groups = {'vehicle': {'treatment': 'Vehicle control', 'n': 10}}

        for tx in self.treatments:
            groups[tx] = {
                'treatment': tx,
                'n': 10,
                'dosing': self._get_dosing(tx),
                'schedule': self._get_schedule(tx)
            }

        # Combination groups
        if len(self.treatments) >= 2:
            from itertools import combinations
            for combo in combinations(self.treatments, 2):
                key = f'{combo[0]}_{combo[1]}_combo'
                groups[key] = {
                    'treatment': f'{combo[0]} + {combo[1]}',
                    'n': 10,
                    'rationale': 'Test synergy/additivity'
                }

        return {
            'groups': groups,
            'readouts': self.readouts,
            'timeline': self._get_timeline(),
            'analysis': self._get_analysis_plan()
        }

    def _get_dosing(self, treatment: str) -> dict:
        """Get dosing information for treatment."""
        # Query ChEMBL/literature for standard dosing
        dosing_db = {
            'DAPT': {'dose': '10 mg/kg', 'route': 'IP', 'frequency': 'daily'},
            'anti-VEGF': {'dose': '5 mg/kg', 'route': 'IP', 'frequency': '2x/week'},
            'Sunitinib': {'dose': '40 mg/kg', 'route': 'oral', 'frequency': 'daily'},
            'Sorafenib': {'dose': '30 mg/kg', 'route': 'oral', 'frequency': 'daily'}
        }
        return dosing_db.get(treatment, {'dose': 'Determine empirically'})

    def _get_schedule(self, treatment: str) -> str:
        """Get treatment schedule."""
        return 'Start at tumor volume ~100mm³, treat for 14-21 days'

    def _get_timeline(self) -> dict:
        """Generate timeline."""
        return {
            'day_-7': 'Tumor implantation',
            'day_0': 'Randomization when tumors reach ~100mm³',
            'day_1-21': 'Treatment period',
            'day_7,14,21': 'Interim measurements',
            'day_21': 'Endpoint analysis'
        }

    def _get_analysis_plan(self) -> dict:
        """Statistical analysis plan."""
        return {
            'tumor_growth': 'Mixed-effects model for repeated measures',
            'survival': 'Kaplan-Meier with log-rank test',
            'combination': 'Bliss independence model for synergy',
            'biomarkers': "Two-way ANOVA with Tukey's HSD"
        }
```

## Phenotype Assessment Templates

### Angiogenesis Phenotypes

```python
ANGIOGENESIS_PHENOTYPE_TEMPLATE = {
    'in_vivo_mouse': {
        'retinal_angiogenesis': {
            'timepoints': ['P3', 'P5', 'P7'],
            'measurements': [
                'vascular_density',
                'tip_cell_number',
                'filopodia_count',
                'radial_outgrowth',
                'vessel_branching'
            ],
            'staining': ['IB4', 'CollagenIV', 'ERG'],
            'analysis': 'AngioTool, ImageJ'
        },
        'tumor_angiogenesis': {
            'measurements': [
                'vessel_density_CD31',
                'vessel_maturity_aSMA',
                'hypoxia_pimonidazole',
                'perfusion_lectin'
            ],
            'imaging': ['IHC', 'IF', 'intravital']
        },
        'matrigel_plug': {
            'measurements': [
                'hemoglobin_content',
                'vessel_invasion',
                'inflammatory_infiltrate'
            ]
        }
    },

    'in_vitro': {
        'tube_formation': {
            'measurements': [
                'tube_length',
                'branch_points',
                'mesh_area',
                'tube_stability'
            ],
            'timepoints': ['4h', '8h', '16h', '24h'],
            'analysis': 'Angiogenesis Analyzer (ImageJ)'
        },
        'sprouting_assay': {
            'measurements': [
                'sprout_number',
                'sprout_length',
                'tip_cell_morphology'
            ]
        },
        'migration_assay': {
            'types': ['scratch_wound', 'transwell', 'chemotaxis'],
            'measurements': ['wound_closure', 'migrated_cells', 'migration_speed']
        }
    },

    'zebrafish': {
        'ISV_formation': {
            'measurements': [
                'ISV_length',
                'ISV_number',
                'DLAV_connection',
                'circulation'
            ],
            'timepoints': ['24hpf', '48hpf', '72hpf'],
            'imaging': 'Confocal live imaging'
        }
    }
}
```

### Generic Phenotype Templates

```python
PHENOTYPE_TEMPLATES = {
    'proliferation': {
        'assays': ['MTT/MTS', 'BrdU', 'Ki67', 'EdU', 'cell_counting'],
        'timepoints': ['24h', '48h', '72h', '96h'],
        'controls': ['serum_free', 'growth_factor_stimulation']
    },

    'apoptosis': {
        'assays': ['AnnexinV/PI', 'TUNEL', 'caspase3', 'SubG1'],
        'timepoints': ['6h', '12h', '24h', '48h'],
        'controls': ['staurosporine', 'serum_starvation']
    },

    'migration': {
        'assays': ['scratch_wound', 'transwell', 'chemotaxis', 'invasion'],
        'measurements': ['wound_closure%', 'migrated_cells', 'distance'],
        'controls': ['mitomycin_C_to_block_proliferation']
    },

    'differentiation': {
        'assays': ['marker_expression', 'functional_assay', 'morphology'],
        'readouts': ['qPCR', 'Western', 'IHC/IF', 'Flow_cytometry']
    },

    'signaling': {
        'assays': ['Western_blot', 'phospho_arrays', 'ELISA', 'reporter'],
        'readouts': ['phosphorylation', 'protein_level', 'activity'],
        'timepoints': ['0', '5min', '15min', '30min', '1h', '4h', '24h']
    }
}
```

## Experimental Design Generator

```python
class UniversalExperimentDesigner:
    """
    Generate experiment designs for any gene/pathway/treatment.
    """

    def __init__(self):
        self.phenotype_db = UniversalPhenotypeDB()

    def design_ko_experiment(
        self,
        genes: list,
        species: str = 'mouse',
        context: str = None,
        include_double_ko: bool = True,
        include_rescue: bool = False
    ) -> dict:
        """
        Design knockout experiment for any genes.

        Parameters
        ----------
        genes : list
            Genes to knockout
        species : str
            Species
        context : str
            Biological context
        include_double_ko : bool
            Include double knockouts
        include_rescue : bool
            Include rescue experiments

        Returns
        -------
        dict
            Complete experiment design
        """
        design = {
            'title': f'Knockout study: {", ".join(genes)}',
            'hypothesis': self._generate_hypothesis(genes, context),
            'groups': {},
            'controls': {},
            'phenotypes': {},
            'timeline': {},
            'statistics': {}
        }

        # Get expected phenotypes for each gene
        for gene in genes:
            pheno_data = self.phenotype_db.get_phenotypes(gene)
            design['phenotypes'][gene] = pheno_data

        # Generate groups
        exp = KOComparisonExperiment(genes, species, context)
        design['groups'] = exp.groups

        # Add rescue groups if requested
        if include_rescue:
            for gene in genes:
                design['groups'][f'{gene}_rescue'] = {
                    'genotype': f'{gene} -/- + {gene}_transgene',
                    'description': f'Rescue of {gene} knockout',
                    'purpose': 'Confirm phenotype specificity'
                }

        return design

    def design_treatment_experiment(
        self,
        target: str,
        treatments: list = None,
        model_system: str = 'cell_line'
    ) -> dict:
        """
        Design treatment experiment for any target.

        Parameters
        ----------
        target : str
            Target gene/pathway
        treatments : list, optional
            Treatments (auto-discovered if None)
        model_system : str
            Experimental model

        Returns
        -------
        dict
            Complete experiment design
        """
        # Auto-discover treatments if not provided
        if treatments is None:
            treatments = self._discover_treatments(target)

        exp = TreatmentComparisonExperiment(treatments, model_system)
        return exp.generate_design()

    def _generate_hypothesis(self, genes: list, context: str) -> str:
        """Generate hypothesis statement."""
        if len(genes) == 1:
            return f"Loss of {genes[0]} will result in {context}-related phenotypes"
        else:
            return (f"Comparing knockouts of {', '.join(genes)} will reveal "
                   f"shared and unique roles in {context}")

    def _discover_treatments(self, target: str) -> list:
        """Auto-discover treatments for a target."""
        from .molecular_tools_discovery import MolecularToolsDiscovery
        tools = MolecularToolsDiscovery()
        return tools.find_inhibitors(target)[:5]

    def generate_markdown_protocol(self, design: dict) -> str:
        """Generate markdown protocol from design."""
        md = f"""# Experimental Protocol

## {design['title']}

### Hypothesis
{design.get('hypothesis', 'To be defined')}

### Experimental Groups
| Group | Genotype/Treatment | N | Purpose |
|-------|-------------------|---|---------|
"""
        for name, info in design.get('groups', {}).items():
            genotype = info.get('genotype', info.get('treatment', '-'))
            n = info.get('n_recommended', info.get('n', 10))
            purpose = info.get('description', info.get('purpose', '-'))
            md += f"| {name} | {genotype} | {n} | {purpose} |\n"

        md += """
### Controls
- Negative control: Wild-type / Vehicle
- Positive control: Based on context

### Phenotype Assessment
"""
        for gene, pheno in design.get('phenotypes', {}).items():
            md += f"\n#### {gene}\n"
            if isinstance(pheno, dict) and 'summary' in pheno:
                for key, value in pheno['summary'].items():
                    md += f"- {key}: {value}\n"

        md += """
### Statistical Analysis
- Primary analysis: ANOVA with post-hoc tests
- Sample size: Calculated for 80% power at α=0.05
- Multiple comparison correction: Benjamini-Hochberg

### Timeline
See detailed timeline in supplementary materials.
"""
        return md
```

## Quick Design Commands

```python
# Design KO comparison for any genes
def quick_ko_design(genes: list, context: str = None) -> dict:
    """Quick knockout experiment design."""
    designer = UniversalExperimentDesigner()
    return designer.design_ko_experiment(genes, context=context)

# Design treatment experiment for any target
def quick_treatment_design(target: str) -> dict:
    """Quick treatment experiment design."""
    designer = UniversalExperimentDesigner()
    return designer.design_treatment_experiment(target)

# Example usage
design = quick_ko_design(['VEGFR2', 'VEGFR3', 'DLL4'], context='angiogenesis')
protocol = UniversalExperimentDesigner().generate_markdown_protocol(design)
print(protocol)
```
