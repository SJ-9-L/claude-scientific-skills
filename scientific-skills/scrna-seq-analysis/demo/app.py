"""
scRNA-seq Analysis & Pathway Discovery - Interactive Demo
Glassmorphism UI Style
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
import json

# Page configuration
st.set_page_config(
    page_title="scRNA-seq Pathway Discovery",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Glassmorphism CSS
st.markdown("""
<style>
    /* Import Google Fonts */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

    /* Global styles */
    * {
        font-family: 'Inter', sans-serif;
    }

    /* Gradient background */
    .stApp {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        background-attachment: fixed;
    }

    /* Animated gradient orbs */
    .stApp::before {
        content: '';
        position: fixed;
        top: -50%;
        left: -50%;
        width: 200%;
        height: 200%;
        background:
            radial-gradient(circle at 20% 80%, rgba(120, 0, 255, 0.15) 0%, transparent 50%),
            radial-gradient(circle at 80% 20%, rgba(0, 212, 255, 0.15) 0%, transparent 50%),
            radial-gradient(circle at 40% 40%, rgba(255, 0, 128, 0.1) 0%, transparent 40%);
        animation: float 20s ease-in-out infinite;
        pointer-events: none;
        z-index: -1;
    }

    @keyframes float {
        0%, 100% { transform: translate(0, 0) rotate(0deg); }
        33% { transform: translate(30px, -30px) rotate(120deg); }
        66% { transform: translate(-20px, 20px) rotate(240deg); }
    }

    /* Glassmorphism card */
    .glass-card {
        background: rgba(255, 255, 255, 0.05);
        backdrop-filter: blur(20px);
        -webkit-backdrop-filter: blur(20px);
        border-radius: 20px;
        border: 1px solid rgba(255, 255, 255, 0.1);
        padding: 2rem;
        margin: 1rem 0;
        box-shadow:
            0 8px 32px 0 rgba(0, 0, 0, 0.37),
            inset 0 1px 0 0 rgba(255, 255, 255, 0.1);
    }

    .glass-card-highlight {
        background: rgba(99, 102, 241, 0.1);
        backdrop-filter: blur(20px);
        -webkit-backdrop-filter: blur(20px);
        border-radius: 20px;
        border: 1px solid rgba(99, 102, 241, 0.3);
        padding: 2rem;
        margin: 1rem 0;
        box-shadow:
            0 8px 32px 0 rgba(99, 102, 241, 0.2),
            inset 0 1px 0 0 rgba(255, 255, 255, 0.1);
    }

    /* Sidebar styling */
    [data-testid="stSidebar"] {
        background: rgba(255, 255, 255, 0.03);
        backdrop-filter: blur(20px);
        -webkit-backdrop-filter: blur(20px);
        border-right: 1px solid rgba(255, 255, 255, 0.1);
    }

    [data-testid="stSidebar"] > div:first-child {
        background: transparent;
    }

    /* Headers */
    h1, h2, h3 {
        color: #ffffff !important;
        font-weight: 600;
    }

    h1 {
        background: linear-gradient(135deg, #00d4ff, #7c3aed, #f472b6);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        font-size: 2.5rem !important;
    }

    /* Text */
    p, span, label, .stMarkdown {
        color: rgba(255, 255, 255, 0.9) !important;
    }

    /* Input fields */
    .stTextInput > div > div > input,
    .stSelectbox > div > div > div,
    .stMultiSelect > div > div > div {
        background: rgba(255, 255, 255, 0.05) !important;
        border: 1px solid rgba(255, 255, 255, 0.1) !important;
        border-radius: 12px !important;
        color: white !important;
        backdrop-filter: blur(10px);
    }

    .stTextInput > div > div > input:focus {
        border-color: rgba(99, 102, 241, 0.5) !important;
        box-shadow: 0 0 20px rgba(99, 102, 241, 0.3) !important;
    }

    /* Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #6366f1, #8b5cf6) !important;
        color: white !important;
        border: none !important;
        border-radius: 12px !important;
        padding: 0.75rem 2rem !important;
        font-weight: 600 !important;
        transition: all 0.3s ease !important;
        box-shadow: 0 4px 15px rgba(99, 102, 241, 0.4) !important;
    }

    .stButton > button:hover {
        transform: translateY(-2px) !important;
        box-shadow: 0 6px 20px rgba(99, 102, 241, 0.6) !important;
    }

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        background: rgba(255, 255, 255, 0.05);
        border-radius: 15px;
        padding: 0.5rem;
        gap: 0.5rem;
    }

    .stTabs [data-baseweb="tab"] {
        background: transparent;
        color: rgba(255, 255, 255, 0.7);
        border-radius: 10px;
        padding: 0.5rem 1rem;
    }

    .stTabs [aria-selected="true"] {
        background: rgba(99, 102, 241, 0.3) !important;
        color: white !important;
    }

    /* Metrics */
    [data-testid="stMetricValue"] {
        color: #00d4ff !important;
        font-weight: 700;
    }

    [data-testid="stMetricLabel"] {
        color: rgba(255, 255, 255, 0.7) !important;
    }

    /* DataFrames */
    .dataframe {
        background: rgba(255, 255, 255, 0.05) !important;
        border-radius: 12px !important;
        overflow: hidden;
    }

    .dataframe th {
        background: rgba(99, 102, 241, 0.2) !important;
        color: white !important;
    }

    .dataframe td {
        color: rgba(255, 255, 255, 0.9) !important;
        border-color: rgba(255, 255, 255, 0.1) !important;
    }

    /* Expander */
    .streamlit-expanderHeader {
        background: rgba(255, 255, 255, 0.05) !important;
        border-radius: 12px !important;
        color: white !important;
    }

    /* Success/Info/Warning boxes */
    .stAlert {
        background: rgba(255, 255, 255, 0.05) !important;
        backdrop-filter: blur(10px);
        border-radius: 12px !important;
        border: 1px solid rgba(255, 255, 255, 0.1) !important;
    }

    /* Code blocks */
    .stCodeBlock {
        background: rgba(0, 0, 0, 0.3) !important;
        border-radius: 12px !important;
    }

    /* Plotly charts background */
    .js-plotly-plot .plotly .bg {
        fill: transparent !important;
    }

    /* Hide Streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}

    /* Custom scrollbar */
    ::-webkit-scrollbar {
        width: 8px;
        height: 8px;
    }

    ::-webkit-scrollbar-track {
        background: rgba(255, 255, 255, 0.05);
    }

    ::-webkit-scrollbar-thumb {
        background: rgba(99, 102, 241, 0.5);
        border-radius: 4px;
    }

    ::-webkit-scrollbar-thumb:hover {
        background: rgba(99, 102, 241, 0.7);
    }

    /* Feature card */
    .feature-card {
        background: rgba(255, 255, 255, 0.03);
        backdrop-filter: blur(10px);
        border-radius: 16px;
        border: 1px solid rgba(255, 255, 255, 0.08);
        padding: 1.5rem;
        margin: 0.5rem 0;
        transition: all 0.3s ease;
    }

    .feature-card:hover {
        background: rgba(255, 255, 255, 0.06);
        border-color: rgba(99, 102, 241, 0.3);
        transform: translateY(-2px);
    }

    /* Stat pill */
    .stat-pill {
        display: inline-block;
        background: linear-gradient(135deg, rgba(99, 102, 241, 0.2), rgba(139, 92, 246, 0.2));
        border: 1px solid rgba(99, 102, 241, 0.3);
        border-radius: 20px;
        padding: 0.3rem 1rem;
        margin: 0.2rem;
        font-size: 0.85rem;
        color: #a5b4fc;
    }

    /* Glow effect */
    .glow {
        text-shadow: 0 0 20px rgba(99, 102, 241, 0.5);
    }
</style>
""", unsafe_allow_html=True)


# ============== DATA FUNCTIONS ==============

def get_sample_genes():
    """Sample gene data for demonstration"""
    return {
        "VEGF Pathway": ["VEGFA", "VEGFR2", "VEGFR3", "NRP1", "DLL4", "NOTCH1", "HES1", "HEY1"],
        "EGFR Pathway": ["EGFR", "ERBB2", "ERBB3", "GRB2", "SOS1", "KRAS", "BRAF", "MEK1", "ERK1"],
        "Notch Pathway": ["NOTCH1", "NOTCH2", "DLL1", "DLL4", "JAG1", "JAG2", "HES1", "HEY1"],
        "Wnt Pathway": ["WNT3A", "FZD1", "LRP5", "CTNNB1", "APC", "AXIN1", "TCF7", "LEF1"],
        "TGF-β Pathway": ["TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7"],
    }

def get_phenotype_data(gene: str):
    """Get mock phenotype data"""
    phenotypes = {
        "VEGFR2": {
            "mouse_ko": [
                {"phenotype": "Embryonic lethality (E8.5-9.5)", "severity": "Lethal", "evidence": "PMID:9651508"},
                {"phenotype": "Absent vascular development", "severity": "Severe", "evidence": "PMID:9651508"},
                {"phenotype": "No endothelial cell differentiation", "severity": "Severe", "evidence": "PMID:9651507"},
            ],
            "cell_effects": [
                {"cell_type": "Tip cells", "effect": "Abolished migration", "fold_change": -100},
                {"cell_type": "Stalk cells", "effect": "Reduced proliferation", "fold_change": -75},
                {"cell_type": "Pericytes", "effect": "Failed recruitment", "fold_change": -60},
            ]
        },
        "DLL4": {
            "mouse_ko": [
                {"phenotype": "Embryonic lethality (E10.5)", "severity": "Lethal", "evidence": "PMID:15652354"},
                {"phenotype": "Hypervascularization", "severity": "Severe", "evidence": "PMID:15652354"},
                {"phenotype": "Arteriovenous malformation", "severity": "Severe", "evidence": "PMID:15652354"},
            ],
            "cell_effects": [
                {"cell_type": "Tip cells", "effect": "Excess formation", "fold_change": 300},
                {"cell_type": "Stalk cells", "effect": "Reduced specification", "fold_change": -50},
            ]
        }
    }
    return phenotypes.get(gene, {
        "mouse_ko": [{"phenotype": "Data not available", "severity": "-", "evidence": "-"}],
        "cell_effects": [{"cell_type": "-", "effect": "-", "fold_change": 0}]
    })

def get_molecular_tools(gene: str):
    """Get molecular tools for a gene"""
    tools = {
        "antibodies": [
            {"name": f"Anti-{gene} (Clone A1)", "application": "WB, IHC, IF", "species": "Human, Mouse", "vendor": "Cell Signaling"},
            {"name": f"Anti-{gene} (Clone B2)", "application": "Flow, IP", "species": "Human", "vendor": "BioLegend"},
            {"name": f"Anti-{gene} blocking Ab", "application": "Functional", "species": "Human, Mouse", "vendor": "R&D Systems"},
        ],
        "inhibitors": [
            {"name": "Compound-A", "IC50": "12 nM", "selectivity": "High", "stage": "Approved"},
            {"name": "Compound-B", "IC50": "45 nM", "selectivity": "Medium", "stage": "Phase 3"},
            {"name": "Compound-C", "IC50": "8 nM", "selectivity": "High", "stage": "Phase 2"},
        ],
        "siRNA": [
            {"sequence": "GCUUAAGCUAGGCUAAUUA", "efficiency": "92%", "vendor": "Dharmacon"},
            {"sequence": "UUACGGAUUCCAAGGCCUA", "efficiency": "87%", "vendor": "Ambion"},
        ],
        "crispr": [
            {"guide": f"sg{gene}-1", "efficiency": "89%", "off_target": "Low"},
            {"guide": f"sg{gene}-2", "efficiency": "94%", "off_target": "Very Low"},
        ]
    }
    return tools

def generate_umap_data(n_cells=500):
    """Generate mock UMAP data"""
    np.random.seed(42)

    clusters = {
        "Tip cells": (2, 3, 80),
        "Stalk cells": (-2, 2, 120),
        "Phalanx cells": (0, -2, 100),
        "Arterial EC": (-3, -1, 90),
        "Venous EC": (3, -2, 110),
    }

    data = []
    for cluster, (cx, cy, n) in clusters.items():
        x = np.random.normal(cx, 0.8, n)
        y = np.random.normal(cy, 0.8, n)
        for i in range(n):
            data.append({
                "UMAP1": x[i],
                "UMAP2": y[i],
                "Cluster": cluster,
                "Gene_Expression": np.random.exponential(2) if cluster == "Tip cells" else np.random.exponential(0.5)
            })

    return pd.DataFrame(data)


# ============== MAIN APP ==============

def main():
    # Sidebar
    with st.sidebar:
        st.markdown("""
        <div style="text-align: center; padding: 1rem 0;">
            <h2 style="background: linear-gradient(135deg, #00d4ff, #7c3aed); -webkit-background-clip: text; -webkit-text-fill-color: transparent; font-size: 1.5rem;">🧬 scRNA-seq</h2>
            <p style="color: rgba(255,255,255,0.6); font-size: 0.85rem;">Pathway Discovery Tool</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("---")

        page = st.radio(
            "Navigation",
            ["🏠 Home", "🔬 Novel Gene Discovery", "🧪 Experiment Design", "🗺️ Pathway Analysis", "📊 Competitive Analysis", "🎨 Figure Generator"],
            label_visibility="collapsed"
        )

        st.markdown("---")

        st.markdown("""
        <div class="feature-card">
            <h4 style="color: #a5b4fc; margin-bottom: 0.5rem;">📚 Databases</h4>
            <span class="stat-pill">28 Bioinformatics</span>
            <span class="stat-pill">12 Clinical</span>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div style="margin-top: 2rem; text-align: center; color: rgba(255,255,255,0.4); font-size: 0.75rem;">
            Built with Claude Scientific Skills<br>
            v2.17.0
        </div>
        """, unsafe_allow_html=True)

    # Main content
    if "🏠 Home" in page:
        render_home()
    elif "🔬 Novel Gene" in page:
        render_novel_discovery()
    elif "🧪 Experiment" in page:
        render_experiment_design()
    elif "🗺️ Pathway" in page:
        render_pathway_analysis()
    elif "📊 Competitive" in page:
        render_competitive()
    elif "🎨 Figure" in page:
        render_figure_generator()


def render_home():
    """Render home page"""
    st.markdown("""
    <h1 style="text-align: center; margin-bottom: 0;">scRNA-seq Analysis</h1>
    <p style="text-align: center; color: rgba(255,255,255,0.6); font-size: 1.1rem; margin-bottom: 2rem;">
        & Pathway Discovery Platform
    </p>
    """, unsafe_allow_html=True)

    # Hero metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown("""
        <div class="glass-card" style="text-align: center;">
            <h1 style="font-size: 2.5rem; margin: 0; color: #00d4ff;">40+</h1>
            <p style="margin: 0; color: rgba(255,255,255,0.7);">Databases</p>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="glass-card" style="text-align: center;">
            <h1 style="font-size: 2.5rem; margin: 0; color: #7c3aed;">5</h1>
            <p style="margin: 0; color: rgba(255,255,255,0.7);">Analysis Pipelines</p>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        st.markdown("""
        <div class="glass-card" style="text-align: center;">
            <h1 style="font-size: 2.5rem; margin: 0; color: #f472b6;">∞</h1>
            <p style="margin: 0; color: rgba(255,255,255,0.7);">Pathways Supported</p>
        </div>
        """, unsafe_allow_html=True)
    with col4:
        st.markdown("""
        <div class="glass-card" style="text-align: center;">
            <h1 style="font-size: 2.5rem; margin: 0; color: #34d399;">100%</h1>
            <p style="margin: 0; color: rgba(255,255,255,0.7);">Universal</p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # Features
    st.markdown("""
    <div class="glass-card">
        <h2 style="margin-bottom: 1.5rem;">✨ Core Features</h2>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div class="feature-card">
            <h3>🔬 Novel Gene Discovery</h3>
            <p>scRNA-seq 데이터에서 새로운 마커 유전자를 발굴하고 40+ 데이터베이스와 비교하여 novelty 검증</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="feature-card">
            <h3>🧪 Experiment Design</h3>
            <p>KO 비교 실험 자동 설계: VEGFR2 KO vs VEGFR3 KO vs WT, phenotype 예측, 필요 reagents 리스트</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="feature-card">
            <h3>🗺️ Deep Pathway Analysis</h3>
            <p>단순 up/downregulation이 아닌 실험적 증거 기반의 상세 pathway 분석</p>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class="feature-card">
            <h3>📊 Competitive Intelligence</h3>
            <p>경쟁 연구 그룹, 최신 동향, 특허 현황, research gap 분석</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="feature-card">
            <h3>🎨 Figure Generation</h3>
            <p>Nature Reviews 수준의 pathway schematic 자동 생성 (SVG/PDF)</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="feature-card">
            <h3>💊 Molecular Tools</h3>
            <p>Antibodies, inhibitors, siRNAs, CRISPR guides 자동 검색</p>
        </div>
        """, unsafe_allow_html=True)


def render_novel_discovery():
    """Render novel gene discovery page"""
    st.markdown("<h1>🔬 Novel Gene Discovery</h1>", unsafe_allow_html=True)

    st.markdown("""
    <div class="glass-card">
        <p>scRNA-seq 데이터에서 발굴한 후보 유전자의 novelty를 검증합니다.</p>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.markdown("### 🎯 Input Parameters")

        gene = st.text_input("Gene Symbol", value="GENE_X", placeholder="Enter gene symbol")

        cell_type = st.selectbox(
            "Cell Type Context",
            ["Tip cells", "Stalk cells", "Endothelial cells", "Macrophages", "T cells"]
        )

        databases = st.multiselect(
            "Reference Databases",
            ["PubMed", "Human Cell Atlas", "CellxGene", "GTEx", "TCGA"],
            default=["PubMed", "Human Cell Atlas", "CellxGene"]
        )

        if st.button("🔍 Analyze Novelty", use_container_width=True):
            st.session_state.analyzed = True

    with col2:
        if st.session_state.get('analyzed'):
            st.markdown("### 📊 Novelty Analysis Results")

            # Novelty score visualization
            novelty_score = 0.92

            st.markdown(f"""
            <div class="glass-card-highlight" style="text-align: center;">
                <h2 style="font-size: 3rem; margin: 0; color: #34d399;">{novelty_score:.0%}</h2>
                <p style="color: rgba(255,255,255,0.7); margin: 0;">Novelty Score</p>
                <p style="color: #34d399; margin-top: 0.5rem;">✓ TRULY NOVEL</p>
            </div>
            """, unsafe_allow_html=True)

            # Database search results
            results = pd.DataFrame({
                "Database": databases + ["STRING", "Reactome"],
                "Hits": [0, 0, 0, 2, 1] if len(databases) >= 3 else [0, 0],
                "Relevant": ["No", "No", "No", "No", "No"][:len(databases)+2],
                "Status": ["✓ Novel", "✓ Novel", "✓ Novel", "✓ Not related", "✓ Not related"][:len(databases)+2]
            })

            st.dataframe(results, use_container_width=True, hide_index=True)

            st.markdown("### 📈 Expression Profile")

            # Mock expression data
            umap_data = generate_umap_data()

            fig = px.scatter(
                umap_data, x="UMAP1", y="UMAP2",
                color="Gene_Expression",
                color_continuous_scale="Viridis",
                title=f"{gene} Expression in {cell_type}",
                template="plotly_dark"
            )
            fig.update_layout(
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
                font_color='white'
            )
            st.plotly_chart(fig, use_container_width=True)


def render_experiment_design():
    """Render experiment design page"""
    st.markdown("<h1>🧪 Experiment Design</h1>", unsafe_allow_html=True)

    st.markdown("""
    <div class="glass-card">
        <p>KO 비교 실험을 자동으로 설계하고 필요한 모든 정보를 제공합니다.</p>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.markdown("### 🎯 Experiment Setup")

        genes = st.multiselect(
            "Target Genes (KO)",
            ["VEGFR2", "VEGFR3", "DLL4", "NOTCH1", "EGFR", "BRAF", "KRAS"],
            default=["VEGFR2", "VEGFR3"]
        )

        species = st.selectbox("Species", ["Mouse", "Human", "Zebrafish"])

        context = st.selectbox(
            "Research Context",
            ["Angiogenesis", "Tumor", "Development", "Inflammation"]
        )

        treatment = st.text_input("Treatment (optional)", placeholder="e.g., DAPT")

        if st.button("🧬 Generate Design", use_container_width=True):
            st.session_state.designed = True

    with col2:
        if st.session_state.get('designed') and genes:
            st.markdown("### 📋 Experiment Design")

            # Groups
            st.markdown("""
            <div class="glass-card">
                <h4 style="color: #a5b4fc;">Experimental Groups</h4>
            </div>
            """, unsafe_allow_html=True)

            groups = ["Wild-type (WT)"] + [f"{g} KO" for g in genes]
            if treatment:
                groups = groups + [f"WT + {treatment}"] + [f"{g} KO + {treatment}" for g in genes]

            groups_df = pd.DataFrame({
                "Group": groups,
                "N": [8] * len(groups),
                "Controls": ["—"] + ["WT"] * (len(groups)-1)
            })
            st.dataframe(groups_df, use_container_width=True, hide_index=True)

            # Phenotypes
            st.markdown("""
            <div class="glass-card">
                <h4 style="color: #a5b4fc;">Phenotypes to Measure</h4>
            </div>
            """, unsafe_allow_html=True)

            if context == "Angiogenesis":
                phenotypes = [
                    "Tip cell formation (retina P5)",
                    "Vessel density (CD31 staining)",
                    "Sprouting angiogenesis (aortic ring)",
                    "Filopodia count per tip cell",
                    "Vascular branching index"
                ]
            else:
                phenotypes = ["Phenotype 1", "Phenotype 2", "Phenotype 3"]

            for p in phenotypes:
                st.markdown(f"- {p}")

            # Predicted outcomes
            if genes:
                st.markdown("""
                <div class="glass-card-highlight">
                    <h4 style="color: #f472b6;">📊 Predicted Phenotypes (from database)</h4>
                </div>
                """, unsafe_allow_html=True)

                for gene in genes:
                    pheno_data = get_phenotype_data(gene)
                    with st.expander(f"🧬 {gene} KO Predictions", expanded=True):
                        st.dataframe(
                            pd.DataFrame(pheno_data["mouse_ko"]),
                            use_container_width=True,
                            hide_index=True
                        )


def render_pathway_analysis():
    """Render pathway analysis page"""
    st.markdown("<h1>🗺️ Deep Pathway Analysis</h1>", unsafe_allow_html=True)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.markdown("### 🎯 Select Pathway")

        pathway = st.selectbox(
            "Pathway",
            list(get_sample_genes().keys())
        )

        gene = st.selectbox(
            "Focus Gene",
            get_sample_genes()[pathway]
        )

        analysis_type = st.multiselect(
            "Analysis Types",
            ["KO Studies", "Treatments", "Cell-type specificity", "Crosstalk"],
            default=["KO Studies", "Cell-type specificity"]
        )

        if st.button("🔬 Analyze Pathway", use_container_width=True):
            st.session_state.pathway_analyzed = True

    with col2:
        if st.session_state.get('pathway_analyzed'):
            tabs = st.tabs(["📊 Overview", "🧬 KO Effects", "💊 Molecular Tools", "🔗 Interactions"])

            with tabs[0]:
                st.markdown(f"""
                <div class="glass-card">
                    <h3>{gene} in {pathway}</h3>
                    <p>Comprehensive pathway analysis with experimental evidence</p>
                </div>
                """, unsafe_allow_html=True)

                # Network visualization
                fig = go.Figure()

                genes = get_sample_genes()[pathway]
                n = len(genes)
                angles = np.linspace(0, 2*np.pi, n, endpoint=False)
                x = np.cos(angles) * 2
                y = np.sin(angles) * 2

                # Add edges
                for i in range(n-1):
                    fig.add_trace(go.Scatter(
                        x=[x[i], x[i+1]], y=[y[i], y[i+1]],
                        mode='lines',
                        line=dict(color='rgba(99, 102, 241, 0.3)', width=2),
                        hoverinfo='none'
                    ))

                # Add nodes
                colors = ['#f472b6' if g == gene else '#6366f1' for g in genes]
                sizes = [30 if g == gene else 20 for g in genes]

                fig.add_trace(go.Scatter(
                    x=x, y=y,
                    mode='markers+text',
                    marker=dict(size=sizes, color=colors),
                    text=genes,
                    textposition='top center',
                    textfont=dict(color='white'),
                    hoverinfo='text'
                ))

                fig.update_layout(
                    showlegend=False,
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)',
                    xaxis=dict(visible=False),
                    yaxis=dict(visible=False),
                    height=400
                )

                st.plotly_chart(fig, use_container_width=True)

            with tabs[1]:
                pheno_data = get_phenotype_data(gene)

                st.markdown("#### Mouse KO Phenotypes")
                st.dataframe(
                    pd.DataFrame(pheno_data["mouse_ko"]),
                    use_container_width=True,
                    hide_index=True
                )

                st.markdown("#### Cell-Type Specific Effects")

                effects_df = pd.DataFrame(pheno_data["cell_effects"])

                fig = px.bar(
                    effects_df,
                    x="cell_type",
                    y="fold_change",
                    color="fold_change",
                    color_continuous_scale=["#ef4444", "#f59e0b", "#22c55e"],
                    title=f"{gene} KO Effects by Cell Type",
                    template="plotly_dark"
                )
                fig.update_layout(
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)'
                )
                st.plotly_chart(fig, use_container_width=True)

            with tabs[2]:
                tools = get_molecular_tools(gene)

                st.markdown("#### 🧫 Antibodies")
                st.dataframe(pd.DataFrame(tools["antibodies"]), use_container_width=True, hide_index=True)

                st.markdown("#### 💊 Small Molecule Inhibitors")
                st.dataframe(pd.DataFrame(tools["inhibitors"]), use_container_width=True, hide_index=True)

                col_a, col_b = st.columns(2)
                with col_a:
                    st.markdown("#### 🔬 siRNA")
                    st.dataframe(pd.DataFrame(tools["siRNA"]), use_container_width=True, hide_index=True)
                with col_b:
                    st.markdown("#### ✂️ CRISPR")
                    st.dataframe(pd.DataFrame(tools["crispr"]), use_container_width=True, hide_index=True)

            with tabs[3]:
                st.markdown("#### Pathway Interactions")

                interactions = pd.DataFrame({
                    "Source": [gene, gene, genes[1] if len(genes) > 1 else gene],
                    "Target": [genes[1] if len(genes) > 1 else gene, genes[2] if len(genes) > 2 else gene, genes[3] if len(genes) > 3 else gene],
                    "Type": ["Activation", "Inhibition", "Binding"],
                    "Evidence": ["Experimental", "Database", "Predicted"],
                    "Confidence": [0.95, 0.87, 0.72]
                })

                st.dataframe(interactions, use_container_width=True, hide_index=True)


def render_competitive():
    """Render competitive analysis page"""
    st.markdown("<h1>📊 Competitive Analysis</h1>", unsafe_allow_html=True)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.markdown("### 🎯 Search Parameters")

        topic = st.text_input("Research Topic", value="VEGF signaling tumor angiogenesis")

        timeframe = st.slider("Timeframe", 2015, 2024, (2020, 2024))

        sources = st.multiselect(
            "Sources",
            ["PubMed", "bioRxiv", "Patents", "Clinical Trials"],
            default=["PubMed", "bioRxiv"]
        )

        if st.button("🔍 Analyze Landscape", use_container_width=True):
            st.session_state.competitive_done = True

    with col2:
        if st.session_state.get('competitive_done'):
            tabs = st.tabs(["📈 Trends", "👥 Key Groups", "📑 Patents", "🎯 Gaps"])

            with tabs[0]:
                # Publication trends
                years = list(range(timeframe[0], timeframe[1] + 1))
                pubs = [45, 62, 78, 95, 120] if len(years) >= 5 else [45 + i*15 for i in range(len(years))]

                fig = px.area(
                    x=years[:len(pubs)], y=pubs,
                    title="Publication Trend",
                    template="plotly_dark",
                    color_discrete_sequence=["#6366f1"]
                )
                fig.update_layout(
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)',
                    xaxis_title="Year",
                    yaxis_title="Publications"
                )
                st.plotly_chart(fig, use_container_width=True)

                # Hot topics
                st.markdown("#### 🔥 Hot Topics")
                topics = {
                    "Anti-VEGF resistance": 89,
                    "Tumor microenvironment": 76,
                    "Immunotherapy combination": 72,
                    "Single-cell analysis": 65,
                    "Biomarkers": 58
                }

                fig = px.bar(
                    x=list(topics.values()),
                    y=list(topics.keys()),
                    orientation='h',
                    template="plotly_dark",
                    color=list(topics.values()),
                    color_continuous_scale="Viridis"
                )
                fig.update_layout(
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)',
                    showlegend=False,
                    yaxis_title=""
                )
                st.plotly_chart(fig, use_container_width=True)

            with tabs[1]:
                groups = pd.DataFrame({
                    "Institution": ["Harvard Medical School", "Stanford University", "NIH", "Genentech", "Memorial Sloan Kettering"],
                    "Publications": [45, 38, 32, 28, 25],
                    "Citations": [2340, 1890, 1560, 1450, 1280],
                    "Focus": ["Resistance mechanisms", "Biomarkers", "Basic biology", "Drug development", "Clinical trials"]
                })

                st.dataframe(groups, use_container_width=True, hide_index=True)

            with tabs[2]:
                patents = pd.DataFrame({
                    "Title": ["Anti-VEGF antibody composition", "VEGFR inhibitor combination", "Biomarker panel for response"],
                    "Assignee": ["Genentech", "Novartis", "Foundation Medicine"],
                    "Year": [2022, 2023, 2023],
                    "Status": ["Granted", "Pending", "Pending"]
                })

                st.dataframe(patents, use_container_width=True, hide_index=True)

            with tabs[3]:
                st.markdown("""
                <div class="glass-card-highlight">
                    <h4>🎯 Identified Research Gaps</h4>
                </div>
                """, unsafe_allow_html=True)

                gaps = [
                    {"gap": "Mechanism of acquired resistance in specific tumor types", "opportunity": "High", "competition": "Low"},
                    {"gap": "Novel biomarkers for patient selection", "opportunity": "High", "competition": "Medium"},
                    {"gap": "Combination with immunotherapy optimization", "opportunity": "Medium", "competition": "High"},
                    {"gap": "Single-cell level understanding of TME", "opportunity": "High", "competition": "Medium"},
                ]

                st.dataframe(pd.DataFrame(gaps), use_container_width=True, hide_index=True)


def render_figure_generator():
    """Render figure generator page"""
    st.markdown("<h1>🎨 Figure Generator</h1>", unsafe_allow_html=True)

    st.markdown("""
    <div class="glass-card">
        <p>Nature Reviews 수준의 pathway schematic을 자동 생성합니다.</p>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.markdown("### 🎯 Figure Settings")

        pathway_type = st.selectbox(
            "Pathway Template",
            ["VEGF-Notch Crosstalk", "EGFR-RAS-MAPK", "Wnt/β-catenin", "TGF-β/SMAD", "Custom"]
        )

        style = st.selectbox(
            "Style",
            ["Nature Reviews", "Cell", "Science", "Minimal"]
        )

        show_options = st.multiselect(
            "Include",
            ["Interactions", "Inhibitors", "KO effects", "Cell compartments"],
            default=["Interactions", "Inhibitors"]
        )

        dimensions = st.slider("Width (mm)", 80, 180, 150)

        if st.button("🎨 Generate Figure", use_container_width=True):
            st.session_state.figure_done = True

    with col2:
        if st.session_state.get('figure_done'):
            st.markdown("### 📊 Generated Figure Preview")

            # Create a mock pathway figure
            fig = go.Figure()

            # Define pathway components
            if "VEGF" in pathway_type:
                nodes = {
                    "VEGF-A": (0, 3, "Ligand"),
                    "VEGFR2": (0, 2, "Receptor"),
                    "PLCγ": (-1, 1, "Kinase"),
                    "ERK": (-1, 0, "Kinase"),
                    "DLL4": (1, 2, "Ligand"),
                    "Notch1": (2, 2, "Receptor"),
                    "HES1": (2, 1, "TF"),
                }
            else:
                nodes = {
                    "Ligand": (0, 3, "Ligand"),
                    "Receptor": (0, 2, "Receptor"),
                    "Kinase1": (0, 1, "Kinase"),
                    "TF": (0, 0, "TF"),
                }

            # Colors by type
            type_colors = {
                "Ligand": "#22c55e",
                "Receptor": "#6366f1",
                "Kinase": "#f59e0b",
                "TF": "#ef4444"
            }

            # Draw connections
            connections = [
                ("VEGF-A", "VEGFR2"),
                ("VEGFR2", "PLCγ"),
                ("PLCγ", "ERK"),
                ("VEGFR2", "DLL4"),
                ("DLL4", "Notch1"),
                ("Notch1", "HES1"),
            ]

            for src, tgt in connections:
                if src in nodes and tgt in nodes:
                    fig.add_trace(go.Scatter(
                        x=[nodes[src][0], nodes[tgt][0]],
                        y=[nodes[src][1], nodes[tgt][1]],
                        mode='lines',
                        line=dict(color='rgba(255,255,255,0.3)', width=2),
                        hoverinfo='none',
                        showlegend=False
                    ))

            # Draw nodes
            for name, (x, y, ntype) in nodes.items():
                fig.add_trace(go.Scatter(
                    x=[x], y=[y],
                    mode='markers+text',
                    marker=dict(size=40, color=type_colors.get(ntype, "#6366f1"), line=dict(color='white', width=2)),
                    text=[name],
                    textposition='top center',
                    textfont=dict(color='white', size=10),
                    name=ntype,
                    hovertemplate=f"<b>{name}</b><br>Type: {ntype}<extra></extra>"
                ))

            fig.update_layout(
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(26,26,46,0.8)',
                showlegend=True,
                legend=dict(
                    bgcolor='rgba(0,0,0,0.3)',
                    font=dict(color='white')
                ),
                xaxis=dict(visible=False, range=[-2, 4]),
                yaxis=dict(visible=False, range=[-0.5, 4]),
                height=500,
                margin=dict(l=20, r=20, t=20, b=20)
            )

            st.plotly_chart(fig, use_container_width=True)

            # Export options
            st.markdown("### 📥 Export Options")

            col_a, col_b, col_c = st.columns(3)
            with col_a:
                st.button("📄 Download SVG", use_container_width=True)
            with col_b:
                st.button("📑 Download PDF", use_container_width=True)
            with col_c:
                st.button("🖼️ Download PNG (600 DPI)", use_container_width=True)


if __name__ == "__main__":
    main()
