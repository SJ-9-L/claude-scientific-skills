# Competitive Landscape Analysis Reference

## Research Group Identification

### Literature-Based Analysis
```python
from collections import defaultdict
import requests
from Bio import Entrez

Entrez.email = "your@email.com"

def identify_research_groups(topic, years=5, min_publications=3):
    """
    Identify key research groups working on a topic.

    Returns:
    - Top authors/groups
    - Institution affiliations
    - Publication counts
    - Citation metrics
    - Collaboration networks
    """

    # Search PubMed
    query = f"{topic} AND {years}[PDAT]"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1000)
    results = Entrez.read(handle)

    # Fetch paper details
    papers = []
    for i in range(0, len(results['IdList']), 200):
        batch = results['IdList'][i:i+200]
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=','.join(batch),
            rettype="xml"
        )
        papers.extend(Entrez.read(fetch_handle)['PubmedArticle'])

    # Extract authors and affiliations
    author_stats = defaultdict(lambda: {
        'papers': [],
        'affiliations': set(),
        'first_author': 0,
        'last_author': 0
    })

    for paper in papers:
        article = paper['MedlineCitation']['Article']
        pmid = paper['MedlineCitation']['PMID']

        if 'AuthorList' in article:
            authors = article['AuthorList']
            for i, author in enumerate(authors):
                if 'LastName' in author and 'ForeName' in author:
                    name = f"{author['LastName']} {author['ForeName']}"
                    author_stats[name]['papers'].append(pmid)

                    if 'AffiliationInfo' in author:
                        for aff in author['AffiliationInfo']:
                            author_stats[name]['affiliations'].add(
                                aff.get('Affiliation', '')
                            )

                    if i == 0:
                        author_stats[name]['first_author'] += 1
                    if i == len(authors) - 1:
                        author_stats[name]['last_author'] += 1

    # Filter and rank
    key_researchers = []
    for name, stats in author_stats.items():
        if len(stats['papers']) >= min_publications:
            key_researchers.append({
                'name': name,
                'publication_count': len(stats['papers']),
                'first_author': stats['first_author'],
                'last_author': stats['last_author'],
                'affiliations': list(stats['affiliations']),
                'pmids': stats['papers']
            })

    # Sort by publication count
    key_researchers.sort(key=lambda x: x['publication_count'], reverse=True)

    return key_researchers
```

### Institution Analysis
```python
def analyze_institutions(researchers):
    """Group researchers by institution."""
    institution_map = defaultdict(list)

    # Common institution patterns
    keywords = {
        'harvard': 'Harvard University',
        'mit': 'MIT',
        'stanford': 'Stanford University',
        'nih': 'NIH',
        'max planck': 'Max Planck Institute',
        'karolinska': 'Karolinska Institute',
        'cambridge': 'University of Cambridge',
        'oxford': 'University of Oxford',
        'yale': 'Yale University',
        'johns hopkins': 'Johns Hopkins University'
    }

    for researcher in researchers:
        for aff in researcher['affiliations']:
            aff_lower = aff.lower()
            for keyword, institution in keywords.items():
                if keyword in aff_lower:
                    institution_map[institution].append(researcher['name'])
                    break

    return dict(institution_map)
```

## Trend Analysis

### Publication Trends
```python
def analyze_publication_trends(topic, start_year=2015, end_year=2024):
    """Analyze publication trends over time."""
    trends = {}

    for year in range(start_year, end_year + 1):
        query = f"{topic} AND {year}[PDAT]"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        result = Entrez.read(handle)
        trends[year] = int(result['Count'])

    return trends

def analyze_keyword_emergence(topic, keywords, years=5):
    """Track emergence of new keywords/concepts."""
    from datetime import datetime
    current_year = datetime.now().year

    keyword_trends = {}

    for keyword in keywords:
        keyword_trends[keyword] = {}
        for year in range(current_year - years, current_year + 1):
            query = f"{topic} AND {keyword} AND {year}[PDAT]"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
            result = Entrez.read(handle)
            keyword_trends[keyword][year] = int(result['Count'])

    return keyword_trends
```

### Technology/Method Trends
```python
technology_trends = {
    'single_cell': {
        '2015-2017': ['Smart-seq', 'Drop-seq', 'inDrop'],
        '2018-2020': ['10X Genomics', 'sci-RNA-seq', 'CITE-seq'],
        '2021-2023': ['Spatial transcriptomics', 'Multiome', 'Parse Biosciences'],
        '2024+': ['Perturb-seq at scale', 'Long-read scRNA-seq']
    },
    'analysis_methods': {
        '2015-2017': ['t-SNE', 'Seurat v1-2', 'Monocle'],
        '2018-2020': ['UMAP', 'Seurat v3', 'Scanpy', 'scVI'],
        '2021-2023': ['CellTypist', 'scArches', 'CellRank', 'Squidpy'],
        '2024+': ['Foundation models', 'LLM integration', 'Multi-modal']
    }
}
```

## Patent Analysis

### USPTO Patent Search
```python
def search_patents(query, max_results=100):
    """Search USPTO patents."""
    # Using PatentsView API
    base_url = "https://api.patentsview.org/patents/query"

    params = {
        'q': {"_text_any": {"patent_abstract": query}},
        'f': [
            "patent_number",
            "patent_title",
            "patent_abstract",
            "patent_date",
            "assignee_organization",
            "inventor_first_name",
            "inventor_last_name"
        ],
        'o': {"per_page": max_results}
    }

    response = requests.post(base_url, json=params)
    return response.json()

def analyze_patent_landscape(topic):
    """Analyze patent landscape for a topic."""
    patents = search_patents(topic)

    # Analyze assignees (companies/institutions)
    assignee_counts = defaultdict(int)
    for patent in patents.get('patents', []):
        for assignee in patent.get('assignees', []):
            org = assignee.get('assignee_organization', 'Unknown')
            assignee_counts[org] += 1

    # Analyze trends by year
    year_counts = defaultdict(int)
    for patent in patents.get('patents', []):
        date = patent.get('patent_date', '')
        if date:
            year = date[:4]
            year_counts[year] += 1

    return {
        'total_patents': len(patents.get('patents', [])),
        'top_assignees': dict(sorted(
            assignee_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:20]),
        'yearly_trend': dict(year_counts)
    }
```

## Research Gap Identification

### Systematic Gap Analysis
```python
def identify_research_gaps(topic, known_aspects):
    """
    Identify under-researched areas within a topic.

    Parameters:
    -----------
    topic : str
        Main research topic
    known_aspects : list
        Known sub-topics to check coverage

    Returns:
    --------
    dict with:
    - well_covered: High publication count
    - emerging: Recent increase in publications
    - gaps: Under-researched areas
    """

    results = {
        'well_covered': [],
        'emerging': [],
        'gaps': []
    }

    for aspect in known_aspects:
        query = f"{topic} AND {aspect}"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        count = int(Entrez.read(handle)['Count'])

        # Check recent trend
        recent_query = f"{topic} AND {aspect} AND 2022:2024[PDAT]"
        recent_handle = Entrez.esearch(db="pubmed", term=recent_query, retmax=0)
        recent_count = int(Entrez.read(recent_handle)['Count'])

        if count > 100:
            results['well_covered'].append({
                'aspect': aspect,
                'total': count,
                'recent': recent_count
            })
        elif count < 20:
            results['gaps'].append({
                'aspect': aspect,
                'total': count,
                'recent': recent_count
            })
        elif recent_count > count * 0.4:  # >40% in last 2 years
            results['emerging'].append({
                'aspect': aspect,
                'total': count,
                'recent': recent_count
            })

    return results
```

### Unmet Needs Analysis
```python
unmet_needs_framework = {
    'categories': {
        'mechanism': 'Unknown molecular mechanisms',
        'target': 'Undruggable targets',
        'biomarker': 'Lack of predictive biomarkers',
        'resistance': 'Treatment resistance mechanisms',
        'translation': 'Preclinical to clinical gap',
        'heterogeneity': 'Cellular/tumor heterogeneity',
        'delivery': 'Drug delivery challenges',
        'toxicity': 'Off-target toxicity concerns'
    },

    'assessment_queries': {
        'mechanism': '"{topic}" AND (mechanism OR pathway) AND (unknown OR unclear OR elusive)',
        'target': '"{topic}" AND (target OR therapeutic) AND (undruggable OR challenging)',
        'biomarker': '"{topic}" AND biomarker AND (needed OR lacking OR predictive)',
        'resistance': '"{topic}" AND (resistance OR refractory OR relapse)',
        'translation': '"{topic}" AND (clinical trial OR phase) AND (failure OR discontinued)',
        'heterogeneity': '"{topic}" AND (heterogeneity OR subpopulation OR variability)'
    }
}

def assess_unmet_needs(topic):
    """Assess unmet needs in a research area."""
    needs = {}

    for category, query_template in unmet_needs_framework['assessment_queries'].items():
        query = query_template.format(topic=topic)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        count = int(Entrez.read(handle)['Count'])

        needs[category] = {
            'description': unmet_needs_framework['categories'][category],
            'literature_count': count,
            'priority': 'high' if count > 50 else 'medium' if count > 10 else 'low'
        }

    return needs
```

## Competitive Report Generation

### Report Template
```python
def generate_competitive_report(topic, analysis_results):
    """Generate comprehensive competitive landscape report."""

    report = f"""
# Competitive Landscape Report: {topic}

## Executive Summary
- Total publications analyzed: {analysis_results['total_papers']}
- Key research groups identified: {len(analysis_results['key_researchers'])}
- Active institutions: {len(analysis_results['institutions'])}
- Patent filings: {analysis_results['patents']['total_patents']}

## Key Research Groups

### Top 10 Most Productive Authors
| Rank | Researcher | Publications | First Author | Last Author | Institution |
|------|-----------|--------------|--------------|-------------|-------------|
"""

    for i, researcher in enumerate(analysis_results['key_researchers'][:10], 1):
        institution = researcher['affiliations'][0][:50] if researcher['affiliations'] else 'N/A'
        report += f"| {i} | {researcher['name']} | {researcher['publication_count']} | {researcher['first_author']} | {researcher['last_author']} | {institution} |\n"

    report += f"""
## Publication Trends

### Annual Publication Count
"""

    for year, count in sorted(analysis_results['trends'].items()):
        report += f"- {year}: {count} publications\n"

    report += f"""
## Research Gaps and Opportunities

### Identified Gaps
"""

    for gap in analysis_results['gaps']['gaps']:
        report += f"- **{gap['aspect']}**: Only {gap['total']} publications (potential opportunity)\n"

    report += f"""
### Emerging Areas
"""

    for emerging in analysis_results['gaps']['emerging']:
        report += f"- **{emerging['aspect']}**: {emerging['recent']} recent publications (growing interest)\n"

    report += f"""
## Patent Landscape

### Top Patent Holders
"""

    for assignee, count in list(analysis_results['patents']['top_assignees'].items())[:10]:
        report += f"- {assignee}: {count} patents\n"

    report += f"""
## Unmet Needs

| Category | Description | Priority |
|----------|-------------|----------|
"""

    for category, details in analysis_results['unmet_needs'].items():
        report += f"| {category} | {details['description']} | {details['priority']} |\n"

    return report
```

## Citation Network Analysis

### Build Citation Network
```python
import networkx as nx

def build_citation_network(seed_pmids, depth=2):
    """
    Build citation network from seed papers.

    Parameters:
    -----------
    seed_pmids : list
        Starting paper PMIDs
    depth : int
        How many levels of citations to follow

    Returns:
    --------
    NetworkX graph with papers as nodes and citations as edges
    """

    G = nx.DiGraph()

    # Add seed papers
    for pmid in seed_pmids:
        G.add_node(pmid, level=0)

    # Fetch citations iteratively
    for level in range(depth):
        current_level_nodes = [
            n for n, d in G.nodes(data=True)
            if d.get('level') == level
        ]

        for pmid in current_level_nodes:
            # Get papers citing this paper
            link_handle = Entrez.elink(
                dbfrom="pubmed",
                db="pubmed",
                id=pmid,
                linkname="pubmed_pubmed_citedin"
            )
            links = Entrez.read(link_handle)

            if links[0]['LinkSetDb']:
                cited_by = [
                    link['Id']
                    for link in links[0]['LinkSetDb'][0]['Link']
                ]

                for citing_pmid in cited_by[:50]:  # Limit to top 50
                    if citing_pmid not in G:
                        G.add_node(citing_pmid, level=level + 1)
                    G.add_edge(pmid, citing_pmid)

    return G

def identify_hub_papers(citation_graph):
    """Identify influential papers in citation network."""

    # Calculate centrality metrics
    pagerank = nx.pagerank(citation_graph)
    in_degree = dict(citation_graph.in_degree())
    out_degree = dict(citation_graph.out_degree())

    # Combine metrics
    hub_scores = {}
    for node in citation_graph.nodes():
        hub_scores[node] = {
            'pagerank': pagerank.get(node, 0),
            'citations': in_degree.get(node, 0),
            'references': out_degree.get(node, 0)
        }

    # Rank by PageRank
    ranked = sorted(
        hub_scores.items(),
        key=lambda x: x[1]['pagerank'],
        reverse=True
    )

    return ranked[:20]
```

## Collaboration Network Analysis

```python
def build_collaboration_network(papers):
    """Build co-authorship network from papers."""

    G = nx.Graph()

    for paper in papers:
        article = paper['MedlineCitation']['Article']

        if 'AuthorList' in article:
            authors = article['AuthorList']
            author_names = []

            for author in authors:
                if 'LastName' in author and 'ForeName' in author:
                    name = f"{author['LastName']} {author['ForeName']}"
                    author_names.append(name)

                    if name not in G:
                        G.add_node(name, paper_count=0)
                    G.nodes[name]['paper_count'] += 1

            # Add edges between co-authors
            for i, author1 in enumerate(author_names):
                for author2 in author_names[i+1:]:
                    if G.has_edge(author1, author2):
                        G[author1][author2]['weight'] += 1
                    else:
                        G.add_edge(author1, author2, weight=1)

    return G

def find_research_communities(collaboration_graph):
    """Identify research communities/clusters."""
    from networkx.algorithms import community

    # Detect communities
    communities = community.louvain_communities(collaboration_graph)

    # Analyze each community
    community_info = []
    for i, comm in enumerate(communities):
        members = list(comm)
        subgraph = collaboration_graph.subgraph(members)

        # Find most central member
        centrality = nx.degree_centrality(subgraph)
        leader = max(centrality.items(), key=lambda x: x[1])[0]

        community_info.append({
            'id': i,
            'size': len(members),
            'leader': leader,
            'members': members[:10],  # Top 10
            'density': nx.density(subgraph)
        })

    return community_info
```
