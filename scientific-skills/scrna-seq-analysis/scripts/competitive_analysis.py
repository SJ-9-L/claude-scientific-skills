#!/usr/bin/env python3
"""
Competitive Landscape Analysis Module

Research group identification, trend analysis, patent landscape,
and research gap identification for scientific competitive intelligence.
"""

import time
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict
from datetime import datetime, timedelta
import json

import requests


@dataclass
class ResearchGroup:
    """Represents a research group/lab."""
    name: str
    institution: str
    publication_count: int = 0
    first_author_count: int = 0
    last_author_count: int = 0
    h_index: Optional[int] = None
    recent_focus: List[str] = field(default_factory=list)
    collaborators: List[str] = field(default_factory=list)
    pmids: List[str] = field(default_factory=list)


@dataclass
class TrendData:
    """Publication and keyword trends."""
    yearly_counts: Dict[int, int] = field(default_factory=dict)
    keyword_emergence: Dict[str, Dict[int, int]] = field(default_factory=dict)
    hot_topics: List[str] = field(default_factory=list)
    declining_topics: List[str] = field(default_factory=list)


@dataclass
class PatentLandscape:
    """Patent analysis results."""
    total_patents: int = 0
    top_assignees: Dict[str, int] = field(default_factory=dict)
    yearly_trend: Dict[int, int] = field(default_factory=dict)
    key_patents: List[Dict] = field(default_factory=list)


@dataclass
class ResearchGap:
    """Identified research gap."""
    topic: str
    description: str
    publication_count: int
    priority: str  # high, medium, low
    suggested_experiments: List[str] = field(default_factory=list)


class CompetitiveLandscape:
    """
    Competitive landscape analysis for research topics.

    Features:
    - Key research group identification
    - Publication trend analysis
    - Patent landscape analysis
    - Research gap identification
    - Report generation
    """

    def __init__(
        self,
        topic: str,
        timeframe: str = '2020-2024'
    ):
        """
        Initialize competitive analysis.

        Parameters
        ----------
        topic : str
            Research topic to analyze
        timeframe : str
            Analysis timeframe (e.g., '2020-2024')
        """
        self.topic = topic
        self.timeframe = timeframe

        # Parse timeframe
        years = timeframe.split('-')
        self.start_year = int(years[0])
        self.end_year = int(years[1]) if len(years) > 1 else datetime.now().year

        self.research_groups: List[ResearchGroup] = []
        self.trends: Optional[TrendData] = None
        self.patents: Optional[PatentLandscape] = None
        self.gaps: List[ResearchGap] = []

    def identify_key_groups(
        self,
        min_publications: int = 5,
        include_industry: bool = True
    ) -> List[ResearchGroup]:
        """
        Identify key research groups working on the topic.

        Parameters
        ----------
        min_publications : int
            Minimum publications to be considered key group
        include_industry : bool
            Include industry affiliations

        Returns
        -------
        list
            Key research groups
        """
        # Query PubMed
        papers = self._search_pubmed(
            f"{self.topic} {self.start_year}:{self.end_year}[PDAT]",
            max_results=500
        )

        # Extract author statistics
        author_stats = self._extract_author_stats(papers)

        # Filter and create research groups
        self.research_groups = []
        for author, stats in author_stats.items():
            if stats['total'] >= min_publications:
                group = ResearchGroup(
                    name=author,
                    institution=stats['affiliations'][0] if stats['affiliations'] else 'Unknown',
                    publication_count=stats['total'],
                    first_author_count=stats['first_author'],
                    last_author_count=stats['last_author'],
                    pmids=stats['pmids']
                )
                self.research_groups.append(group)

        # Sort by publication count
        self.research_groups.sort(key=lambda x: x.publication_count, reverse=True)

        return self.research_groups[:50]

    def _search_pubmed(self, query: str, max_results: int = 100) -> List[Dict]:
        """Search PubMed using Entrez."""
        try:
            from Bio import Entrez
            Entrez.email = "analysis@research.com"

            # Search
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            results = Entrez.read(handle)
            pmids = results.get('IdList', [])

            if not pmids:
                return []

            # Fetch details
            papers = []
            batch_size = 100
            for i in range(0, len(pmids), batch_size):
                batch = pmids[i:i + batch_size]
                fetch_handle = Entrez.efetch(
                    db="pubmed",
                    id=','.join(batch),
                    rettype="xml"
                )
                batch_results = Entrez.read(fetch_handle)
                papers.extend(batch_results.get('PubmedArticle', []))
                time.sleep(0.5)  # Rate limiting

            return papers

        except ImportError:
            print("Biopython not available. Using mock data.")
            return []
        except Exception as e:
            print(f"PubMed search error: {e}")
            return []

    def _extract_author_stats(self, papers: List[Dict]) -> Dict[str, Dict]:
        """Extract author statistics from papers."""
        author_stats = defaultdict(lambda: {
            'total': 0,
            'first_author': 0,
            'last_author': 0,
            'affiliations': [],
            'pmids': []
        })

        for paper in papers:
            try:
                article = paper.get('MedlineCitation', {}).get('Article', {})
                pmid = str(paper.get('MedlineCitation', {}).get('PMID', ''))
                author_list = article.get('AuthorList', [])

                for i, author in enumerate(author_list):
                    if 'LastName' in author and 'ForeName' in author:
                        name = f"{author['LastName']} {author['ForeName']}"
                        author_stats[name]['total'] += 1
                        author_stats[name]['pmids'].append(pmid)

                        # Track position
                        if i == 0:
                            author_stats[name]['first_author'] += 1
                        if i == len(author_list) - 1:
                            author_stats[name]['last_author'] += 1

                        # Extract affiliation
                        if 'AffiliationInfo' in author:
                            for aff in author['AffiliationInfo']:
                                aff_text = aff.get('Affiliation', '')
                                if aff_text and aff_text not in author_stats[name]['affiliations']:
                                    author_stats[name]['affiliations'].append(aff_text)

            except Exception as e:
                continue

        return dict(author_stats)

    def analyze_trends(
        self,
        include_preprints: bool = True,
        include_patents: bool = True,
        keywords: Optional[List[str]] = None
    ) -> TrendData:
        """
        Analyze publication and keyword trends.

        Parameters
        ----------
        include_preprints : bool
            Include bioRxiv/medRxiv
        include_patents : bool
            Include patent analysis
        keywords : list, optional
            Specific keywords to track

        Returns
        -------
        TrendData
            Trend analysis results
        """
        self.trends = TrendData()

        # Yearly publication counts
        for year in range(self.start_year, self.end_year + 1):
            count = self._get_pubmed_count(f"{self.topic} {year}[PDAT]")
            self.trends.yearly_counts[year] = count

        # Keyword emergence
        if keywords:
            for keyword in keywords:
                self.trends.keyword_emergence[keyword] = {}
                for year in range(self.start_year, self.end_year + 1):
                    count = self._get_pubmed_count(
                        f"{self.topic} {keyword} {year}[PDAT]"
                    )
                    self.trends.keyword_emergence[keyword][year] = count

        # Identify hot topics (increasing trend)
        self._identify_hot_topics()

        # Preprints (bioRxiv)
        if include_preprints:
            preprint_data = self._search_biorxiv()
            # Merge with trends

        # Patents
        if include_patents:
            self.patents = self._analyze_patents()

        return self.trends

    def _get_pubmed_count(self, query: str) -> int:
        """Get publication count for a query."""
        try:
            from Bio import Entrez
            Entrez.email = "analysis@research.com"

            handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
            result = Entrez.read(handle)
            return int(result.get('Count', 0))

        except ImportError:
            return 0
        except Exception:
            return 0

    def _search_biorxiv(self) -> List[Dict]:
        """Search bioRxiv for preprints."""
        try:
            end_date = datetime.now().strftime('%Y-%m-%d')
            start_date = (datetime.now() - timedelta(days=365)).strftime('%Y-%m-%d')

            url = f"https://api.biorxiv.org/details/biorxiv/{start_date}/{end_date}"
            response = requests.get(url, timeout=30)

            if response.status_code == 200:
                data = response.json()
                papers = data.get('collection', [])

                # Filter by topic
                filtered = [
                    p for p in papers
                    if self.topic.lower() in p.get('title', '').lower() or
                       self.topic.lower() in p.get('abstract', '').lower()
                ]
                return filtered

        except Exception as e:
            print(f"bioRxiv search error: {e}")

        return []

    def _identify_hot_topics(self):
        """Identify trending and declining topics."""
        if not self.trends or not self.trends.yearly_counts:
            return

        years = sorted(self.trends.yearly_counts.keys())
        if len(years) < 3:
            return

        # Calculate growth rate
        recent = sum(self.trends.yearly_counts.get(y, 0) for y in years[-2:])
        early = sum(self.trends.yearly_counts.get(y, 0) for y in years[:2])

        if early > 0:
            growth_rate = (recent - early) / early
            if growth_rate > 0.5:
                self.trends.hot_topics.append(self.topic)
            elif growth_rate < -0.2:
                self.trends.declining_topics.append(self.topic)

        # Analyze keyword trends
        for keyword, yearly in self.trends.keyword_emergence.items():
            if not yearly:
                continue

            kw_years = sorted(yearly.keys())
            if len(kw_years) >= 2:
                recent_kw = yearly.get(kw_years[-1], 0)
                early_kw = yearly.get(kw_years[0], 0)

                if early_kw > 0 and recent_kw / early_kw > 2:
                    self.trends.hot_topics.append(keyword)
                elif recent_kw < early_kw * 0.5:
                    self.trends.declining_topics.append(keyword)

    def _analyze_patents(self) -> PatentLandscape:
        """Analyze patent landscape."""
        landscape = PatentLandscape()

        try:
            # Using PatentsView API
            base_url = "https://api.patentsview.org/patents/query"

            query = {
                "q": {"_text_any": {"patent_abstract": self.topic}},
                "f": [
                    "patent_number",
                    "patent_title",
                    "patent_date",
                    "assignee_organization"
                ],
                "o": {"per_page": 100}
            }

            response = requests.post(base_url, json=query, timeout=30)

            if response.status_code == 200:
                data = response.json()
                patents = data.get('patents', [])
                landscape.total_patents = data.get('count', len(patents))

                # Analyze assignees
                assignee_counts = defaultdict(int)
                year_counts = defaultdict(int)

                for patent in patents:
                    # Assignee
                    for assignee in patent.get('assignees', []):
                        org = assignee.get('assignee_organization', 'Unknown')
                        if org:
                            assignee_counts[org] += 1

                    # Year
                    date = patent.get('patent_date', '')
                    if date:
                        year = int(date[:4])
                        year_counts[year] += 1

                landscape.top_assignees = dict(sorted(
                    assignee_counts.items(),
                    key=lambda x: x[1],
                    reverse=True
                )[:20])

                landscape.yearly_trend = dict(year_counts)

                # Key patents
                landscape.key_patents = [
                    {
                        'number': p.get('patent_number', ''),
                        'title': p.get('patent_title', ''),
                        'date': p.get('patent_date', '')
                    }
                    for p in patents[:10]
                ]

        except Exception as e:
            print(f"Patent analysis error: {e}")

        return landscape

    def find_research_gaps(
        self,
        current_knowledge: Optional[Any] = None,
        identify_opportunities: bool = True
    ) -> List[ResearchGap]:
        """
        Identify under-researched areas and opportunities.

        Parameters
        ----------
        current_knowledge : Any, optional
            Current pathway/knowledge data
        identify_opportunities : bool
            Suggest experimental opportunities

        Returns
        -------
        list
            Research gaps
        """
        self.gaps = []

        # Define sub-topics to check
        sub_topics = [
            'mechanism',
            'biomarker',
            'resistance',
            'combination therapy',
            'prediction',
            'clinical translation',
            'heterogeneity',
            'single-cell',
            'spatial',
            'AI/machine learning'
        ]

        for sub_topic in sub_topics:
            query = f"{self.topic} {sub_topic}"
            count = self._get_pubmed_count(f"{query} {self.start_year}:{self.end_year}[PDAT]")

            # Determine priority based on count
            if count < 10:
                priority = 'high'
            elif count < 50:
                priority = 'medium'
            else:
                priority = 'low'

            gap = ResearchGap(
                topic=sub_topic,
                description=f"Research on {sub_topic} in {self.topic}",
                publication_count=count,
                priority=priority
            )

            # Suggest experiments for high-priority gaps
            if identify_opportunities and priority == 'high':
                gap.suggested_experiments = self._suggest_experiments(sub_topic)

            if count < 50:  # Only include if relatively under-researched
                self.gaps.append(gap)

        # Sort by priority and count
        priority_order = {'high': 0, 'medium': 1, 'low': 2}
        self.gaps.sort(key=lambda x: (priority_order[x.priority], x.publication_count))

        return self.gaps

    def _suggest_experiments(self, topic: str) -> List[str]:
        """Suggest experiments for a research gap."""
        suggestions = {
            'mechanism': [
                'CRISPR screen for pathway components',
                'Proteomics analysis of protein interactions',
                'Single-cell perturbation studies'
            ],
            'biomarker': [
                'Multi-omics biomarker discovery',
                'Clinical sample validation study',
                'Longitudinal sampling analysis'
            ],
            'resistance': [
                'Drug-resistant cell line generation',
                'Combination screening',
                'Compensatory pathway analysis'
            ],
            'single-cell': [
                'scRNA-seq time course',
                'Multi-modal single-cell profiling',
                'Spatial transcriptomics'
            ],
            'heterogeneity': [
                'Single-cell heterogeneity mapping',
                'Clonal evolution tracking',
                'Microenvironment analysis'
            ]
        }

        return suggestions.get(topic, [
            'Literature deep-dive and hypothesis generation',
            'Pilot experimental validation',
            'Collaborative study design'
        ])

    def generate_report(
        self,
        format: str = 'detailed',
        include_visualizations: bool = True,
        output_path: str = 'competitive_report.md'
    ) -> str:
        """
        Generate competitive landscape report.

        Parameters
        ----------
        format : str
            Report format: 'detailed', 'summary', 'executive'
        include_visualizations : bool
            Include visualization data
        output_path : str
            Output file path

        Returns
        -------
        str
            Report content or file path
        """
        report = f"""# Competitive Landscape Report
## Topic: {self.topic}
## Timeframe: {self.timeframe}
## Generated: {datetime.now().strftime('%Y-%m-%d')}

---

## Executive Summary

This report analyzes the competitive landscape for research on **{self.topic}** from {self.start_year} to {self.end_year}.

### Key Findings
- **Active Research Groups**: {len(self.research_groups)}
- **Total Patents**: {self.patents.total_patents if self.patents else 'N/A'}
- **Research Gaps Identified**: {len(self.gaps)}

---

## Key Research Groups

| Rank | Researcher | Publications | First Author | Last Author | Institution |
|------|-----------|--------------|--------------|-------------|-------------|
"""
        for i, group in enumerate(self.research_groups[:15], 1):
            institution = group.institution[:40] + '...' if len(group.institution) > 40 else group.institution
            report += f"| {i} | {group.name} | {group.publication_count} | {group.first_author_count} | {group.last_author_count} | {institution} |\n"

        report += f"""
---

## Publication Trends

### Annual Publication Counts
"""
        if self.trends:
            for year in sorted(self.trends.yearly_counts.keys()):
                count = self.trends.yearly_counts[year]
                bar = '█' * (count // 10)
                report += f"- **{year}**: {count} publications {bar}\n"

            if self.trends.hot_topics:
                report += f"\n### Hot Topics (Increasing Trend)\n"
                for topic in self.trends.hot_topics:
                    report += f"- {topic}\n"

            if self.trends.declining_topics:
                report += f"\n### Declining Topics\n"
                for topic in self.trends.declining_topics:
                    report += f"- {topic}\n"

        report += f"""
---

## Patent Landscape
"""
        if self.patents:
            report += f"\n**Total Patents**: {self.patents.total_patents}\n\n"
            report += "### Top Patent Holders\n"
            for org, count in list(self.patents.top_assignees.items())[:10]:
                report += f"- {org}: {count} patents\n"

        report += f"""
---

## Research Gaps and Opportunities

| Priority | Topic | Publications | Description |
|----------|-------|--------------|-------------|
"""
        for gap in self.gaps:
            report += f"| {gap.priority.upper()} | {gap.topic} | {gap.publication_count} | {gap.description} |\n"

        report += f"""
### High-Priority Experimental Opportunities
"""
        for gap in self.gaps:
            if gap.priority == 'high' and gap.suggested_experiments:
                report += f"\n**{gap.topic}**:\n"
                for exp in gap.suggested_experiments:
                    report += f"- {exp}\n"

        report += f"""
---

## Recommendations

1. **Collaboration Opportunities**: Consider reaching out to top research groups for potential collaboration
2. **Gap Exploitation**: Focus on high-priority research gaps for novel contributions
3. **IP Strategy**: Monitor patent landscape for freedom-to-operate analysis
4. **Trend Surfing**: Invest in hot topics for maximum impact

---

*Report generated by Competitive Landscape Analysis Module*
"""

        # Save report
        with open(output_path, 'w') as f:
            f.write(report)

        print(f"Report saved to {output_path}")
        return report


class LiteratureAnalyzer:
    """Detailed literature analysis."""

    def __init__(self):
        self.papers = []

    def search(
        self,
        query: str,
        databases: List[str] = None,
        date_range: str = '2020-2024',
        filters: Optional[Dict] = None
    ) -> List[Dict]:
        """Search multiple literature databases."""
        if databases is None:
            databases = ['pubmed']

        all_papers = []

        if 'pubmed' in databases:
            pubmed_papers = self._search_pubmed(query, date_range)
            all_papers.extend(pubmed_papers)

        self.papers = all_papers
        return all_papers

    def _search_pubmed(self, query: str, date_range: str) -> List[Dict]:
        """Search PubMed."""
        try:
            from Bio import Entrez
            Entrez.email = "analysis@research.com"

            years = date_range.split('-')
            full_query = f"{query} {years[0]}:{years[1]}[PDAT]"

            handle = Entrez.esearch(db="pubmed", term=full_query, retmax=100)
            results = Entrez.read(handle)

            return [{'pmid': pmid, 'source': 'pubmed'} for pmid in results.get('IdList', [])]

        except ImportError:
            return []

    def identify_research_groups(
        self,
        papers: List[Dict],
        min_publications: int = 3,
        include_affiliations: bool = True
    ) -> List[Dict]:
        """Identify research groups from papers."""
        # Implementation similar to CompetitiveLandscape
        return []

    def analyze_citations(
        self,
        seed_papers: List[Dict],
        depth: int = 2
    ) -> Dict:
        """Analyze citation network."""
        return {'nodes': [], 'edges': []}

    def identify_gaps(
        self,
        current_literature: List[Dict],
        known_pathways: List[str],
        suggest_experiments: bool = True
    ) -> List[Dict]:
        """Identify research gaps."""
        return []


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Competitive Landscape Analysis')
    parser.add_argument('topic', help='Research topic')
    parser.add_argument('--timeframe', default='2020-2024', help='Analysis timeframe')
    parser.add_argument('--output', '-o', default='competitive_report.md', help='Output file')

    args = parser.parse_args()

    analyzer = CompetitiveLandscape(
        topic=args.topic,
        timeframe=args.timeframe
    )

    print("Identifying key research groups...")
    groups = analyzer.identify_key_groups()

    print("Analyzing trends...")
    trends = analyzer.analyze_trends()

    print("Finding research gaps...")
    gaps = analyzer.find_research_gaps()

    print("Generating report...")
    analyzer.generate_report(output_path=args.output)

    print(f"Analysis complete. Report saved to {args.output}")
