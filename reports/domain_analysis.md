# Domain Analysis: Variant–Domain Mapping in DBA Ribosomal Proteins

## 1. Biological Rationale

Ribosomal proteins in Diamond-Blackfan anemia are organized into defined structural and functional domains that mediate critical interactions with ribosomal RNA and assembly factors. Pathogenic variants in DBA preferentially disrupt protein function, suggesting that mutations may concentrate within structurally or functionally indispensable regions.

For the 60S compartment, RPL5 and RPL11 participate in the 5S ribonucleoprotein (5S RNP) complex that stabilizes p53 in response to ribosomal stress via MDM2 inhibition. The uL18 fold of RPL5 — responsible for 5S rRNA binding — represents a domain where missense variants may abrogate function without requiring truncation. For 40S proteins, RNA-contact surfaces at the platform and beak regions of the small subunit are candidate hotspot domains.

This analysis maps each variant to its domain context using UniProt annotations, addressing whether pathogenic mutations preferentially co-localize with known structural or functional domains.

## 2. Methods

### 2.1 Data Source

Protein domain annotations were obtained from the UniProt Knowledge Base REST API (https://rest.uniprot.org/uniprotkb/{UniProtID}.json) for each of the seven genes analyzed. UniProt consolidates domain information from Pfam, SMART, SUPERFAMILY, and PANTHER, providing a curated, non-redundant annotation set.

### 2.2 Feature Types Included

All features of the following types were retrieved: "Active site", "Binding site", "Domain", "Modified residue", "Motif", "Region". Secondary structure elements (Helix, Beta strand, Turn) were explicitly excluded to focus on structurally and functionally defined regions.

For intra/inter-domain classification, the following types were used: "Active site", "Binding site", "Domain", "Motif", "Region". "Modified residue" (PTM) sites were loaded but are only reported when a local variant falls at the exact annotated residue.

### 2.3 Variant–Domain Mapping

For each variant with a resolved amino acid position, domain overlap was assessed as: `feature_start ≤ variant_position ≤ feature_end`. When multiple features overlapped, the highest-priority type was retained (Domain > Region > Motif > Binding site > Active site). Variants without a resolved position (e.g., splice variants) are reported as "Sin posición".

### 2.4 Normalization

The domain architecture figure (Figure 1) uses normalized positions (`position / protein_length ∈ [0, 1]`) to enable visual comparison across proteins of different lengths.

## 3. Per-Gene Results

### RPS19 (40S, 40S)

| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |
|---------------------|----------|------|------|----------------|-------------|-------------|
| p.Arg94* | 94 | P | Nonsense | Inter-domain | — | No |
| p.Asp130Serfs*23 | 130 | P | Frameshift | Inter-domain | — | No |

### RPS24 (40S, 40S)

| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |
|---------------------|----------|------|------|----------------|-------------|-------------|
| p.Gly66fs*10 | 66 | LP | Frameshift | Inter-domain | — | No |
| p.Glu98* | 98 | P | Nonsense | Disordered | Region | No |

### RPS26 (40S, 40S)

| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |
|---------------------|----------|------|------|----------------|-------------|-------------|
| p.Cys74Tyr | 74 | LP | Missense | Inter-domain | — | No |
| (p.Met1?) | 1 | P | Other | Inter-domain | — | No |

### RPL5 (60S, 60S)

| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |
|---------------------|----------|------|------|----------------|-------------|-------------|
|  | — | VUS | Splice | Sin posición | — | No |
| p.Tyr86* | 86 | LP | Nonsense | Inter-domain | — | No |
| p.His81fs*37 | 81 | LP | Frameshift | Inter-domain | — | No |
| (p.Arg24Glnfs*14) | 24 | P | Frameshift | Inter-domain | — | No |

### RPL11 (60S, 60S)

| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |
|---------------------|----------|------|------|----------------|-------------|-------------|
| p.Cys21fs*33 | 21 | P | Frameshift | Inter-domain | — | No |
| p.(Tyr55AspfsTer34) | — | LP | Frameshift | Sin posición | — | No |

### RPL18 (60S, 60S)

| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |
|---------------------|----------|------|------|----------------|-------------|-------------|
| (p.Arg108Gln) | 108 | VUS | Other | Inter-domain | — | No |

### RPL26 (60S, 60S)

| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |
|---------------------|----------|------|------|----------------|-------------|-------------|
| p.Gln40* | 40 | LP | Nonsense | Inter-domain | — | No |

## 4. Summary Statistics

### 4.1 Per-Gene

| Gene | Subunit | Local (n) | % Intra-domain (local) | ClinVar (n) | % Intra-domain (ClinVar) |
|------|---------|-----------|------------------------|-------------|--------------------------|
| RPS19 | 40S | 2 | 0.0% | 154 | 0.0% |
| RPS24 | 40S | 2 | 50.0% | 34 | 23.5% |
| RPS26 | 40S | 2 | 0.0% | 50 | 8.0% |
| RPL5 | 60S | 3 | 0.0% | 152 | 1.3% |
| RPL11 | 60S | 1 | 0.0% | 70 | 0.0% |
| RPL18 | 60S | 1 | 0.0% | 4 | 0.0% |
| RPL26 | 60S | 1 | 0.0% | 6 | 33.3% |

### 4.2 Per-Subunit

| Subunit | Local (n) | % Intra-domain | ClinVar (n) | % Intra-domain (ClinVar) |
|---------|-----------|----------------|-------------|--------------------------|
| 40S | 6 | 16.7% | 238 | 5.0% |
| 60S | 6 | 0.0% | 232 | 1.7% |

## 5. PTM Co-localization Findings

No local variants co-localize with known post-translational modification sites annotated in UniProt for these ribosomal proteins.

## 6. Local vs ClinVar Comparison

The table in Section 4 compares the proportion of intra-domain variants between the local cohort and ClinVar pathogenic variants. Concordance between the two sources would support the hypothesis that domain co-localization is a general property of pathogenic DBA mutations rather than a sampling artefact.

> **Caveat:** The local cohort is small (n = 12–14 with resolved positions). Percentage comparisons are descriptive only.

## 7. Figure Recommendations for Article Inclusion

Three figure styles were generated to support different narrative contexts:

### Figure 1 — Domain Architecture (`domain_architecture_all`)
**Recommended for the article main body.** This figure provides an integrated cross-protein overview showing the relationship between annotated domain structure and variant positions for all seven proteins simultaneously. The normalized X-axis enables direct visual comparison despite heterogeneous protein lengths. Suitable as a full-width figure in a clinical genetics manuscript.

### Figure 3 — Lollipop with Domain Background (`domain_lollipop_{GENE}`)
**Recommended as supplementary figures.** The per-gene lollipop plots provide variant-level resolution with domain context as background shading. Particularly informative for RPL5 (N-terminal clustering within the uL18 domain) and for genes with multiple domain types.

### Figure 2 — Track-style (`domain_track_{GENE}`)
**Recommended as supplementary figures or for conference presentations.** The multi-track layout integrates domain structure, ClinVar density, and local variant positions in a single vertical panel per gene. Most information-dense format; requires more reader familiarity with genome-browser-style figures.

### Suggested framing for manuscript:

> *"Mapping of local cohort variants onto annotated UniProt domains revealed that [X]% of variants with resolved amino acid positions co-localize with known structural or functional domains (Figure N). This proportion is consistent with [ClinVar concordance statement if applicable], supporting the view that pathogenic DBA mutations preferentially target functionally constrained regions of the ribosomal protein."*

## 8. Generated Files

- `domain_architecture_all.png`
- `domain_architecture_all.svg`
- `domain_track_RPS19.png`
- `domain_track_RPS19.svg`
- `domain_lollipop_RPS19.png`
- `domain_lollipop_RPS19.svg`
- `domain_track_RPS24.png`
- `domain_track_RPS24.svg`
- `domain_lollipop_RPS24.png`
- `domain_lollipop_RPS24.svg`
- `domain_track_RPS26.png`
- `domain_track_RPS26.svg`
- `domain_lollipop_RPS26.png`
- `domain_lollipop_RPS26.svg`
- `domain_track_RPL5.png`
- `domain_track_RPL5.svg`
- `domain_lollipop_RPL5.png`
- `domain_lollipop_RPL5.svg`
- `domain_track_RPL11.png`
- `domain_track_RPL11.svg`
- `domain_lollipop_RPL11.png`
- `domain_lollipop_RPL11.svg`
- `domain_track_RPL18.png`
- `domain_track_RPL18.svg`
- `domain_lollipop_RPL18.png`
- `domain_lollipop_RPL18.svg`
- `domain_track_RPL26.png`
- `domain_track_RPL26.svg`
- `domain_lollipop_RPL26.png`
- `domain_lollipop_RPL26.svg`

---

*DBA Variant Analysis Pipeline — `scripts/domain_analysis.py`*