# Stat-Microbiome-Metabolites
Statistical analysis of gut microbiota and metabolite profiles

Gut Microbiome and Metabolite Profiles Analysis
│
├── 1. Study Design & Metadata
│   ├── Sample metadata
│   ├── Clinical / anthropometric variables
│   ├── Experimental groups
│   └── Statistical test selection
│       ├── Normality assessment
│       └── Homogeneity of variance (Levene's test)
│
├── 2. Microbiome Data Processing
│   ├── Quality control
│   ├── Taxonomic assignment
│   ├── Phyloseq object
│   └── Filtering / normalization
│
├── 3. Microbiome Diversity
│   │
│   ├── Alpha Diversity
│   │   ├── Observed
│   │   ├── Chao1
│   │   └── Shannon
│   │
│   └── Beta Diversity
│       ├── Bray-Curtis
│       ├── Unweighted UniFrac
│       ├── Weighted UniFrac
│       └── PERMANOVA / PERMDISP
│
├── 4. Taxonomic Composition
│   ├── Phylum
│   ├── Class
│   ├── Order
│   ├── Family
│   ├── Genus
│   
│
├── 5. Differential Abundance
│   ├── Kruskal-Wallis test
│   ├── Dunn's post-hoc test
│   └── BH correction
│
├── 6. Metabolomics
│   ├── Data preprocessing
│   ├── Compound mapping
│   │   └── PubChem
│   ├── Data transformation
│   └── Differential metabolite analysis
│       ├── Normality assessment
│       ├── Homogeneity of variance
│       └── Appropriate statistical test
│
├── 7. Functional Profiling
│   ├── PICRUSt2
│   ├── KEGG Orthologs (KOs)
│   ├── Enzyme Commission (EC) numbers
│   ├── KO / EC mapping
│   │   └── KEGG database
│   ├── Functional annotation
│   │   └── KEGG BRITE
│   └── Differential analysis
│       └── ALDEx2
│           ├── Pathways
│           ├── KOs
│           └── ECs
│
├── 8. Microbiome–Metabolite Integration
│   ├── Multiple Factor Analysis (MFA)
│   ├── Redundancy Analysis (RDA)
│   
│
├── 9. Diet Analysis (CM group)
│   ├── Statistical analysis
│   │   ├── Kruskal-Wallis test
│   │   ├── Dunn's post-hoc test
│   │   └── BH correction
│   │
│   └── Diet–Omics Associations
│       └── HAllA
│           ├── Diet vs Genus
│           └── Diet vs Metabolites
│
└── 10. Integrated Functional Analysis
    ├── Functional annotation
    │   └── KEGG BRITE
    │
    ├── Differential functional analysis
    │   └── ALDEx2
    │
    └── HAllA Association Analysis
        ├── Genus vs Function
        ├── Metabolites vs Function
        └── Metabolites vs Genus




# Analysis Overview
- DA = Differential Analysis.

                    STUDY DATA
                        │
          ┌─────────────┼─────────────┐
          │             │             │
          ▼             ▼             ▼
      Metadata      Microbiome    Metabolomics
          │             │             │
          │             ▼             │
          │         Diversity         │
          │             │             │
          │             ▼             ▼
          │       Taxonomic DA     Metabolite DA
          │             │             │
          │             └──────┬──────┘
          │                    │
          │                    ▼
          │             Multi-omics
          │              Integration
          │             ┌──────┴──────┐
          │             │             │
          │            MFA           RDA
          │             │             │
          └─────────────┼─────────────┘
                        │
                        ▼
                Functional Profiling
                        │
             ┌──────────┴──────────┐
             │                     │
             ▼                     ▼
           ALDEx2                 HAllA
             │                     │
             ▼                     ▼
      Differential             Association
       Functions                Networks





# Statistical Analysis Framework
                   Dataset
                      │
                      ▼
             Check data structure
                      │
                      ▼
              Normality assessment
                      │
             ┌────────┴────────┐
             │                 │
          Normal          Non-normal
             │                 │
             ▼                 ▼
       Parametric       Non-parametric
          tests              tests
             │                 │
             └────────┬────────┘
                      ▼
             Multiple-testing
                correction
                      │
                      ▼
                 Interpretation




# Statistical Analysis Approach  

├── 1. Define Study Design
│   │
│   ├── 2 Groups
│   │
│   └── ≥3 Groups
│
├── 2. Assess Distributional Assumptions
│   │
│   ├── Normality
│   │   └── Shapiro-Wilk test
│   │
│   └── Homogeneity of Variance
│       └── Levene's test
│
├── 3. Two-Group Comparison
│   │
│   ├── Approximately Normal
│   │   │
│   │   ├── Homogeneous variance
│   │   │   └── Student's t-test
│   │   │
│   │   └── Heterogeneous variance
│   │       └── Welch's t-test
│   │
│   └── Non-normal
│       └── Wilcoxon rank-sum test
│           (Mann-Whitney U test)
│
├── 4. ≥3-Group Comparison
│   │
│   ├── Approximately Normal
│   │   │
│   │   ├── Homogeneous variance
│   │   │   └── One-way ANOVA
│   │   │       │
│   │   │       └── Significant?
│   │   │           └── Tukey HSD post-hoc test
│   │   │
│   │   └── Heterogeneous variance
│   │       └── Welch's ANOVA
│   │           │
│   │           └── Significant?
│   │               └── Games-Howell post-hoc test
│   │
│   └── Non-normal
│       └── Kruskal-Wallis test
│           │
│           └── Significant?
│               └── Dunn's post-hoc test
│
└── 5. Multiple-Testing Correction
    │
    ├── Single outcome / limited comparisons
    │   └── Appropriate post-hoc correction
    │
    └── High-dimensional features
        └── Benjamini-Hochberg (BH/FDR)






# Statistical Analysis for Cross-Sectional / Between-Group and Paired / Longitudinal Studies
│
├── 1. Study Design
│   │
│   ├── Independent / Between-Group
│   │   ├── 2 Groups
│   │   └── ≥3 Groups
│   │
│   ├── Paired / Longitudinal
│   │   ├── 2 Time Points
│   │   └── ≥3 Time Points
│   │
│   └── Multivariate / Community Data
│       └── Microbiome β-diversity
│
├── 2. Independent Between-Group Analysis
│   │
│   ├── 2 Groups
│   │   ├── Normal + Homogeneous
│   │   │   └── Student's t-test
│   │   ├── Normal + Heterogeneous
│   │   │   └── Welch's t-test
│   │   └── Non-normal
│   │       └── Wilcoxon rank-sum test
│   │
│   └── ≥3 Groups
│       ├── Normal + Homogeneous
│       │   └── One-way ANOVA
│       │       └── Tukey HSD
│       ├── Normal + Heterogeneous
│       │   └── Welch's ANOVA
│       │       └── Games-Howell
│       └── Non-normal
│           └── Kruskal-Wallis
│               └── Dunn's test
│
├── 3. Paired / Longitudinal Analysis
│   │
│   ├── 2 Time Points
│   │   ├── Normal
│   │   │   └── Paired t-test
│   │   └── Non-normal
│   │       └── Wilcoxon signed-rank test
│   │
│   └── ≥3 Time Points
│       ├── Parametric
│       │   └── Repeated-measures ANOVA
│       └── Non-parametric
│           └── Friedman test
│
├── 4. Microbiome Community Analysis
│   │
│   ├── β-diversity
│   │   ├── Bray-Curtis
│   │   ├── Unweighted UniFrac
│   │   └── Weighted UniFrac
│   │
│   ├── PERMANOVA
│   └── PERMDISP
│
└── 5. Multiple-Testing Correction
    ├── Post-hoc comparisons
    └── High-dimensional features
        └── BH / FDR




                Normality
                    ↓
       Shapiro-Wilk in EACH group
                    ↓
       ┌────────────┴────────────┐
       ↓                         ↓
 all groups normal          ≥1 group non-normal
       ↓                         ↓
   Levene's test            Kruskal-Wallis
       ↓
 ┌─────┴─────┐
 ↓           ↓
equal      unequal
variance   variance
 ↓           ↓
ANOVA     Welch ANOVA
             ↓
       Games-Howell
