# Longitudinal Program Impact & Heterogeneous Treatment Effect Modeling

![Econometric Diagnostics](econometric_diagnostics.png)

## Overview
This repository presents a quantitative statistical modeling framework designed to evaluate longitudinal intervention impacts and conditional treatment heterogeneity in observational panel data. Originally formulated as an empirical capstone study on multi-year student assessment cohorts ($N = 376$ analytical sample drawn from $1,100+$ student-year administrative records across Grades PK–5), this implementation reconstructs the full econometric workflow using reproducible Python pipelines and privacy-preserving synthetic panel data.

---

## Core Quantitative Methodology

### 1. Panel Cohort Engineering & Selection Diagnostics
* **Cohort Reconstruction:** Tracked longitudinal student records across three academic cycles (2011–2014) under incomplete administrative test coverage and non-random program assignment.
* **Attrition & Balance Testing:** Executed Pearson’s chi-square ($\chi^2$) diagnostics comparing the analytical cohort ($N = 376$) against excluded records ($N = 728$) to verify demographic representativeness ($p > 0.05$ across covariates), while quantifying baseline covariate imbalances between treatment ($n = 57$) and control ($n = 319$) groups.

### 2. Multi-Scale Psychometric Standardization
* Normalized heterogeneous state-level (PSSA) and local (DRP) assessment scores into within-cell $z$-scores stratified by assessment instrument, grade level, and academic year:
  $$z_{i,g,t} = \frac{X_{i,g,t} - \mu_{g,t}}{\sigma_{g,t}}$$
* Enabled scale-invariant longitudinal growth comparisons across distinct student cohorts.

### 3. Progressive OLS Specifications & Subgroup Interaction Modeling

- *Sequential Specification ($M_1$–$M_5$):* Estimated nested linear regression models controlling for baseline pretest scores and demographic/institutional covariates (`IEP`, `ELL`, `Gender`, `Ethnicity`):

$$
Y_{i} = \beta_{0} + \beta_{1} \text{Treatment}_{i} + \beta_{2} Y_{i,\text{pre}} + \mathbf{X}_{i}'\boldsymbol{\gamma} + \varepsilon_{i}
$$

  Incorporating full covariate controls increased explained variance ($R^2$) from $0.186$ to $0.273$ in local assessment models ($N = 131$) and from $0.204$ to $0.308$ in state assessment models ($N = 52$).
- *Heterogeneous Treatment Effects (CATE):* Introduced multiplicative interaction terms ($\text{Treatment}_{i} \times \text{Subgroup}_{i}$) across baseline achievement tiers and demographic indicators to test for differential subgroup gains.
- *Early-Grade Paired Inference:* Conducted paired-sample $t$-tests on early literacy cohorts (Grades K–1, $N = 160$) to isolate short-term skill acquisition.

### 4. Statistical Power & Minimum Detectable Effect (MDE) Diagnostics

- Quantified sample-size constraints in upper-grade observational subsamples by deriving the Minimum Detectable Effect (MDE) at $\alpha = 0.05$ (two-sided) and statistical power $1 - \beta = 0.80$:

$$
\text{MDE} = (z_{1-\alpha/2} + z_{1-\beta}) \sqrt{(1 - R^2)\left(\frac{1}{n_{T}} + \frac{1}{n_{C}}\right)}
$$

- Demonstrated mathematically why small analytical subsamples require large effect sizes to achieve statistical significance, motivating covariate-matching (PSM/IPW) and randomized experimental designs (RCT).

---

## Repository Structure
```text
├── README.md                              # Technical methodology & econometric documentation
├── longitudinal_impact_modeling.ipynb     # Interactive Jupyter Notebook (data simulation, OLS, MDE & plots)
├── synthetic_panel_data.csv               # Privacy-preserving synthetic panel dataset (N = 376)
└── econometric_diagnostics.png            # Coefficient forest plot & statistical power MDE sensitivity curve
