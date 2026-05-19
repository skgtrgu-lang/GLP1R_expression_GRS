# GLP1R_expression_GRS
Code used to generate GLP1R expression GRS scores using PLINK2

# GLP1R Expression Genetic Risk Score (GRS)
**Exenatide-PD3 Pharmacogenomics | Rowan Gurney, UCL Queen Square Institute of Neurology**
---
## Overview

This repository contains the score file and analysis scripts used to generate a
**GLP1R expression genetic risk score (GRS)** in participants from the
Exenatide-PD3 randomised controlled trial.

The GRS is a weighted sum of 15 cis-eQTL variants associated with GLP1R gene
expression in whole blood (GTEx v8, European ancestry), as reported by
Triozzi et al. A higher score indicates higher genetically predicted GLP1R
expression.

This score was used as a pharmacogenetic stratification instrument in:
> Gurney R et al. Genetically proxied GLP1R expression does not predict exenatide response in the Exenatide-PD3 trial. Brain, 2026.
---
## Score Description
| Property | Value |
|---|---|
| Gene | GLP1R (chromosome 6) |
| Variants | 15 cis-eQTL variants |
| Weights | eQTL effect sizes (betas) from GTEx v8 whole blood EUR |
| Source | Triozzi et al., Supplementary Table S3 (eTable 2) |
| Genome build | hg38 (GRCh38) |
| Score type | Weighted sum: sum(dosage * beta) |
| Direction | Higher score = higher predicted GLP1R expression |

---

## Repository Contents

```
