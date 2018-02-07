---
title: CFE HCV Algorithm
description: version cfe-hcv 1.5 - 07-Dec-2016
---

Rules are applied for each drug and each genotype found in a sample. Each rule
contains a mutation description, the `=>` symbol, and an action. A mutation
description is the position number of an amino acid within the gene region,
plus a variant amino acid. An exclamation mark in the mutation means any
variant except those listed. A "TRUE" mutation is always applied for that
genotype. An action is either a number that adds to the score, or a name in
quotes. The program adds up all the scores and collects all the names, before
deciding on a resistance report.

* If the score is 4, then it reports "Resistance Possible".
* If the score is 8 or more, then it reports "Resistance Likely".
* If the score is 0 with the "Not available" name, then it reports "Resistance
    Interpretation Not Available".
* If the score is 0 with the "Not indicated" name, then it reports
    "Not Indicated".
* If the score is 0 with the "Effect unknown" name, then it reports "Mutations
    Detected; Effect Unknown".
* If there are gaps in coverage for some of the mutation locations, then it
    reports "Sequence does not meet quality-control standards".
* Otherwise, it reports "Likely Susceptible".

Here are the rules used for each of the drugs.

| Gene | Drug | Genotype | Rules |
|------|------|----------|-------|
| NS3  | Boceprevir (Victrelis™)                                   | 1A | SCORE FROM |