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
| NS3  | Boceprevir (Victrelis™)                                   | 1A | SCORE FROM ( 36AM=>8, 54AS=>8, 55AI=>8, 107I=>4, 155K=>8, 156ST=>8, 158I=>4, 168N=>4, 170FTV=>8, R155!KR=>"Effect unknown", V36!AMV=>"Effect unknown", T54!AST=>"Effect unknown", V55!AIV=>"Effect unknown", D168!ND=>"Effect unknown", I170!FTVI=>"Effect unknown", V107!IV=>"Effect unknown", A156!STA=>"Effect unknown", V158!IV=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 36MA=>8, 54CGSA=>8, 55A=>8, 107I=>4, 155CK=>8, 156TVS=>8, 158I=>4, 170TA=>8, 175L=>4, 170I=>0, R155!CKR=>"Effect unknown", V36!MAV=>"Effect unknown", T54!CGSAT=>"Effect unknown", V55!AV=>"Effect unknown", V170!TAIV=>"Effect unknown", V107!IV=>"Effect unknown", A156!TVSA=>"Effect unknown", V158!IV=>"Effect unknown", M175!LM=>"Effect unknown" ) |
