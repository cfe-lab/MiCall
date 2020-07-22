---
title: CFE HCV Algorithm
subtitle: version cfe-hcv 1.5 - 07-Dec-2016
---

Rules are applied for each drug and each genotype found in a sample. Each rule
contains a score and a list of the mutation descriptions that receive that
score. A mutation description is the position number of an amino acid within
the gene region, plus a variant amino acid. An exclamation mark in the mutation
means any variant except those listed. A "TRUE" mutation is always applied for
that genotype. A score is either a number, or a name in quotes. The program
adds up all the numerical scores and collects all the names, before deciding on
a resistance report.

* If the score is 4, then it reports "Resistance Possible".
* If the score is 8 or more, then it reports "Resistance Likely".
* If the score is 0 with the "Not available" name, then it reports "Resistance
    Interpretation Not Available", which means that resistance interpretation
    is not available for the specified genotype and drug.
* If the score is 0 with the "Not indicated" name, then it reports
    "Not Indicated", which means that the drug is not indicated for use in
    Canada.
* If the score is 0 with the "Effect unknown" name, then it reports "Mutations
    Detected; Effect Unknown".
* If there is insufficient coverage due to poor amplification, sequencing, or
    mapping, then it reports "Sequence does not meet quality-control standards".
* Otherwise, it reports "Likely Susceptible".

Here are the rules used for each of the drugs. Below that is a description of
how the scores are chosen.

| Gene | Drug | Genotype | Score | Mutations |
|------|------|----------|-------|-----------|
| NS3  | Boceprevir (Victrelis™)                                   | 1A | 8 | 36AM, 54AS, 55AI, 155K, 156ST, 170FTV |
|      |                                                           |    | 4 | 107I, 158I, 168N |
|      |                                                           |    | Effect unknown | R155!KR, V36!AMV, T54!AST, V55!AIV, D168!DN, I170!FITV, V107!IV, A156!AST, V158!IV |
|      |                                                           | 1B | 8 | 36AM, 54ACGS, 55A, 155CK, 156STV, 170AT |
|      |                                                           |    | 4 | 107I, 158I, 175L |
|      |                                                           |    | Effect unknown | R155!CKR, V36!AMV, T54!ACGST, V55!AV, V170!AITV, V107!IV, A156!ASTV, V158!IV, M175!LM |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Grazoprevir (a component of Zepatier™)                    | 1A | 8 | 56H, 156GLTV, 168!CDMNPQRSW |
|      |                                                           |    | 4 | 36ALM, 56F, 107I, 122GT, 155GIKS, 156M, 158A, 168CNS, 170V |
|      |                                                           |    | Effect unknown | D168MPQRW, R155!GIKRS, V36!ALMV, Y56!FHY, S122!GRST, V107!IV, A156!AGLMSTV, I170!ITV, V158!AV |
|      |                                                           | 1B | 8 | 56H, 155GTW, 156TV, 168AFGHIKLTV |
|      |                                                           |    | 4 | 56F, 107I, 156G, 168ENY |
|      |                                                           |    | Effect unknown | Y56!FHY, A156!AGTV, V107!IV, R155!GKRTW, D168CMPQRSW |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | 8 | 77S, 168QR, 178R |
|      |                                                           |    | Effect unknown | Q168!QR, Q178!QR, V170!V, Y56!Y, S122!S, R155!R, A156!A, N77!NS |
|      |                                                           | 4  | 8 | 168AV |
|      |                                                           |    | 4 | 107I, 122G, 156GMTV, 158I, 168CEGNY, 170I |
|      |                                                           |    | Effect unknown | D168!ACDEGNVY, R155!R, Y56!Y, T122!GT, V107!IV, A156!AGMTV, V170!IV, V158!IV |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Paritaprevir (a component of Technivie™ and Holkira Pak™) | 1A | 8 | 43L, 155GKTW, 156T, 168AEFHILNTVY |
|      |                                                           |    | 4 | 36AMT, 55I, 56H, 80KL, 132V, 155S, 156GSV, 334S, 342P, 357K, 406AI, 449I, 470S |
|      |                                                           |    | Effect unknown | T449!IT, V36!ALMTV, E357!EK, V406!AIV, D168CGKMPQRSW, P470!PS, F43!FL, P334!PS, Q80!KLQR, S342!PS, V55!IV, Y56!HY, I132!IV, R155!GKRSTW, A156!AGSTV |
|      |                                                           | 1B | 8 | 155K, 156V, 168AFHKTVY |
|      |                                                           |    | 4 | 56H, 156T, 168ILN, 357K |
|      |                                                           |    | Effect unknown | Y56!HY, D168CGMPQRSW, R155!KR, A156!ASTV, E357!EK |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | 8 | 56H, 155C, 156TV, 168HV |
|      |                                                           |    | Effect unknown | Y56!HY, D168!DHV, R155!CR, A156!ATV |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Simeprevir (Galexos™)                                     | 1A | 8 | 80KR, 122GR, 155KT, 156GTV, 168AEHV, 170T |
|      |                                                           |    | 4 | 36L, 43ILV, 54S, 80GHILM, 122AINT, 155GIMQSW, 156LMS, 168CFGIKLNSTY, 170V |
|      |                                                           |    | Effect unknown | Q80!GHIKLMQR, R155!GIKMQRSTW, V36!LV, T54!ST, D168MPQRW, S122!AGINRST, F43!FILV, A156!AGLMSTV, I170!ITV |
|      |                                                           | 1B | 8 | 41R, 43ISV, 80KR, 122AIRT, 155GKQTW, 156GTV, 168!CDKLMPRSW |
|      |                                                           |    | 4 | 43L, 80HIM, 155C, 156S, 168KL, 170T |
|      |                                                           |    | Effect unknown | Q80!GHIKLMQR, R155!CGIKMQRTW, D168CMPRSW, Q41!QR, S122!AIRST, F43!FILSV, A156!AGSTV, V170!ATV |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | 8 | 155K, 156G |
|      |                                                           |    | Effect unknown | R155!KR, A156!AGM |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Telaprevir (Incivek™)                                     | 1A | 8 | 36AGILM, 54AS, 155GKMT, 156STV |
|      |                                                           |    | 4 | 132V, 168N |
|      |                                                           |    | Effect unknown | I132!IV, T54!AST, D168!DN, V36!AGILMV, R155!GKMRT, A156!ASTV |
|      |                                                           | 1B | 8 | 36AGILM, 54AS, 155K, 156FSTV |
|      |                                                           |    | Effect unknown | A156!AFSTV, R155!KR, V36!AGILMV, T54!AST |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
| NS5a | Daclatasvir (Daklinza™)                                   | 1A | 8 | 28AT, 30DEGHKNRY, 31IMV, 58D, 93CHNS |
|      |                                                           |    | 4 | 24GR, 28SV, 30T, 54R, 58PR |
|      |                                                           |    | Effect unknown | H54!HR, Y93!CFHNSY, K24!GKNR, A92!AK, K26!EK, M28!AMSTV, H58!DHPR, Q30ACFIMPSVW, L31!ILMV |
|      |                                                           | 1B | 8 | 31FV, 32X, 93H |
|      |                                                           |    | 4 | 28MT, 29SX, 30GHPQ, 31IM, 32L, 58S, 62D, 92K, 93N |
|      |                                                           |    | Effect unknown | P32!LP, Q62!DQ, Y93!HNY, A92!AK, P58!PS, L28!LMT, P29!PS, R30!GHPQR, L31!FILMV |
|      |                                                           | 2  | 8 | 28CS, 31M, 92R, 93H |
|      |                                                           |    | Effect unknown | C92!CR, F28!CFS, Y93!HY, L31!LM |
|      |                                                           | 3  | 8 | 30K, 31FIMV, 93H |
|      |                                                           |    | 4 | 30S, 31P, 62AILPRT |
|      |                                                           |    | Effect unknown | Y93!HY, M28!MTV, S62!AILPRST, A30!AKSTV, L31!FILMPV |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Elbasvir (a component of Zepatier™)                       | 1A | 8 | 28AGT, 30DEGHR, 31FIMV, 58D, 93CHNS |
|      |                                                           |    | 4 | 30KY |
|      |                                                           |    | Effect unknown | H58!DH, M28!AGMTV, Y93!CHNSY, Q30!DEGHKLQRY, L31!FILMV |
|      |                                                           | 1B | 8 | 31FMV, 93H |
|      |                                                           |    | 4 | 28M |
|      |                                                           |    | Effect unknown | L28!LM, Y93!HY, L31!FLMV |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | 8 | 30DK, 31FM, 93H |
|      |                                                           |    | Effect unknown | Y93!HY, A30!ADK, L31!FLM |
|      |                                                           | 4  | 8 | 30FH, 31V, 58D, 93H |
|      |                                                           |    | 4 | 31I, 32L |
|      |                                                           |    | Effect unknown | P32!LP, P58!DP, Y93!HY, L30!FHLS, M31!IMV |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Ledipasvir (a component of Harvoni™)                      | 1A | 8 | 28AGT, 30EGHKNRY, 31IMV, 32L, 38F, 58D, 92P, 93CHNS |
|      |                                                           |    | 4 | 24GNR, 28V, 30LT, 31P, 58P, 92T, 93F |
|      |                                                           |    | Effect unknown | P32!LP, S38!FS, K24!GKNR, A92!APT, H58!DHP, M28!AGMTV, Y93!CFHNSY, Q30ACDFIMPSVW, L31!ILMPV |
|      |                                                           | 1B | 8 | 31IV, 58D, 92K, 93H |
|      |                                                           |    | 4 | 28M, 31M, 92T, 93C |
|      |                                                           |    | Effect unknown | P32!LP, L28!LM, P58!DP, A92!AKT, Y93!CHSY, L31!FILMV |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Ombitasvir (a component of Technivie™ and Holkira Pak™)   | 1A | 8 | 28TV, 30EKRY, 31V, 58D, 93CFHLNS |
|      |                                                           |    | 4 | 24R, 28A, 32L, 54Y, 58PR |
|      |                                                           |    | Effect unknown | P32!LP, H54!HY, K24!KR, H58!DHPR, M28!AMTV, Y93!CFHLNSY, Q30!EHKLNQRTY, L31!ILMV |
|      |                                                           | 1B | 8 | 28T, 29X, 31FV, 93HNS |
|      |                                                           |    | 4 | 28M, 30GHPQ, 31MT, 32X, 54Y, 58AS, 92E |
|      |                                                           |    | Effect unknown | P32!P, Q54!QY, Y93!HNSY, A92!AE, P58!APS, L28!LMT, P29!P, R30!GHPQR, L31!FILMTV |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | 8 | 28V |
|      |                                                           |    | 4 | 28S, 31I, 58S |
|      |                                                           |    | Effect unknown | P58!APS, L28!LSV, M31!ILM |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Velpatasvir (a component of Epclusa™)                     | 1A | 8 | 28GT, 30EK, 31MV, 32L, 58D, 93CHNRSW |
|      |                                                           |    | 4 | 24MRT, 28V, 30HLR, 31I, 58P, 92K |
|      |                                                           |    | Effect unknown | P32!LP, K24!KMRT, A92!AK, H58!DHP, M28!GMTV, Y93!CHNRSWY, Q30!EHKLQR, L31!ILMV |
|      |                                                           | 1B | 8 | 92K |
|      |                                                           |    | 4 | 31IMV, 93CHNS |
|      |                                                           |    | Effect unknown | Q24!KQ, P58!PRT, A92!AK, Y93!CHNSY, R30!RS, L31!ILMV |
|      |                                                           | 2  | 8 | 28S, 31V, 92T, 93HN |
|      |                                                           |    | 4 | 31IM |
|      |                                                           |    | Effect unknown | C92!CT, P58!APS, F28!FS, Y93!HNY, L31!ILMV |
|      |                                                           | 3  | 8 | 30K, 31M, 93HS |
|      |                                                           |    | 4 | 30V, 31PV, 38P, 58L, 92K, 93NR |
|      |                                                           |    | Effect unknown | S38!PS, P58!LP, E92!EK, Y93!HNRSY, A30!AKV, L31!LMPV |
|      |                                                           | 4  | 4 | 28T, 30HS, 31V, 32L, 93CHNSW |
|      |                                                           |    | Effect unknown | P32!LP, L28!LT, Y93!CHNSWY, L30!HLS, M31!MV |
|      |                                                           | 5  | Effect unknown | L28!L, P58!P, A92!A, T93!T, Q30!Q, L31!L |
|      |                                                           | 6  | 8 | 31V, 32ALQR |
|      |                                                           |    | Effect unknown | P32!ALPQR, A92!A, T58!T, F28!F, T93!T, R30!R, L31!LV |
|      |                                                           | 6E | 8 | TRUE |
| NS5b | Dasabuvir (a component of Holkira Pak™)                   | 1A | 8 | 314H, 316Y, 395G, 414ITV, 444K, 446KQ, 448CH, 553TV, 554S, 556GR, 561H |
|      |                                                           |    | 4 | 307R, 450V, 553I, 558R, 559GINV |
|      |                                                           |    | Effect unknown | Y448!CHY, A450!AV, A553!AITV, G554!GS, A395!AG, S556!GRS, G558!GR, D559!DGINV, Y561!HY, G307!GR, M414!IMTV, C316!CY, L314!HL, N444!DKN, E446!EKQ |
|      |                                                           | 1B | 8 | 316HNY, 368T, 411S, 414ITV, 448CH, 553V, 556G, 559G |
|      |                                                           |    | 4 | 316W, 445FY, 556R |
|      |                                                           |    | Effect unknown | S368!ST, Y448!CHY, K307!KR, D559!DG, A553!AV, N411!NS, C316!CHNWY, C445!CFY, M414!IMTV, S556!GRS |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not indicated | TRUE |
|      | Sofosbuvir (a component of Epclusa™)                      | 1A | 8 | 282T |
|      |                                                           |    | 4 | 61G, 112T, 159F, 237G, 320V, 321I, 321A, 355H, 473T |
|      |                                                           |    | Effect unknown | A112!AT, L320!LV, Q355!HQ, V321!AIV, S473!ST, S282!ST, L159!FL, C316!C, E237!EG, D61!DG |
|      |                                                           | 1B | 8 | 282T |
|      |                                                           |    | 4 | 142S, 321I |
|      |                                                           |    | Effect unknown | L320!L, V321!IV, S282!ST, C316!C, N142!NS, L159!L |
|      |                                                           | 2  | 8 | 282T |
|      |                                                           |    | Effect unknown | L320!L, V321!V, S282!ST, C316!C, L159!L |
|      |                                                           | 3  | 8 | 282T |
|      |                                                           |    | 4 | 142T, 159F, 237G, 314FIP, 321A, 355H |
|      |                                                           |    | Effect unknown | L320!L, V321!AV, Q355!HQ, S282!ST, L159!FL, C316!C, L314!FILP, N142!NT, E237!EG |
|      |                                                           | 4  | 8 | 282T |
|      |                                                           |    | 4 | 237G, 321I |
|      |                                                           |    | Effect unknown | L320!L, V321!IV, S282!ST, C316!C, E237!EG, L159!L |
|      |                                                           | 5  | 8 | 282T |
|      |                                                           |    | 4 | 289I |
|      |                                                           |    | Effect unknown | L320!L, M289!IM, V321!V, S282!ST, C316!C, L159!L |
|      |                                                           | 6  | 8 | 282T |
|      |                                                           |    | Effect unknown | L320!L, V321!V, S282!ST, C316!C, L159!L |
|      | Sofosbuvir (a component of Harvoni™)                      | 1A | 8 | 282T |
|      |                                                           |    | 4 | 61G, 112T, 159F, 237G, 320V, 321I, 321A, 355H, 473T |
|      |                                                           |    | Effect unknown | A112!AT, L320!LV, Q355!HQ, V321!AIV, S473!ST, S282!ST, L159!FL, C316!C, E237!EG, D61!DG |
|      |                                                           | 1B | 8 | 282T |
|      |                                                           |    | 4 | 142S, 321I |
|      |                                                           |    | Effect unknown | L320!L, V321!IV, S282!ST, C316!C, N142!NS, L159!L |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      |                                                           | 6E | Not available | TRUE |

The overall resistance scores of a mutant strain are derived from drug
susceptibility and clinical observations.

The drug susceptibility of a mutant is compared with the drug susceptibility of
the wild type, and the fold change is reported. For example, a mutant with five
times lower drug susceptibility than the wild type reports a fold change of 5.

Some mutants have a single amino acid changed from the wild type, and some
mutants have more than one change.

Another type of study reports mutations that are observed in clinical patients
with virological failure.

Based on all these study results, we assign the scores shown in the table
above. For each drug and genotype combination, we choose a middle range of fold
change values for drug susceptibility. Any mutations that fall in that middle
range get a score of 4 (resistance possible), mutations above that range get a
score of 8 (resistance likely), and mutations below that range get no score
(susceptible). Mutations observed in clinical patients with virological failure
get a score of 4 (resistance possible).

If reports for mutants with more than one change disagree with the reports for
the individual changes, then the combination takes precedence.

As an example, here are the fold change ranges for NS3 drugs, genotype 1a:

| Drugs | Middle Range |
|-------|--------------|
|Grazoprevir, Paritaprevir, Glecaprevir, Simeprevir, Asunaprevir| 5 - 10 |
|Voxilaprevir| 2.6 - 10* |

`*` Some studies only reported ranges of fold changes, so we used the 2.5 - 20
range as a middle range for those studies.
