---
title: CFE HCV Algorithm
description: version cfe-hcv 1.6 - 23-Apr-2018
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
| NS3  | Asunaprevir (a component of Sunvepra™)                    | 1A | 8 | R155K, D168AEGVY |
|      |                                                           |    | 4 | V36GLM, T54S, V55A, Q80KLR, S122GINT, R155GQ, D168FHNT, I170T |
|      |                                                           |    | Effect unknown | V36!GLMV, Q80!KLQR, S122!GINST, R155!GKQR, A156!A, D168CIKLMPQRSW |
|      |                                                           | 1B | 8 | Q80K, R155K, A156V, D168ACEFGHVY |
|      |                                                           |    | 4 | V36GM, T54S, Y56HL, N77S, Q80LR, S122DGINT, R155GQ, A156ST, D168NT, V170A |
|      |                                                           |    | Effect unknown | V36!GMV, Y56!HLY, Q80!KLQR, S122!DGINST, R155!GKQR, A156!ASTV, D168IKLMPQRSW |
|      |                                                           | 2  | 8 | TRUE |
|      |                                                           | 3  | 8 | TRUE |
|      | Boceprevir (Victrelis™)                                   | 1A | 8 | V36AM, T54AS, V55AI, V107I, R155K, A156ST, D168N, I170FTV |
|      |                                                           |    | 4 | V158I |
|      |                                                           |    | Effect unknown | V36!AMV, T54!AST, V55!AIV, V107!IV, R155!KR, A156!AST, V158!IV, D168!DN, I170!FITV |
|      |                                                           | 1B | 8 | V36AM, T54ACGS, V55A, V107I, R155CK, A156STV, V158I, V170AT, M175L |
|      |                                                           |    | Effect unknown | V36!AMV, T54!ACGST, V55!AV, V107!IV, R155!CKR, A156!ASTV, V158!IV, V170!ATV, M175!LM |
|      | Glecaprevir (a component of Mavyret™)                     | 1A | 8 | A156T, D168AFVY |
|      |                                                           |    | 4 | V36AFM, V55AI, Y56H, Q80KL, I132V, R155KT, A156GV, D168EHT |
|      |                                                           |    | Effect unknown | V36!AFMV, V55!AIV, Q80!KLQ, R155!KRT, A156!AGTV, D168!ADEFHTVY |
|      |                                                           | 1B | 8 | A156TV, D168K |
|      |                                                           |    | 4 | Y56F, D168FT, E176G |
|      |                                                           |    | Effect unknown | A156!ATV, D168!DFKT |
|      |                                                           | 2  | 8 | A156MTV |
|      |                                                           |    | 4 | I132L, D168EFHT |
|      |                                                           |    | Effect unknown | A156!AMTV, D168!DEFHT |
|      |                                                           | 3  | 8 | Q80KR, A156G, Q168LR |
|      |                                                           |    | 4 | T54S, Y56HN, V107I, L132I, A166STY, Q168K, D186LR, I366V |
|      |                                                           |    | Effect unknown | Y56!HNY, Q80!KQR, A156!AG, A166!ASTY, Q168!KLQR, D186!DLR |
|      |                                                           | 4  | 8 | A156TV, D168H |
|      |                                                           |    | 4 | Y56H, G90R, D168V |
|      |                                                           |    | Effect unknown | A156!ATV, D168!DHV |
|      |                                                           | 5  | 4 | D168E |
|      |                                                           |    | Effect unknown | A156!A, D168!DE |
|      |                                                           | 6  | 8 | L80R, D168AGHVY |
|      |                                                           |    | Effect unknown | L80!LR, A156!A, D168!ADGHVY |
|      | Grazoprevir (a component of Zepatier™)                    | 1A | 8 | Y56H, Q80KL, R155KQTW, A156GLTV, D168!CDMPQRSW |
|      |                                                           |    | 4 | V36ALM, Q41R, I48V, V55A, Y56F, Q80S, V107I, S122GT, I132V, R155GIS, V158A, D168C, I170V, T185S, E357GK |
|      |                                                           |    | Effect unknown | V36!ALMV, Y56!FHY, Q80!KLQS, S122!GST, R155!GIKQRSTW, A156!AGLTV, D168MPQRSW |
|      |                                                           | 1B | 8 | Y56H, R155GTW, A156TV, D168AEFGHIKLTV |
|      |                                                           |    | 4 | T54S, Y56F, Q80L, V107I, S122NRT, V132IL, R155Q, D168NY, V170I, T185S, E357K |
|      |                                                           |    | Effect unknown | Y56!FHY, S122!NRST, V132!ILV, R155!GQRTW, A156!ATV, D168CMPQRSW |
|      |                                                           | 3  | 8 | A156TV |
|      |                                                           |    | 4 | Y56H, Q178R |
|      |                                                           |    | Effect unknown | A156!ATV, Q168!Q |
|      |                                                           | 4  | 8 | D168AV |
|      |                                                           |    | 4 | V107I, A156MTV, D168G, V170I |
|      |                                                           |    | Effect unknown | A156!AMTV, D168!ADGV |
|      |                                                           | 6  | 4 | V36IL, Y56H, L80KQ, I132L, A156MTV, D168CEY, I170V |
|      |                                                           |    | Effect unknown | V36!ILV, L80!KLQ, A156!AMTV, D168!CDEY |
|      | Paritaprevir (a component of Technivie™ and Holkira Pak™) | 1A | 8 | F43L, R155GKTW, A156T, D168!CDGMPQRSW |
|      |                                                           |    | 4 | V36AMT, V55I, Y56H, Q80KL, I132V, R155SV, A156G, P334S, S342P, E357K, V406AI, T449I, P470S |
|      |                                                           |    | Effect unknown | V36!AMTV, F43!FL, Q80!KLQ, R155!GKRSTVW, A156!AGT, D168CGMPQRSW, V406!AIV |
|      |                                                           | 1B | 8 | R155K, A156V, D168AFHKTVY |
|      |                                                           |    | 4 | Y56FH, A156T, D168ILN, E357K |
|      |                                                           |    | Effect unknown | R155!KR, A156!ATV, D168CEGMPQRSW |
|      |                                                           | 3  | 8 | TRUE, Q168V |
|      |                                                           |    | 4 | Y56H |
|      |                                                           | 4  | 8 | Y56H, R155C, A156TV, D168HV |
|      |                                                           |    | Effect unknown | Y56!HY, R155!CR, A156!ATV, D168!DHV |
|      | Simeprevir (Galexos™)                                     | 1A | 8 | Q80KR, S122NR, R155KT, A156GTV, D168AEHV, I170T |
|      |                                                           |    | 4 | A5V, V36L, Q41R, T54S, Q80L, G120S, S122ACGIT, R123M, I132L, R155GQ, D168FT, I170V, N174RS, T178I, P264S, T343N, T344I, L356F, A477T, V511I, G534S, E555D |
|      |                                                           |    | Effect unknown | Q80!KLQR, S122!ACGINRST, R155!GKQRT, A156!AGTV, D168!ADEFHTV, I170!ITV |
|      |                                                           | 1B | 8 | F43ISV, Q80KR, S122R, R155GKQTW, A156GTV, D168AEHINQTVY |
|      |                                                           |    | 4 | A39T, T46S, V48I, S61L, I71I, Q80H, Q86R, S122AGIT, D168F, V170IT, S174F, T344A, T358F, G383S, I426V, T448A, I472V, V535A, P574AS |
|      |                                                           |    | Effect unknown | F43!FISV, Q80!HKQR, S122!AGIRST, R155!GKQRTW, A156!AGTV, D168CGKLMPRSW, V170!ITV |
|      |                                                           | 3  | 4 | Q168V |
|      |                                                           | 4  | 8 | R155K, A156G |
|      |                                                           |    | Effect unknown | Q80!Q, T122!T, R155!KR, A156!AG, D168!D, V170!V |
|      | Telaprevir (Incivek™)                                     | 1A | 8 | V36AGILM, T54AS, R155GKMT, A156STV |
|      |                                                           |    | 4 | I132V, D168N |
|      |                                                           |    | Effect unknown | V36!AGILMV, T54!AST, I132!IV, R155!GKMRT, A156!ASTV, D168!DN |
|      |                                                           | 1B | 8 | V36AGILM, T54AS, R155K, A156FSTV |
|      |                                                           |    | Effect unknown | V36!AGILMV, T54!AST, R155!KR, A156!AFSTV |
|      | Voxilaprevir (a component of Vosevi™)                     | 1A | 8 | R155W, A156LT, D168KLR |
|      |                                                           |    | 4 | V36AG, Q41R, F43S, Q80K, R155G, A156IMPV, D168AFITV |
|      |                                                           |    | Effect unknown | V36!AGV, R155!GRW, A156!AILMPTV, D168!ADFIKLRTV |
|      |                                                           | 1B | 8 | R155W, A156TV |
|      |                                                           |    | 4 | V36AM, S122D, A156S, D168VY, V170A |
|      |                                                           |    | Effect unknown | V36!AMV, R155!RW, A156!ASTV, D168!DVY |
|      |                                                           | 2  | 8 | A156LTV |
|      |                                                           |    | 4 | F43V |
|      |                                                           |    | Effect unknown | A156!ALTV |
|      |                                                           | 3  | 8 | Q41K, Q80K, A156TV, Q168R |
|      |                                                           |    | 4 | Q41R, F43Y, V55A, L175M, A294V |
|      |                                                           |    | Effect unknown | Q41!KQR, Q80!KQ, A156!ATV, Q168!QR |
|      |                                                           | 4  | 8 | A156LTV |
|      |                                                           |    | 4 | Q41R, D168ETV, T185S, T382S |
|      |                                                           |    | Effect unknown | A156!ALTV, D168!DETV |
|      |                                                           | 5  | 8 | D168HR |
|      |                                                           |    | 4 | D168AKY |
|      |                                                           |    | Effect unknown | A156!A, D168!ADHKRY |
|      |                                                           | 6  | 4 | Q41KR, Y56H, D168AEH |
|      |                                                           |    | Effect unknown | Q41!KQR, A156!A, D168!ADEH |
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
|      | Velpatasvir (a component of Epclusa™ and Vosevi™)         | 1A | 8 | 28GT, 30EK, 31MV, 32L, 58D, 93CHNRSW |
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
|      | Sofosbuvir (a component of Epclusa™ and Vosevi™)          | 1A | 8 | 282T |
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
