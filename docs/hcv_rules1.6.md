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
| NS5a | Daclatasvir (Daklinza™)                                   | 1A | 8 | K24T, M28AGT, Q30DEGHKNRY, L31IMV, P32L, H58D, A92K, Y93CHNST |
|      |                                                           |    | 4 | K24R, M28V, Q30T, H54RY, H58P, A92P, Y93FR |
|      |                                                           |    | Effect unknown | M28!AGMTV, Q30ACFILMPSVW, L31!ILMV, P32!LP, H54!HRY, H58!DHP, A92!AKP, Y93!CFHNRSTY |
|      |                                                           | 1B | 8 | P29d, L31V, P32d, A92K, Y93H |
|      |                                                           |    | 4 | L28M, R30GHQ, L31FIMW, Q54H, P58S, Q62E, A92E, Y93N |
|      |                                                           |    | Effect unknown | L28!LM, R30!GHQR, L31!FILMVW, P32!P, A92!AEK, Y93!HNY |
|      |                                                           | 2  | 8 | TRUE, F28CS, L31M, C92R, Y93H |
|      |                                                           |    | Effect unknown | F28!CFS, L31!LM, C92!CR, Y93!HY |
|      |                                                           | 3  | 8 | M28K, A30EK, L31FIMV, Y93H |
|      |                                                           |    | 4 | M28T, A30SV, L31P, S62AILPRT |
|      |                                                           |    | Effect unknown | M28!KMT, A30!AEKSV, L31!FILMPV, S62!AILPRST, Y93!HY |
|      |                                                           | 4  | 8 | L28M, L30GHRS, Y93HRSW |
|      |                                                           |    | 4 | L30CPQY, M31V |
|      |                                                           | 5  | 8 | L31FV |
|      |                                                           | 6  | 8 | Q24H, L31M, P32LS, T58ANS |
|      | Elbasvir (a component of Zepatier™)                       | 1A | 8 | M28AGT, Q30DEGHKR, L31FIMV, H58D, Y93CHN |
|      |                                                           |    | 4 | K24R, M28V, Q30LY, E62D, Y93S, D427N |
|      |                                                           |    | Effect unknown | M28!AGMTV, Q30!DEGHKLQRY, L31!FILMV, H58!DH, Y93!CHNSY |
|      |                                                           | 1B | 8 | L31FMV, Y93H |
|      |                                                           |    | 4 | L28M, R30Q, P58S |
|      |                                                           |    | Effect unknown | L31!FLMV, Y93!HY |
|      |                                                           | 2  | 4 | T24AS, F28CL, L31IM, P58S |
|      |                                                           | 3  | 8 | A30DK, L31FM, Y93H |
|      |                                                           |    | 4 | P58S, S62APT |
|      |                                                           |    | Effect unknown | A30!ADK, L31!FLM, P58!PS, S62!APST, Y93!HY |
|      |                                                           | 4  | 8 | L30FH, P58D, Y93H |
|      |                                                           |    | 4 | L28MST, L30RS, M31IV, P32L, H54Y, P58TY, Y93C |
|      |                                                           |    | Effect unknown | L30!FHLRS, M31!IMV, P58!DPTY, Y93!CHY |
|      |                                                           | 6  | 4 | F28LM, L31M |
|      | Ledipasvir (a component of Harvoni™)                      | 1A | 8 | M28AGT, Q30EGHKNRY, L31IMV, P32L, S38F, H58D, A92K, Y93CHLNRSTW |
|      |                                                           |    | 4 | K24R, M28V, Q30LT, L31P, H58PR, Y93F |
|      |                                                           |    | Effect unknown | M28!AGMTV, Q30ACDFIMPSVW, L31!ILMPV, P32!LP, S38!FS, H58!DHPR, A92!AK, Y93ADEGIKMPQV |
|      |                                                           | 1B | 8 | P58D, A92K, Y93HNS |
|      |                                                           |    | 4 | L28M, L31IV, A92T, Y93C |
|      |                                                           |    | Effect unknown | L31!ILV, P58!DP, A92!AKT, Y93!CHNSY |
|      |                                                           | 2  | 8 | TRUE, L31M |
|      |                                                           | 3  | 8 | TRUE, A30KSV, Y93H |
|      |                                                           | 4  | 8 | L28M, L30R, Y93HSW |
|      |                                                           |    | 4 | L30S, M31V, Y93C |
|      |                                                           | 5  | 4 | L31M |
|      |                                                           | 6  | 8 | TRUE |
|      |                                                           |    | 4 | Q24K, F28V, R30A, T58P |
|      | Ombitasvir (a component of Technivie™ and Holkira Pak™)   | 1A | 8 | M28TV, Q30EKRY, L31V, H58D, Y93CFHLNS |
|      |                                                           |    | 4 | K24R, M28A, P32L, H54Y, H58PR, E62D |
|      |                                                           |    | Effect unknown | M28!AMTV, Q30!EKQRY, L31!LV, H58!DHPR, Y93!CFHLNSY |
|      |                                                           | 1B | 8 | L28T, P29d, Y93HN |
|      |                                                           |    | 4 | L28M, R30Q, L31FMV, Q54Y, P58S, Y93S |
|      |                                                           |    | Effect unknown | L28!LMT, P29!P, L31!FLMV, Y93!HNSY |
|      |                                                           | 2  | 8 | L31MV |
|      |                                                           |    | 4 | T24A, F28F, Y93H |
|      |                                                           | 3  | 8 | M28T, Y93H |
|      |                                                           |    | 4 | L31F |
|      |                                                           | 4  | 4 | L28S, M31I, P58S |
|      |                                                           |    | Effect unknown | L28!LS |
|      |                                                           | 5  | 8 | L28I, L31FV |
|      |                                                           | 6  | 8 | L31V, T58N |
|      |                                                           |    | 4 | T58AS |
|      | Pibrentasvir (a component of Mavyret™)                    | 1A | 8 | M28G, Q30Dd, E62A, Y93N |
|      |                                                           |    | 4 | K24EQR, M28ATV, P29QR, Q30EGHKLR, L31M, H58CDY, Y93H |
|      |                                                           |    | Effect unknown | K24!EKQR, M28!AGMTV, P29!PQR, Q30!DEGHKLQR, H58!CDHY, E62!AE, Y93!HNY |
|      |                                                           | 1B | 8 | P32d |
|      |                                                           |    | 4 | L28M, P29d, R30Q, L31K, P32KNRT |
|      |                                                           |    | Effect unknown | P32!KNPRT |
|      |                                                           | 2  | 4 | L31M |
|      |                                                           | 3  | 4 | S24FG, M28GK, P29Q, A30GK, L31FIM, Y93H |
|      |                                                           |    | Effect unknown | S24!FGS, M28!GKM, A30!AGK, L31!FILM, E92!E, Y93!HY |
|      | Velpatasvir (a component of Epclusa™ and Vosevi™)         | 1A | 8 | M28GT, Q30EK, L31IMV, P32L, H58D, A92K, Y93CHNRSW |
|      |                                                           |    | 4 | K24MT, M28AV, Q30GHLRT, L31F, S38Y, H58NPR, S85N, Y93LT, L153V, A213A, T377T |
|      |                                                           |    | Effect unknown | K24!KMT, M28!AGMTV, Q30!EGHKLQRT, L31!FILMV, P32!LP, H58!DHNPR, A92!AK, Y93!CHLNRSTWY |
|      |                                                           | 1B | 8 | A92K, Y93HN |
|      |                                                           |    | 4 | Q24K, L31IMV, P58T, Y93CRST |
|      |                                                           |    | Effect unknown | L31!ILMV, A92!AK, Y93!CHNRSTY |
|      |                                                           | 2  | 8 | F28S, Y93HN |
|      |                                                           |    | 4 | F28CF, L31IMV, C92R |
|      |                                                           |    | Effect unknown | F28!CFS, L31!ILMV, C92!CR, Y93!HNY |
|      |                                                           | 3  | 8 | M28T, A30K, L31FM, Y93HS |
|      |                                                           |    | 4 | S14E, A30HRV, L31PV, S38PY, K41K, V52M, P58GLT, H85Y, E92K, Y93NR, S103DP, A291P, S379P, S385Q, R404K, P407S |
|      |                                                           |    | Effect unknown | M28!MT, A30!AHKRV, L31!FLMPV, S38!PSY, P58!GLPT, Y93!HNRSY |
|      |                                                           | 4  | 8 | L28T, Y93S |
|      |                                                           |    | 4 | L30R, M31V, P32L, H54R, Y93HNW |
|      |                                                           |    | Effect unknown | L28!LT, Y93!HNSWY |
|      |                                                           | 5  | 4 | L31I |
|      |                                                           | 6  | 8 | TRUE, L31V, P32ALQR |
|      |                                                           |    | 4 | F28MV, L31IM, T58GH, A92T, T93AHNS |
|      |                                                           |    | Effect unknown | F28!FMV, L31!ILMV, P32!ALPQR, T58!GHT, T93!AHNST |
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
