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
|      |                                                           |    | Effect unknown | V36!AGLMV, Q80!KLQR, S122!GINRST, R155!GKQR, A156!A, D168CIKLMPQRSW |
|      |                                                           | 1B | 8 | Q80K, R155K, A156V, D168ACEGHVY |
|      |                                                           |    | 4 | V36GM, T54S, Y56HL, N77S, Q80LR, S122DGINT, R155GQ, A156ST, D168FNT, V170A |
|      |                                                           |    | Effect unknown | V36!AGLMV, Y56!HLY, Q80!KLQR, S122!DGINST, R155!GKQR, A156!ASTV, D168IKLMPQRSW |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Glecaprevir (a component of Mavyret™)                     | 1A | 8 | A156T, D168AFY |
|      |                                                           |    | 4 | V36AFM, V55AI, Y56H, Q80KL, I132V, R155KT, A156GV, D168EHTV |
|      |                                                           |    | Effect unknown | V36!AFLMV, V55!AIV, Q80!KLQR, R155!KMRSTV, A156!AGTV, D168!ADEFHNTVY |
|      |                                                           | 1B | 8 | A156TV, D168K |
|      |                                                           |    | 4 | Y56F, D168FTV, E176G |
|      |                                                           |    | Effect unknown | A156!ASTV, D168!ADEFHKTVY |
|      |                                                           | 2  | 8 | A156MTV |
|      |                                                           |    | 4 | I132L, D168EFHT |
|      |                                                           |    | Effect unknown | A156!AMTV, D168!ADEFHSTVY |
|      |                                                           | 3  | 8 | Q80KR, A156G, A166T, Q168LR |
|      |                                                           |    | 4 | T54S, Y56HN, V107I, L132I, A166SY, Q168K, D186LR, I366V |
|      |                                                           |    | Effect unknown | Y56!HNY, Q80!KQR, A156!AG, A166!ASTY, Q168!HKLQR, D186!DLR |
|      |                                                           | 4  | 8 | A156TV, D168H |
|      |                                                           |    | 4 | Y56H, G90R, R155C, D168V |
|      |                                                           |    | Effect unknown | A156!ATV, D168!DHV |
|      |                                                           | 5  | 4 | D168E |
|      |                                                           |    | Effect unknown | A156!A, D168!DE |
|      |                                                           | 6  | 8 | L80R, D168AGHVY |
|      |                                                           |    | Effect unknown | L80!LR, A156!A, D168!ADGHVY |
|      | Grazoprevir (a component of Zepatier™)                    | 1A | 8 | Y56H, Q80KL, R155KQTW, A156GLTV, D168!CDMPQRSW |
|      |                                                           |    | 4 | V36ALM, Q41R, I48V, V55A, Y56F, Q80S, V107I, S122GT, I132V, R155GIS, V158A, D168C, I170V, T185S, E357GK |
|      |                                                           |    | Effect unknown | V36!AILMV, Y56!FHY, Q80!KLQRS, S122!AGRST, R155!GIKQRSTW, A156!AGLSTV, D168MPQRW |
|      |                                                           | 1B | 8 | Y56H, R155GTW, A156TV, D168AEFGHIKLTV |
|      |                                                           |    | 4 | T54S, Y56F, Q80L, V107I, S122NRT, V132IL, R155Q, D168NY, V170I, T185S, E357K |
|      |                                                           |    | Effect unknown | Y56!FHY, S122!AGNRST, V132!ILV, R155!EGKNQRSTW, A156!AGSTV, D168CMPQRW |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | 8 | A156TV |
|      |                                                           |    | 4 | Y56H, Q178R |
|      |                                                           |    | Effect unknown | A156!ATV, Q168!HQR |
|      |                                                           | 4  | 8 | D168AV |
|      |                                                           |    | 4 | V107I, A156MTV, D168EG, V170I |
|      |                                                           |    | Effect unknown | A156!AMTV, D168!ADEGV |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | 4 | V36IL, Y56H, L80KQ, I132L, A156MTV, D168CEY, I170V |
|      |                                                           |    | Effect unknown | V36!ILV, L80!KLQ, A156!AMTV, D168!CDEY |
|      | Paritaprevir (a component of Technivie™ and Holkira Pak™) | 1A | 8 | F43L, R155GKTW, A156T, D168!CDGMPQRSW |
|      |                                                           |    | 4 | V36AMT, V55I, Y56H, Q80KL, I132V, R155SV, A156G, P334S, S342P, E357K, V406AI, T449I, P470S |
|      |                                                           |    | Effect unknown | V36!ALMTV, F43!FL, Q80!KLQR, R155!GKMRSTVW, A156!AGT, D168CGMPQRSW, V406!AIV |
|      |                                                           | 1B | 8 | R155K, A156V, D168AFHKTVY |
|      |                                                           |    | 4 | Y56FH, A156T, D168ILN, E357K |
|      |                                                           |    | Effect unknown | R155!KR, A156!ASTV, D168CGMPQRSW |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | 8 | Y56H, R155C, A156TV, D168HV |
|      |                                                           |    | Effect unknown | Y56!HY, R155!CR, A156!ATV, D168!DHV |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Simeprevir (Galexos™)                                     | 1A | 8 | Q80KR, S122NR, R155KT, A156GTV, D168AEHV, I170T |
|      |                                                           |    | 4 | A5V, V36L, Q41R, T54S, V71A, Q80L, G120S, S122ACGIT, R123M, I132L, R155GQ, D168FT, I170V, N174RS, T178I, T343N, T344I, L356F, A477T, V511I, G534S, E555D, V629I |
|      |                                                           |    | Effect unknown | Q80!KLQR, S122!ACGINRST, R155!GKQRT, A156!AGTV, D168!ADEFHTV, I170!ITV |
|      |                                                           | 1B | 8 | F43ISV, Q80KR, S122R, R155GKQTW, A156GTV, D168AEHINQTVY |
|      |                                                           |    | 4 | A39T, T46S, V48I, S61L, Q80H, Q86R, R117C, S122AGIT, A156S, D168F, V170IT, S174F, S280P, T344A, T358F, G383S, I426V, T448A, I472V, V535A, P574AS, V629I |
|      |                                                           |    | Effect unknown | F43!FILSV, Q80!GHKLQR, S122!AGIRST, R155!GIKMQRTW, A156!AGSTV, D168CKLMPRSW, V170!AITV |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | 8 | R155K, A156G |
|      |                                                           |    | Effect unknown | Q80!Q, T122!T, R155!KR, A156!AGM, D168!D, V170!V |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Voxilaprevir (a component of Vosevi™)                     | 1A | 8 | R155W, A156ILMPT, D168KLR |
|      |                                                           |    | 4 | V36AG, Q41R, F43S, Q80K, R155G, A156GSV, D168AFITV |
|      |                                                           |    | Effect unknown | V36!AGILMV, R155!AGKRTW, A156!AGILMPSTV, D168CMPQW |
|      |                                                           | 1B | 8 | R155W, A156TV |
|      |                                                           |    | 4 | V36AM, S122D, A156S, D168VY, V170A |
|      |                                                           |    | Effect unknown | V36!AIMSV, R155!KRW, A156!ASTV, D168!ADEGVY |
|      |                                                           | 2  | 8 | A156LTV |
|      |                                                           |    | 4 | F43V |
|      |                                                           |    | Effect unknown | A156!AGLTV |
|      |                                                           | 3  | 8 | Q41K, Q80K, A156TV |
|      |                                                           |    | 4 | Q41R, V55A, Q168R, L175M, A294V |
|      |                                                           |    | Effect unknown | Q41!KQR, Q80!KQR, A156!ATV, Q168!HKQR |
|      |                                                           | 4  | 8 | A156LTV |
|      |                                                           |    | 4 | Q41R, D168ETV, T185S, T382S |
|      |                                                           |    | Effect unknown | A156!ALTV, D168!DEKTV |
|      |                                                           | 5  | 8 | D168HRY |
|      |                                                           |    | 4 | D168AK |
|      |                                                           |    | Effect unknown | A156!A, D168!ADEHKRVY |
|      |                                                           | 6  | 4 | Q41KR, Y56H, D168AEH |
|      |                                                           |    | Effect unknown | Q41!KQR, A156!A, D168!ADEHVY |
| NS5a | Daclatasvir (Daklinza™)                                   | 1A | 8 | K24T, M28AGT, Q30DEGHKNRY, L31IMV, P32L, H58D, A92K, Y93CHNST |
|      |                                                           |    | 4 | K24R, M28V, Q30T, H54RY, H58PV, A92P, Y93FR |
|      |                                                           |    | Effect unknown | M28!AGILMTV, Q30AFMSW, L31!ILMV, P32!LP, H54!CHNRY, H58AEFGIKMTW, A92!AKPTV, Y93!CFHNRSTY |
|      |                                                           | 1B | 8 | P29d, R30G, L31V, P32d, A92K, Y93H |
|      |                                                           |    | 4 | L28AMT, P29S, R30HKQ, L31FIMW, Q54H, P58LS, Q62E, A92E, Y93N |
|      |                                                           |    | Effect unknown | L28!AILMTV, R30!GHKLPQRST, L31!FILMVW, P32!LP, A92!AEKPTV, Y93!CFHLNY |
|      |                                                           | 2  | 8 | F28CS, L31M, C92R, Y93H |
|      |                                                           |    | Effect unknown | F28!CFILS, L31!LM, C92!CR, Y93!HY |
|      |                                                           | 3  | 8 | M28K, A30EK, L31FIMV, Y93H |
|      |                                                           |    | 4 | M28T, A30SV, L31P, P58H, S62AILPRT, E92G |
|      |                                                           |    | Effect unknown | M28!IKMTV, A30!AEKRSTV, L31!FILMPV, S62CEFGHMVWY, Y93!HY |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Elbasvir (a component of Zepatier™)                       | 1A | 8 | M28AGT, Q30DEGHKR, L31FIMV, H58D, Y93CHN |
|      |                                                           |    | 4 | K24R, M28V, Q30LY, E62D, Y93S, D427N |
|      |                                                           |    | Effect unknown | M28!AGMTV, Q30!DEGHKLQRY, L31!FILMV, H58!DH, Y93!CFHNSY |
|      |                                                           | 1B | 8 | R30Q, L31FMV, Y93H |
|      |                                                           |    | 4 | L28M, P58S |
|      |                                                           |    | Effect unknown | L31!FLMV, Y93!CHY |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | 8 | A30DK, L31FM, Y93H |
|      |                                                           |    | 4 | P58S, S62APT |
|      |                                                           |    | Effect unknown | A30!ADK, L31!FLM, P58!PS, S62!APST, Y93!HY |
|      |                                                           | 4  | 8 | L30FH, M31V, P58D, Y93H |
|      |                                                           |    | 4 | L28MST, L30RS, M31I, P32L, H54Y, P58TY, Y93C |
|      |                                                           |    | Effect unknown | L30!FHLPRS, M31!IMV, P58!DPTY, Y93!CHY |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Ledipasvir (a component of Harvoni™)                      | 1A | 8 | M28AGT, Q30EGHKNRY, L31IMV, P32L, S38F, H58D, A92K, Y93CHLNRSTW |
|      |                                                           |    | 4 | K24R, M28V, Q30LT, L31P, H58PR, Y93F |
|      |                                                           |    | Effect unknown | M28!AGILMTV, Q30ADFMPSW, L31!ILMPV, P32!LP, S38!FST, H58!CDHPQRY, A92!AKPTV, Y93ADEGIKMPQV |
|      |                                                           | 1B | 8 | P58D, A92K, Y93HN |
|      |                                                           |    | 4 | L28M, L31IMV, A92T, Y93CS |
|      |                                                           |    | Effect unknown | L31!FILMV, P58!DPQRST, A92!AEKTV, Y93!CFHLNSY |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Ombitasvir (a component of Technivie™ and Holkira Pak™)   | 1A | 8 | M28TV, Q30EKRY, L31V, H58D, Y93CFHLNS |
|      |                                                           |    | 4 | K24R, M28A, P32L, H54Y, H58PR, E62D |
|      |                                                           |    | Effect unknown | M28!AMTV, Q30!EHKLNQRTY, L31!ILMV, H58!DHPR, Y93!CFHLNSY |
|      |                                                           | 1B | 8 | L28T, P29d, Y93HN |
|      |                                                           |    | 4 | L28M, R30Q, L31FMV, Q54Y, P58S, Y93S |
|      |                                                           |    | Effect unknown | L28!LMT, P29!P, L31!FILMTV, Y93!HNSY |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | 8 | L28V |
|      |                                                           |    | 4 | L28S, M31I, P58S |
|      |                                                           |    | Effect unknown | L28!LSV |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Pibrentasvir (a component of Mavyret™)                    | 1A | 8 | M28G, Q30Dd, E62A, Y93HN |
|      |                                                           |    | 4 | K24EQR, M28ATV, P29QR, Q30EGHKLR, L31M, H58CDY |
|      |                                                           |    | Effect unknown | K24!EKQR, M28!AGMTV, P29!PQR, Q30!DEGHKLQRY, H58!CDHPRY, E62!ADE, Y93!CFHLNSY |
|      |                                                           | 1B | 8 | P32d |
|      |                                                           |    | 4 | L28M, P29d, R30Q, L31K, P32KNRT |
|      |                                                           |    | Effect unknown | P32!KNPRT |
|      |                                                           | 2  | 4 | L31M |
|      |                                                           | 3  | 4 | S24FG, M28GK, P29Q, A30GK, L31FIM, Y93H |
|      |                                                           |    | Effect unknown | S24!FGS, M28!GKMT, A30!AGK, L31!FILM, E92!E, Y93!HY |
|      | Velpatasvir (a component of Epclusa™ and Vosevi™)         | 1A | 8 | M28GT, Q30EK, L31IMV, P32L, H58D, A92K, Y93CHNRSW |
|      |                                                           |    | 4 | K24MT, M28AV, Q30GHLRT, L31F, S38Y, H58NPR, S85N, Y93LT, L153V |
|      |                                                           |    | Effect unknown | K24CDFHILPVWY, M28!AGILMTV, Q30ADFMNPW, L31!FILMV, P32!LP, H58!CDHLNPQRY, A92!AKPTV, Y93ADEGIKMPQV |
|      |                                                           | 1B | 8 | A92K, Y93HN |
|      |                                                           |    | 4 | Q24K, L31IMV, P58T, Y93CRST |
|      |                                                           |    | Effect unknown | L31!FILMV, A92!AEKPTV, Y93!CFHLNRSTY |
|      |                                                           | 2  | 8 | F28S, L31V, C92T, Y93HN |
|      |                                                           |    | 4 | F28C, L31IM, P58A, C92RS, Y93F |
|      |                                                           |    | Effect unknown | F28!CFLSV, L31!ILMV, C92!ACKNRST, Y93!CFHLNSTY |
|      |                                                           | 3  | 8 | M28T, A30K, L31FM, Y93HS |
|      |                                                           |    | 4 | S14E, M28V, A30HRV, L31PV, S38PY, V52M, P58GLT, H85Y, E92K, Y93NR, S103DP, A291P, S379P, S385Q, R404K, P407S |
|      |                                                           |    | Effect unknown | M28!LMTV, A30!AHKQRSV, L31!FLMPV, S38!PSY, P58!AGHLPST, Y93!FHNRSY |
|      |                                                           | 4  | 8 | L28T, Y93HNS |
|      |                                                           |    | 4 | L30R, M31V, P32L, H54R, Y93CW |
|      |                                                           |    | Effect unknown | L28!LMTV, Y93!CHNSWY |
|      |                                                           | 5  | 4 | L31I |
|      |                                                           | 6  | 8 | L31V, P32ALQR |
|      |                                                           |    | 4 | F28MV, L31IM, T58GH, A92T, T93AHNS |
|      |                                                           |    | Effect unknown | F28!AFMV, L31!ILMV, P32!ALPQR, T58!AGHPST, T93!ACFHLNSTY |
| NS5b | Dasabuvir (a component of Holkira Pak™)                   | 1A | 8 | L314H, C316Y, A395G, M414ITV, N444K, E446KQ, Y448CH, A553TV, G554S, S556GR, Y561H |
|      |                                                           |    | 4 | G307R, A450V, A553I, G558R, D559GINV, L588F |
|      |                                                           |    | Effect unknown | L314!HL, C316!CY, A395!AG, M414!IMTV, N444!KN, E446!EKQ, Y448!CHY, A553!AITV, G554!DGS, S556!GRS, D559!DGINV, Y561!HY |
|      |                                                           | 1B | 8 | C316HNY, S368T, N411S, M414ITV, Y448CH, A553V, S556G, D559G |
|      |                                                           |    | 4 | C445F, S556R |
|      |                                                           |    | Effect unknown | C316!CHNY, A395!A, M414!IMTV, D444!D, Y448!CHY, A553!AV, S556!GRS, D559!DG |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |
|      | Sofosbuvir (a component of Epclusa™ and Vosevi™)          | 1A | 8 | S282GRT |
|      |                                                           |    | 4 | D61G, S62N, D66G, F101L, A112T, N142T, L159F, E202G, E237G, C316HR, L320IV, V321AFI, Q355H, F415Y, S473T, N590d |
|      |                                                           |    | Effect unknown | L159!FL, S282!GRST, C316!CHR, L320!ILV, V321!AFIV |
|      |                                                           | 1B | 8 | S282GRT |
|      |                                                           |    | 4 | N142S, K355T |
|      |                                                           |    | Effect unknown | L159!L, S282!GRST, C316!C, L320!L, V321!V |
|      |                                                           | 2  | 8 | S282GRT |
|      |                                                           |    | Effect unknown | L159!L, S282!GRST, C316!C, L320!L, V321!V |
|      |                                                           | 3  | 8 | S282GRT |
|      |                                                           |    | 4 | K100R, R120C, N142ST, L159F, E237G, L314FIP, L320FI, V321A, Q355HR |
|      |                                                           |    | Effect unknown | L159!FL, S282!GRST, C316!C, L320!FIL, V321!AV |
|      |                                                           | 4  | 8 | S282GRT |
|      |                                                           |    | 4 | E237G, V321I |
|      |                                                           |    | Effect unknown | L159!L, S282!GRST, C316!C, L320!L, V321!IV |
|      |                                                           | 5  | 8 | S282GRT |
|      |                                                           |    | 4 | M289I |
|      |                                                           |    | Effect unknown | L159!L, S282!GRST, C316!C, L320!L, V321!V |
|      |                                                           | 6  | 8 | S282GRT |
|      |                                                           |    | Effect unknown | L159!L, S282!GRST, C316!C, L320!L, V321!V |
|      |                                                           | 6E | Not available | TRUE |
|      | Sofosbuvir (a component of Harvoni™)                      | 1A | 8 | S282GRT |
|      |                                                           |    | 4 | D61G, S62N, D66G, F101L, A112T, N142T, L159F, E202G, E237G, C316HR, L320IV, V321AFI, Q355H, F415Y, S473T, N590d |
|      |                                                           |    | Effect unknown | L159!FL, S282!GRST, C316!CHR, L320!ILV, V321!AFIV |
|      |                                                           | 1B | 8 | S282GRT |
|      |                                                           |    | 4 | N142S, K355T |
|      |                                                           |    | Effect unknown | L159!L, S282!GRST, C316!C, L320!L, V321!V |
|      |                                                           | 2  | Not indicated | TRUE |
|      |                                                           | 3  | Not indicated | TRUE |
|      |                                                           | 4  | Not indicated | TRUE |
|      |                                                           | 5  | Not indicated | TRUE |
|      |                                                           | 6  | Not indicated | TRUE |

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
