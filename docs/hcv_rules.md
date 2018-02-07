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
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Grazoprevir (a component of Zepatier™)                    | 1A | SCORE FROM ( 36ALM=>4, 56H=>8, 56F=>4, 107I=>4, 122GT=>4, 122R=>0, 155GIKS=>4, 156GLTV=>8, 156M=>4, 156S=>0, 158A=>4, 168AEFGHIKLTVY=>8, 168CNS=>4, 170V=>4, 170T=>0, D168!AEFGHIKLTVYCNSD=>"Effect unknown", R155!GIKSR=>"Effect unknown", V36!ALMV=>"Effect unknown", Y56!HFY=>"Effect unknown", S122!GTRS=>"Effect unknown", V107!IV=>"Effect unknown", A156!GLTVMSA=>"Effect unknown", I170!VTI=>"Effect unknown", V158!AV=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 56H=>8, 56F=>4, 107I=>4, 155GTW=>8, 155K=>0, 156TV=>8, 156G=>4, 168AFGHIKLTV=>8, 168ENY=>4, Y56!HFY=>"Effect unknown", A156!TVGA=>"Effect unknown", V107!IV=>"Effect unknown", R155!GTWKR=>"Effect unknown", D168!AFGHIKLTVENYD=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( 77S=>8, 168QR=>8, 178R=>8, Q168!QRQ=>"Effect unknown", Q178!RQ=>"Effect unknown", V170!V=>"Effect unknown", Y56!Y=>"Effect unknown", S122!S=>"Effect unknown", R155!R=>"Effect unknown", A156!A=>"Effect unknown", N77!SN=>"Effect unknown" ) |
|      |                                                           | 4  | SCORE FROM ( 107I=>4, 122G=>4, 156GMTV=>4, 158I=>4, 168AV=>8, 168CEGNY=>4, 170I=>4, D168!AVCEGNYD=>"Effect unknown", R155!R=>"Effect unknown", Y56!Y=>"Effect unknown", T122!GT=>"Effect unknown", V107!IV=>"Effect unknown", A156!GMTVA=>"Effect unknown", V170!IV=>"Effect unknown", V158!IV=>"Effect unknown" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Paritaprevir (a component of Technivie™ and Holkira Pak™) | 1A | SCORE FROM ( 36AMT=>4, 36L=>0, 43L=>8, 55I=>4, 56H=>4, 80KL=>4, 80R=>0, 132V=>4, 155GKTW=>8, 155S=>4, 156T=>8, 156GVS=>4, 168AEFHILNTVY=>8, 334S=>4, 342P=>4, 357K=>4, 406AI=>4, 449I=>4, 470S=>4, T449!IT=>"Effect unknown", V36!AMTLV=>"Effect unknown", E357!KE=>"Effect unknown", V406!AIV=>"Effect unknown", D168!AEFHILNTVYD=>"Effect unknown", P470!SP=>"Effect unknown", F43!LF=>"Effect unknown", P334!SP=>"Effect unknown", Q80!KLRQ=>"Effect unknown", S342!PS=>"Effect unknown", V55!IV=>"Effect unknown", Y56!HY=>"Effect unknown", I132!VI=>"Effect unknown", R155!GKTWSR=>"Effect unknown", A156!TGVSA=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 56H=>4, 155K=>8, 156V=>8, 156T=>4, 156S=>0, 168AFHKTVY=>8, 168ILN=>4, 168E=>0, 357K=>4, Y56!HY=>"Effect unknown", D168!AFHKTVYILNED=>"Effect unknown", R155!KR=>"Effect unknown", A156!VTSA=>"Effect unknown", E357!KE=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( 56H=>8, 155C=>8, 156TV=>8, 168HV=>8, Y56!HY=>"Effect unknown", D168!HVD=>"Effect unknown", R155!CR=>"Effect unknown", A156!TVA=>"Effect unknown" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Simeprevir (Galexos™)                                     | 1A | SCORE FROM ( 36L=>4, 43LIV=>4, 54S=>4, 80KR=>8, 80GHMIL=>4, 122GR=>8, 122ANIT=>4, 155KT=>8, 155GIMQSW=>4, 156GTV=>8, 156LMS=>4, 168AEHV=>8, 168CFGIKLNSTY=>4, 170T=>8, 170V=>4, Q80!KRGHMILQ=>"Effect unknown", R155!KTGIMQSWR=>"Effect unknown", V36!LV=>"Effect unknown", T54!ST=>"Effect unknown", D168!AEHVCFGIKLNSTYD=>"Effect unknown", S122!GRANITS=>"Effect unknown", F43!LIVF=>"Effect unknown", A156!GTVLMSA=>"Effect unknown", I170!TVI=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 41R=>8, 43SIV=>8, 43L=>4, 80KR=>8, 80HIM=>4, 80GL=>0, 122AIRT=>8, 155KGQTW=>8, 155C=>4, 155IM=>0, 156GTV=>8, 156S=>4, 168AEFGHINQTVY=>8, 168KL=>4, 170T=>4, 170A=>0, Q80!KRHIMGLQ=>"Effect unknown", R155!KGQTWCIMR=>"Effect unknown", D168!AEFGHINQTVYKLD=>"Effect unknown", Q41!RQ=>"Effect unknown", S122!AIRTS=>"Effect unknown", F43!SIVLF=>"Effect unknown", A156!GTVSA=>"Effect unknown", V170!TAV=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( 155K=>8, 156G=>8, 156M=>0, R155!KR=>"Effect unknown", A156!GMA=>"Effect unknown" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Telaprevir (Incivek™)                                     | 1A | SCORE FROM ( 36AGILM=>8, 54AS=>8, 132V=>4, 155GMTK=>8, 156STV=>8, 168N=>4, I132!VI=>"Effect unknown", T54!AST=>"Effect unknown", D168!ND=>"Effect unknown", V36!AGILMV=>"Effect unknown", R155!GMTKR=>"Effect unknown", A156!STVA=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 36GILMA=>8, 54SA=>8, 155K=>8, 156FTVS=>8, A156!FTVSA=>"Effect unknown", R155!KR=>"Effect unknown", V36!GILMAV=>"Effect unknown", T54!SAT=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
| NS5a | Daclatasvir (Daklinza™)                                   | 1A | SCORE FROM ( 24RG=>4, 24N=>0, 26E=>0, 28AT=>8, 28SV=>4, 30EHKNRYDG=>8, 30T=>4, 30L=>0, 31IMV=>8, 54R=>4, 58D=>8, 58PR=>4, 92K=>0, 93CHNS=>8, 93F=>0, H54!RH=>"Effect unknown", Y93!CHNSFY=>"Effect unknown", K24!RGNK=>"Effect unknown", A92!KA=>"Effect unknown", K26!EK=>"Effect unknown", M28!ATSVM=>"Effect unknown", H58!DPRH=>"Effect unknown", Q30!EHKNRYDGTLQ=>"Effect unknown", L31!IMVL=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 28MT=>4, 29XS=>4, 30GHPQ=>4, 31FV=>8, 31IM=>4, 32X=>8, 32L=>4, 58S=>4, 62D=>4, 92K=>4, 93H=>8, 93N=>4, P32!XLP=>"Effect unknown", Q62!DQ=>"Effect unknown", Y93!HNY=>"Effect unknown", A92!KA=>"Effect unknown", P58!SP=>"Effect unknown", L28!MTL=>"Effect unknown", P29!XSP=>"Effect unknown", R30!GHPQR=>"Effect unknown", L31!FVIML=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( 28SC=>8, 31M=>8, 92R=>8, 93H=>8, C92!RC=>"Effect unknown", F28!SCF=>"Effect unknown", Y93!HY=>"Effect unknown", L31!ML=>"Effect unknown" ) |
|      |                                                           | 3  | SCORE FROM ( 28T=>0, 28V=>0, 30K=>8, 30S=>4, 30TV=>0, 31FIMV=>8, 31P=>4, 62ALIRPT=>4, 93H=>8, Y93!HY=>"Effect unknown", M28!TVM=>"Effect unknown", S62!ALIRPTS=>"Effect unknown", A30!KSTVA=>"Effect unknown", L31!FIMVPL=>"Effect unknown" ) |
|      |                                                           | 4  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Elbasvir (a component of Zepatier™)                       | 1A | SCORE FROM ( 28AGT=>8, 28V=>0, 30DEGHR=>8, 30KY=>4, 30L=>0, 31IMFV=>8, 58D=>8, 93CHNS=>8, H58!DH=>"Effect unknown", M28!AGTVM=>"Effect unknown", Y93!CHNSY=>"Effect unknown", Q30!DEGHRKYLQ=>"Effect unknown", L31!IMFVL=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 28M=>4, 31FMV=>8, 93H=>8, L28!ML=>"Effect unknown", Y93!HY=>"Effect unknown", L31!FMVL=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( 30DK=>8, 31FM=>8, 93H=>8, Y93!HY=>"Effect unknown", A30!DKA=>"Effect unknown", L31!FML=>"Effect unknown" ) |
|      |                                                           | 4  | SCORE FROM ( 30HF=>8, 30S=>0, 31V=>8, 31I=>4, 32L=>4, 58D=>8, 93H=>8, P32!LP=>"Effect unknown", P58!DP=>"Effect unknown", Y93!HY=>"Effect unknown", L30!HFSL=>"Effect unknown", M31!VIM=>"Effect unknown" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Ledipasvir (a component of Harvoni™)                      | 1A | SCORE FROM ( 24RGN=>4, 28AGT=>8, 28V=>4, 30EGHKNRY=>8, 30LT=>4, 31IMV=>8, 31P=>4, 32L=>8, 38F=>8, 58D=>8, 58P=>4, 92P=>8, 92T=>4, 93CHNS=>8, 93F=>4, P32!LP=>"Effect unknown", S38!FS=>"Effect unknown", K24!RGNK=>"Effect unknown", A92!PTA=>"Effect unknown", H58!DPH=>"Effect unknown", M28!AGTVM=>"Effect unknown", Y93!CHNSFY=>"Effect unknown", Q30!EGHKNRYLTQ=>"Effect unknown", L31!IMVPL=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 28M=>4, 31IV=>8, 31M=>4, 31F=>0, 32L=>0, 58D=>8, 92K=>8, 92T=>4, 93H=>8, 93C=>4, 93S=>0, P32!LP=>"Effect unknown", L28!ML=>"Effect unknown", P58!DP=>"Effect unknown", A92!KTA=>"Effect unknown", Y93!HCSY=>"Effect unknown", L31!IVMFL=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Ombitasvir (a component of Technivie™ and Holkira Pak™)   | 1A | SCORE FROM ( 24R=>4, 28TV=>8, 28A=>4, 30EKRY=>8, 30HLNT=>0, 31V=>8, 31MI=>0, 32L=>4, 54Y=>4, 58D=>8, 58PR=>4, 93CFHLNS=>8, P32!LP=>"Effect unknown", H54!YH=>"Effect unknown", K24!RK=>"Effect unknown", H58!DPRH=>"Effect unknown", M28!TVAM=>"Effect unknown", Y93!CFHLNSY=>"Effect unknown", Q30!EKRYHLNTQ=>"Effect unknown", L31!VMIL=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 28T=>8, 28M=>4, 29X=>8, 30GHPQ=>4, 31FV=>8, 31MT=>4, 31I=>0, 32X=>4, 54Y=>4, 58AS=>4, 92E=>4, 93HNS=>8, P32!XP=>"Effect unknown", Q54!YQ=>"Effect unknown", Y93!HNSY=>"Effect unknown", A92!EA=>"Effect unknown", P58!ASP=>"Effect unknown", L28!TML=>"Effect unknown", P29!XP=>"Effect unknown", R30!GHPQR=>"Effect unknown", L31!FVMTIL=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( 28V=>8, 28S=>4, 31I=>4, 31L=>0, 58S=>4, 58AP=>0, P58!SAPP=>"Effect unknown", L28!VSL=>"Effect unknown", M31!ILM=>"Effect unknown" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Velpatasvir (a component of Epclusa™)                     | 1A | SCORE FROM ( 24RTM=>4, 28GT=>8, 28V=>4, 30EK=>8, 30HLR=>4, 31MV=>8, 31I=>4, 32L=>8, 58D=>8, 58P=>4, 92K=>4, 93CHNRSW=>8, P32!LP=>"Effect unknown", K24!RTMK=>"Effect unknown", A92!KA=>"Effect unknown", H58!DPH=>"Effect unknown", M28!GTVM=>"Effect unknown", Y93!CHNRSWY=>"Effect unknown", Q30!EKHLRQ=>"Effect unknown", L31!MVIL=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 24K=>0, 30S=>0, 31MVI=>4, 58RT=>0, 92K=>8, 93CHNS=>4, Q24!KQ=>"Effect unknown", P58!RTP=>"Effect unknown", A92!KA=>"Effect unknown", Y93!CHNSY=>"Effect unknown", R30!SR=>"Effect unknown", L31!MVIL=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( 28S=>8, 31V=>8, 31MI=>4, 31L=>0, 58A=>0, 58S=>0, 92T=>8, 93HN=>8, C92!TC=>"Effect unknown", P58!ASP=>"Effect unknown", F28!SF=>"Effect unknown", Y93!HNY=>"Effect unknown", L31!VMILL=>"Effect unknown" ) |
|      |                                                           | 3  | SCORE FROM ( 30K=>8, 30V=>4, 31M=>8, 31VP=>4, 38P=>4, 58L=>4, 92K=>4, 93HS=>8, 93NR=>4, S38!PS=>"Effect unknown", P58!LP=>"Effect unknown", E92!KE=>"Effect unknown", Y93!HSNRY=>"Effect unknown", A30!KVA=>"Effect unknown", L31!MVPL=>"Effect unknown" ) |
|      |                                                           | 4  | SCORE FROM ( 28T=>4, 30HS=>4, 31V=>4, 32L=>4, 93CHNSW=>4, P32!LP=>"Effect unknown", L28!TL=>"Effect unknown", Y93!CHNSWY=>"Effect unknown", L30!HSL=>"Effect unknown", M31!VM=>"Effect unknown" ) |
|      |                                                           | 5  | SCORE FROM ( L28!L=>"Effect unknown", P58!P=>"Effect unknown", A92!A=>"Effect unknown", T93!T=>"Effect unknown", Q30!Q=>"Effect unknown", L31!L=>"Effect unknown" ) |
|      |                                                           | 6  | SCORE FROM ( 31V=>8, 32ALQR=>8, P32!ALQRP=>"Effect unknown", A92!A=>"Effect unknown", T58!T=>"Effect unknown", F28!F=>"Effect unknown", T93!T=>"Effect unknown", R30!R=>"Effect unknown", L31!VL=>"Effect unknown" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>8 ) |
| NS5b | Dasabuvir (a component of Holkira Pak™)                   | 1A | SCORE FROM ( 307R=>4, 314H=>8, 316Y=>8, 395G=>8, 414ITV=>8, 444K=>8, 444D=>0, 446KQ=>8, 448CH=>8, 450V=>4, 553TV=>8, 553I=>4, 554S=>8, 556GR=>8, 558R=>4, 559GINV=>4, 561H=>8, Y448!CHY=>"Effect unknown", A450!VA=>"Effect unknown", A553!TVIA=>"Effect unknown", G554!SG=>"Effect unknown", A395!GA=>"Effect unknown", S556!GRS=>"Effect unknown", G558!RG=>"Effect unknown", D559!GINVD=>"Effect unknown", Y561!HY=>"Effect unknown", G307!RG=>"Effect unknown", M414!ITVM=>"Effect unknown", C316!YC=>"Effect unknown", L314!HL=>"Effect unknown", N444!KDN=>"Effect unknown", E446!KQE=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 307R=>0, 316HNY=>8, 316W=>4, 368T=>8, 411S=>8, 414ITV=>8, 445FY=>4, 448CH=>8, 553V=>8, 556G=>8, 556R=>4, 559G=>8, S368!TS=>"Effect unknown", Y448!CHY=>"Effect unknown", K307!RK=>"Effect unknown", D559!GD=>"Effect unknown", A553!VA=>"Effect unknown", N411!SN=>"Effect unknown", C316!HNYWC=>"Effect unknown", C445!FYC=>"Effect unknown", M414!ITVM=>"Effect unknown", S556!GRS=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not indicated" ) |
|      | Sofosbuvir (a component of Epclusa™)                      | 1A | SCORE FROM ( 61G=>4, 112T=>4, 159F=>4, 237G=>4, 282T=>8, 320V=>4, 321I=>4, 321A=>4, 355H=>4, 473T=>4, A112!TA=>"Effect unknown", L320!VL=>"Effect unknown", Q355!HQ=>"Effect unknown", V321!IAV=>"Effect unknown", S473!TS=>"Effect unknown", S282!TS=>"Effect unknown", L159!FL=>"Effect unknown", C316!C=>"Effect unknown", E237!GE=>"Effect unknown", D61!GD=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 142S=>4, 282T=>8, 321I=>4, L320!L=>"Effect unknown", V321!IV=>"Effect unknown", S282!TS=>"Effect unknown", C316!C=>"Effect unknown", N142!SN=>"Effect unknown", L159!L=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( 282T=>8, L320!L=>"Effect unknown", V321!V=>"Effect unknown", S282!TS=>"Effect unknown", C316!C=>"Effect unknown", L159!L=>"Effect unknown" ) |
|      |                                                           | 3  | SCORE FROM ( 142T=>4, 159F=>4, 237G=>4, 282T=>8, 314IFP=>4, 321A=>4, 355H=>4, L320!L=>"Effect unknown", V321!AV=>"Effect unknown", Q355!HQ=>"Effect unknown", S282!TS=>"Effect unknown", L159!FL=>"Effect unknown", C316!C=>"Effect unknown", L314!IFPL=>"Effect unknown", N142!TN=>"Effect unknown", E237!GE=>"Effect unknown" ) |
|      |                                                           | 4  | SCORE FROM ( 237G=>4, 282T=>8, 321I=>4, L320!L=>"Effect unknown", V321!IV=>"Effect unknown", S282!TS=>"Effect unknown", C316!C=>"Effect unknown", E237!GE=>"Effect unknown", L159!L=>"Effect unknown" ) |
|      |                                                           | 5  | SCORE FROM ( 282T=>8, 289I=>4, L320!L=>"Effect unknown", M289!IM=>"Effect unknown", V321!V=>"Effect unknown", S282!TS=>"Effect unknown", C316!C=>"Effect unknown", L159!L=>"Effect unknown" ) |
|      |                                                           | 6  | SCORE FROM ( 282T=>8, L320!L=>"Effect unknown", V321!V=>"Effect unknown", S282!TS=>"Effect unknown", C316!C=>"Effect unknown", L159!L=>"Effect unknown" ) |
|      | Sofosbuvir (a component of Harvoni™)                      | 1A | SCORE FROM ( 61G=>4, 112T=>4, 159F=>4, 237G=>4, 282T=>8, 320V=>4, 321I=>4, 321A=>4, 355H=>4, 473T=>4, A112!TA=>"Effect unknown", L320!VL=>"Effect unknown", Q355!HQ=>"Effect unknown", V321!IAV=>"Effect unknown", S473!TS=>"Effect unknown", S282!TS=>"Effect unknown", L159!FL=>"Effect unknown", C316!C=>"Effect unknown", E237!GE=>"Effect unknown", D61!GD=>"Effect unknown" ) |
|      |                                                           | 1B | SCORE FROM ( 142S=>4, 282T=>8, 321I=>4, L320!L=>"Effect unknown", V321!IV=>"Effect unknown", S282!TS=>"Effect unknown", C316!C=>"Effect unknown", N142!SN=>"Effect unknown", L159!L=>"Effect unknown" ) |
|      |                                                           | 2  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 3  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 4  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 5  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6  | SCORE FROM ( TRUE=>"Not indicated" ) |
|      |                                                           | 6E | SCORE FROM ( TRUE=>"Not available" ) |
