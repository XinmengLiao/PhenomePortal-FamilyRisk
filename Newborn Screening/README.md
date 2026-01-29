Things need to be improved in the Home Page. 

### --- 1 ---
Integrate the FamilyRisk Knowledgebase into the Home Page. In this knowledgebase, users could check the information of Newborn Screening Projects, and generate their own screening list with the "Screening List Generator" tool. \
The knowledge base is deveploped via Rshinny and currently built on Scilifelab Server. 

### --- 2 ---
For every analysis where PGS calculation is available, more options need to be added: 
1. Experimental Factor Ontology ID (EFO), e.g. EFO_0004838
2. PGS publication ID (PGP), e.g. PGP000128
3. PGS ID (PGS), e.g. PGS000668 (This one already exist)

Both functions for carrier screening also have PGS calculation. Pleas added the abovementioned options. 

### --- 3 ---
For both carrier and newborn screening, add options for users to select the gene-disease screening list. 
1. It is optional for users to upload or select the screening panel. If no screening panel is uploaded or selected, the default screening panels will be 
2. Users could choose one or more screening panels from: 
    Newborn screening projects: ACMGv3.3, BabySeq_GroupA, BabySeq_GroupB, EarlyCheck_Group2, EarlyCheck_Group1, EarlyCheck_Group3, EarlyCheck_Group4, Genomic101, NBScreening, BabyScreen+, BabyDetect, Guardian_Group1, Guardian_Group2 \
    Expanded carrier screening list: ACMG_Carrier_Tier_1, ACMG_Carrier_Tier_2, ACMG_Carrier_Tier_3, ACMG_Carrier_Tier_4, ACMG_Carrier_Outside_gnomAD
3. Users could upload their own screening list in the `.txt` format of: \
    `Genes\tDisease\tInheritance\tMIM\tDisorder_Group\tProject`
    If users generate the screneing list from FamilyRisk knowledgebase, the list is already in the right format and ready to be uploaded.

### --- 4 ---
"Help & Support" section for tutorial documentation and video. 
