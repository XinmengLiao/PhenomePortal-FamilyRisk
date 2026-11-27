## 🔹 Cohort/Population-level Newborn Screening 

### 1. Upload genomic data and set configurations
- From the **Home** page, choose **"Newborn Screening - Cohort"**. Click **"Continue"**. 
- Input sample name, upload genomic data in compressed VCF format (.vcf.gz), upload sample metadata including sample ID and sex.
- Upload customized screening panel (.txt/.tsv) or select one or more preset screening lists. 
- Select reference genome version.
- Select variants with high sequencing quality and/or reported in ClinVar in **“Variant selection"**.
- Set the parameters for filtering variants in **"Advanced variant filters"**: (global minor allele frequency, predicted scores, clinical pathogenicity). 

    For explortaory analyses in the **"Research options"**: 
- Polygenic risk score calculation: enter any of the PGS Catalog Identifiers if needed. 
- Pharmacogenomics: select to run PGx profiling or not. 

<br>

### 2. View Results
The submitted job and job status are presented in **"Analysis Jobs"** on the homepage. \
Once the job is completed, click the job name to view all results.

The newborn screening of population level is for research purpose only, with each section tagged as "RESEARCH". 

Results sections:
- Sample & Variant Summary: This section summarizes the total number of screened newborns, filtered variants, and the screened ClinVar P/LP variants. 
- Screen-positive Results: This section presents the potential monogenic rare diseases identified in all newborns. These diseases are caused by both ClinVar and ACMG P/LP variants. 
- Disease Carrier Status: This section presents the potential disease carriers identified in all newborns. These diseases are caused by both ClinVar and ACMG P/LP variants, which are heterozygous states and associated with recessive diseases.
- Statistics: 
  - Pathogenic variant stats: The total number of pathogenic variants that associated with screen-positive or carrier status findings. 
  - Ontology stats: the total number of diseases that associated with screen-positive or carrier status findings. 
  - Cohort-level allele frequency: the allele frequency of each identified variant within the cohort. 
 - Variant Details: This section presents all details for the filtered variants.
- Pharmacogenomics: This section provides a report for each family member summarizing all identified pharmacogenomics diplotypes in the newborns via PharmCat. Detailed explanations and recommendations are illustrated in the structured report. 
- Polygenic Risk Scores: This section provides polygenic scores of each family member and reference population regarding the required PGS model. The similarity scores are also provided for genetic clustering. Information is also summarized into a strctured report.  


<br>

### 3. Export Data  
- Results in each section could be downloaded separately.
    - <img width="477" height="66" alt="image" src="https://github.com/user-attachments/assets/de80c455-eddf-4a33-84ec-09db955c08cd" />
- All results could be wrapped up and downloaded as a zipped file by clicking the download button in the bottom-right corner.
    - <img width="437" height="69" alt="image" src="https://github.com/user-attachments/assets/08035174-5529-4e9e-b81b-d47ce5eaea93" />