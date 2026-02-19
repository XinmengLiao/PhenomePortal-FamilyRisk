## 🔹 Family-level Newborn Screening 

### 1. Upload genomic data and set configurations
- From the **Home** page, choose **“Newborn Screening - Cohort/Population**. Click **"Continue"**. 
- Input sample name, upload genomic data in compressed VCF format (.vcf.gz), upload sample metadata including sample ID and sex.
- Upload customized screening panel (.txt/.tsv) or select one or more preset screening lists. 
- Select reference genome version.
- Set the parameters for filtering variants (global minor allele frequency, predicted scores, clinical pathogenicity). 
- Enter any of the PGS Catalog Identifiers if needed. 
- Select to run PGx profiling or not. 
- Click **"Submit"** to proceed.  

<br>

### 2. View Results
The submitted job and job status are presented in **"Analysis Jobs"** on the homepage. \
Once the job is completed, click the job name to view all results. \
Results are organized into 8 sections:
- Cohort Information: This section summarizes the total number of screened newborns. 
- Variant Summary: This section summarizes screened clinical P/LP variants (both ClinVar or ACMG P/LP). 
- Screen-positive Results: This section presents the potential monogenic rare diseases identified in all newborns. These diseases are caused by both ClinVar and ACMG P/LP variants. 
- Disease Carrier Status: This section presents the potential disease carriers identified in all newborns. These diseases are caused by both ClinVar and ACMG P/LP variants, which are heterozygous states and associated with recessive diseases.
- Cohort-level allele frequency: the allele frequency of each identified variant within the cohort. 
- Pharmacogenomics: This section provides a report for each family member summarizing all identified pharmacogenomics diplotypes in the newborns via PharmCat. Detailed explanations and recommendations are illustrated in the structured report. 
- Polygenic Risk Scores: This section provides polygenic scores of each family member and reference population regarding the required PGS model. The similarity scores are also provided for genetic clustering. Information is also summarized into a strctured report.  
- Variant Details: This section presents all details for the filtered variants.

<br>

### 3. Export Data  
- Results in each section could be downloaded separately.
    - <img width="477" height="66" alt="image" src="https://github.com/user-attachments/assets/de80c455-eddf-4a33-84ec-09db955c08cd" />
- All results could be wrapped up and downloaded as a zipped file by clicking the download button in the bottom-right corner.
    - <img width="437" height="69" alt="image" src="https://github.com/user-attachments/assets/08035174-5529-4e9e-b81b-d47ce5eaea93" />