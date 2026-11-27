## 🔹 Family-level Carrier Screening 

### 1. Upload genomic data and set configurations
- From the **Home** page, choose **“Carrier Screening - Family**. Click **"Continue"**. 
- Input sample name, upload genomic data in compressed VCF format (.vcf.gz), and upload PED file illustrating the parents' relationship.  
- Upload customized screening panel (.txt/.tsv) or select one or more preset screening lists. 
- Select reference genome version.
- Select variants with high sequencing quality and/or reported in ClinVar in **“Variant selection"**.
- Set the parameters for filtering variants in **"Advanced variant filters"**: (global minor allele frequency, predicted scores, clinical pathogenicity). 

    For explortaory analyses in the **"Research options"**: 
- Polygenic risk score calculation: enter any of the PGS Catalog Identifiers if needed. 
- Pharmacogenomics: select to run PGx profiling or not. 

    Click **"Submit"** to proceed.  

<br>

### 2. View Results
The submitted job and job status are presented in **"Analysis Jobs"** on the homepage. \
Once the job is completed, click the job name to view all results. 

FamilyRisk offers **"Clinical Mode"** for monogenic disease screening and clinical reporting, and **"Research Mode"** for exploratory analyses of complex traits by polygenic risk scores and drug response by pharmacogenomics. The results are clearly seperated with tag **"CLINICAL"** or **"RESEARCH"**. 

Results sections of clinical mode:
- Sample & Variant Summary: This section summarizes 1) the family members' relationships from PED file; 2) the total number of screened variant and clinical P/LP variants (both ClinVar or ACMG P/LP); 3) the total number of positive findings of rare monogenic diseases and carrier status. 
- Screen-positive Results: This section presents the recessive disease that would be potentially transmitted to the offspring. These diseases are caused by both ClinVar and ACMG P/LP variants. 
- Findings by Person: This section summarises how variants are distributed across family members, covering both the variants in screen-positive findings and all screened variants.
- Variant Details: This section presents all details for the filtered variants.

Results sections of research mode:
- Pharmacogenomics: This section provides a report for each family member summarizing all identified pharmacogenomics diplotypes in the newborns via PharmCat. Detailed explanations and recommendations are illustrated in the structured report. 
- Polygenic Risk Scores: This section provides polygenic scores of each family member and reference population regarding the required PGS model. The similarity scores are also provided for genetic clustering. Information is also summarized into a strctured report.  

<br>

### 3. Choose pathogenic variants and potential diseases to be included in the final PDF report. 
FamilyRisk provides a template for clinical reporting, where users can further customize the report by selecting specific results and adding professional advice, suggestions, and comments. 
1) In **"Screen-positive Results"**, **"Disease Carrier Status"**, and **"Variant Details"** sections, switch the button under **"Include"** to include or exclude the specific variant to be presented in the final report. 
<img width="550" height="238" alt="image" src="https://github.com/user-attachments/assets/27b7866c-a698-4b26-8ba1-47c1c3631f31" />

2) Click **"Report Generation"** to generate the report. 
<img width="338" height="118" alt="image" src="https://github.com/user-attachments/assets/10e63348-038e-4eba-af07-502a32fc6765" />

3) Add additional analytical texts, recommendations, and references via the portal. 
4) Set the report language. Options: English, Turkish, Spanish, Korean, Chinese, Japanese. 
5) Click **"Continue"** to view the final digital report.
<img width="535" height="676" alt="image" src="https://github.com/user-attachments/assets/80c10c72-4d4c-49b7-a38a-8b805517ccd9" />

6) Download, print, or re-edit the report following the previous steps.

<br>
<br>

### 4. Export Data  
- Results in each section could be downloaded separately.
    - <img width="477" height="66" alt="image" src="https://github.com/user-attachments/assets/de80c455-eddf-4a33-84ec-09db955c08cd" />
- All results could be wrapped up and downloaded as a zipped file by clicking the download button in the bottom-right corner.
     - <img width="437" height="69" alt="image" src="https://github.com/user-attachments/assets/08035174-5529-4e9e-b81b-d47ce5eaea93" />