## 🔹 Individual-level Newborn Screening 

### 1. Upload genomic data and set configurations
- From the **Home** page, choose **“Newborn Screening - Individual”**. Click **"Continue"**. 
- Input sample name, upload genomic data in compressed VCF format (.vcf.gz). 
- Upload customized screening panel (.txt/.tsv) or select one or more preset screening lists. 
- Select gender, reference genome version.
- Set the parameters for filtering variants (global minor allele frequency, predicted scores, clinical pathogenicity). 
- Enter any of the PGS Catalog Identifiers if needed. 
- Select to run PGx profiling or not. 
- Click **"Submit"** to proceed.  

<br>

### 2. View Results
The submitted job and job status are presented in **"Analysis Jobs"** on the homepage. \
Once the job is completed, click the job name to view all results. \
Results are organized into 6 sections:
- Sample & Variant Summary: This section summarizes the sample ID, gender, number of variants retained, and clinical P/LP variants (both ClinVar and ACMG P/LP). 
- Screen-positive Results: This section presents the potential monogenic rare diseases identified in the newborn. These diseases are caused by both ClinVar and ACMG P/LP variants. 
- Disease Carrier Status: This section presents the potential disease carriers. These diseases are caused by both ClinVar and ACMG P/LP variants, which are heterozygous states and associated with recessive diseases. 
- Pharmacogenomics: This section provides a report summarizing all identified pharmacogenomics diplotypes in the newborns via PharmCat. Detailed explanations and recommendations are illustrated in the structured report. 
- Polygenic Risk Scores: This section 
- Variant Details: This section presents all details for the filtered variants.

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
<img width="558" height="722" alt="image" src="https://github.com/user-attachments/assets/3ff368af-fdb3-4810-a992-102b0c6311e3" />

<br>
<br>

### 4. Export Data  
- Results in each section could be downloaded separately.
    - <img width="477" height="66" alt="image" src="https://github.com/user-attachments/assets/de80c455-eddf-4a33-84ec-09db955c08cd" />
- All results could be wrapped up and downloaded as a zipped file by clicking the download button in the bottom-right corner.
    - <img width="437" height="69" alt="image" src="https://github.com/user-attachments/assets/08035174-5529-4e9e-b81b-d47ce5eaea93" />