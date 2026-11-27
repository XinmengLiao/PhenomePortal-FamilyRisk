# 👨‍👩‍👧‍👦 [FamilyRisk](https://www.phenomeportal.org/familyrisk)

### A Free-Access, Integrated Web Server for Early-Life Disease Screening and Clinical Reporting
FamilyRisk is a free-access, integrated web server for early-life disease screening and clinical reporting. It combines a curated screening knowledge base with genomic analysis pipelines for newborn and carrier screening at the individual, family and population level. The analyses are organised into a **Clinical Mode**, for monogenic disease screening and clinical reporting, and a **Research Mode**, for population-level newborn screening and exploratory analyses of complex traits and drug response. FamilyRisk is available through the web server and a REST API. No specific account or registration is needed. User privacy and personal data are secured and protected.

FamilyRisk also features the **FamilyRiskAI**, an intelligent chatbox assist clinical geneticists, genetic counselors and researchers to quickly explore the screening knowledge base, build customized screening lists and set optional filtering options, pinpoint variants of interest, and summarize the screening findings. 

📍 **Access FamilyRisk:** https://www.phenomeportal.org/familyrisk

<br>
<div align="center">
<img width="75%" alt="image" src="https://github.com/user-attachments/assets/52e08b1a-0de4-4376-9131-389851b94881" />
</div>
<br>


## Contents

- [Key Features](#key-features)
  - [Screening Knowledge Base](#screening-knowledge-base)
  - [Clinical Mode](#clinical-mode)
  - [Research Mode](#research-mode)
  - [Programmatic Access & AI Assistant](#programmatic-access-ai-assistant)
- [Tutorials](#tutorials)


## 🌟 Key Features

FamilyRisk organises its downstream analyses into two modes:

| Mode | Purpose | Analyses |
|:-----|:--------|:---------|
| 🩺 **Clinical Mode** | Monogenic disease screening and clinical reporting | Versatile genomic analysis · Segregation analysis · Clinical reporting |
| 🔬 **Research Mode** | Population-level newborn screening to understand regional rare disease characteristics, and exploratory analyses for complex traits and drug response | Cohort-level newborn screening · Polygenic risk evaluation · Pharmacogenomics profiling |

## Screening Knowledge Base

#### 1️⃣ Comprehensive Screening Database

FamilyRisk integrates curated knowledge on screening strategies, gene–disease associations, and clinical outcomes to support screening study design and clinical decision-making.

The database aggregates evidence from seven established newborn sequencing initiatives (BabySeq, BabyScreen+, BabyDetect, GUARDIAN, EarlyCheck, Genomics 101, NBScreening), ACMG Secondary Findings, and ACMG Expanded Carrier Screening.

<div align="center">
<img width="75%" alt="FamilyRisk database overview" src="https://github.com/user-attachments/assets/c2ff8e49-444d-4c5d-8d96-81c9864a36f1" />
</div>

#### 2️⃣ Flexible Panel Customization

The **Screening List Generator** embedded in the database allows users to:

- Select preset screening lists
- Add or remove gene–disease associations
- Create customized panels tailored to specific clinical or research needs

This feature supports adaptable implementation across diverse screening scenarios.

## Clinical Mode

#### 3️⃣ Versatile Genomic Analysis

FamilyRisk streamlines genomic analysis at the individual and family level. It supports interpretation of small variants (SNVs/indels), copy number variants, and structural variants. Compatible with both GRCh37 and GRCh38 reference genomes.

#### 4️⃣ Segregation Analysis

FamilyRisk provides integrated segregation analysis to evaluate inheritance patterns and assess transmission of pathogenic variants within families.

#### 5️⃣ Flexible Clinical Reporting

FamilyRisk separates results into structured sub-panels for clarity and usability in the user-friendly webpage.

Users can:

- Select specific variants or diseases
- Explore the inherited mode and visualize the segregation pattern
- Generate customizable clinical PDF reports
- Add professional comments, clinical interpretations, and recommendations

This facilitates translation from genomic data to actionable clinical documentation.

## Research Mode

#### 6️⃣ Population-level Newborn Screening

FamilyRisk performs newborn screening on a whole cohort from a multi-sample VCF and a sample sheet. It identifies regionally prevalent monogenic rare diseases and disease carrier status, and calculates cohort-level allele frequencies, helping users understand the rare disease characteristics of a local population.

#### 7️⃣ Polygenic Risk Evaluation

FamilyRisk enables polygenic score (PGS) evaluation for complex diseases and common traits. It also compares individual's PGS against the reference populations to estimate genetic susceptibility.

#### 8️⃣ Pharmacogenomics Profiling

FamilyRisk phases Pharmacogenomic (PGx) diplotype to provide Evidence-based therapeutic recommendations.

## Programmatic Access & AI Assistant

#### 9️⃣ REST API

FamilyRisk is also available as a REST API at https://webserver-api.cloudxomics.com/api/familyrisk.

- **One service, two entry points:** the web page is itself a client of the public API, so everything you can do on the page can also be scripted.
- **No account required:** each visitor or API client holds a random access token. Only its SHA-256 digest is stored, and every analysis is visible only to the token that created it.
- **Validated before running:** each submission is checked against the analysis's published contract (input files, parameters, defaults and allowed values), so invalid requests are refused immediately instead of failing hours later.
- **Direct and private uploads:** VCF files go straight from the client to private object storage through signed upload links. The web server keeps no copy of files or results.
- **Reproducible results:** each analysis runs in an isolated container from a pinned, versioned pipeline release, so every result can be traced to the exact pipeline and reference data that produced it.
- **Secure by design:** TLS 1.3 encryption, per-address request limits, short-lived download links (one hour), and no request bodies, tokens or cookies in logs.

#### 🔟 FamilyRiskAI Chatbox

**FamilyRiskAI** is an AI assistant for the clinical geneticists, genetic counselors and researchers who use FamilyRisk.

- **Explore the knowledge base:** find genes, diseases and preset panels in the curated screening list (1,556 genes, 3,112 gene–disease associations with GenCC validity), and compare the 18 preset panels.
- **Build custom panels:** assemble custom gene panels using only genes from the curated list, then download them as a gene list with one click.
- **Configure analyses:** explain each analysis and draft its configuration, checked against the live pipeline parameters. Proposals arrive as cards the user can apply to the analysis form.
- **Explain results:** summarise and explain the results of your finished analyses and of the public demos. It reads result rows only and never sees sample names or raw VCF fields.
- **Grounded and transparent:** each key fact names its source, for example *(curated list)*. Every step is streamed live with its purpose, and the full trace can be exported.
- **Read-only and privacy-preserving:** the assistant works only through read-only tools and never submits, changes or deletes anything. Conversations are not stored, and the service log keeps no questions, answers or results.


### ❗Note
FamilyRisk and FamilyRiskAI do not provide a diagnosis, prognosis or treatment recommendation. All screening findings need confirmation and interpretation by a genetic counselor or clinical geneticist.


## Tutorials

The full tutorial documentation can be found here:

- 📖 [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Home/)
- 📚 [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki)

Each function of FamilyRisk is documented with step-by-step Markdown guides:

| Function | Description | Documentation | Wiki |
|:---------|:------------|:-------------:|:----:|
| **Database** | Integrated knowledge of carrier and newborn screening projects, gene panels, and clinical outcomes. | [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Database/) | [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki/Database) |
| **Newborn Screening-Individual** | Individual-level newborn screening. Identify potential monogenic rare diseases and disease carrier status. | [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Newborn-Screening-%E2%80%90-Individual/) | [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki/Newborn-Screening-%E2%80%90-Individual) |
| **Newborn Screening-Family** | Family-level newborn screening. Identify potential monogenic rare diseases and disease carrier status in the newborns. Segregation analysis to visualize variant inheritance patterns. | [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Newborn-Screening-%E2%80%90-Family/) | [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki/Newborn-Screening-%E2%80%90-Family) |
| **Newborn Screening-Population** | Cohort-level newborn screening. Identify regional prevalent monogenic rare diseases and disease carrier status. Calculate cohort-level allele frequency. | [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Newborn-Screening-%E2%80%90-Population/) | [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki/Newborn-Screening-%E2%80%90-Population) |
| **Newborn Screening-CNV/SV** | Annotate individual-level copy number variants and structural variants in the newborn. | [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Newborn-Screening-%E2%80%90-CNV-%26-SV/) | [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki/Newborn-Screening-%E2%80%90-CNV-&-SV) |
| **Carrier Screening-Individual** | Individual-level carrier screening. Identify potential recessive diseases that might inherit by the offspring. | [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Carrier-Screening-%E2%80%90-Individual/) | [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki/Carrier-Screening-%E2%80%90-Individual) |
| **Carrier Screening-Family** | Couple-level carrier screening. Identify potential recessive diseases that might inherit by the offspring. | [Documentation](https://xinmengliao.github.io/PhenomePortal-FamilyRisk/Carrier-Screening-%E2%80%90-Family/) | [Wiki](https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/wiki/Carrier-Screening-%E2%80%90-Family) |
