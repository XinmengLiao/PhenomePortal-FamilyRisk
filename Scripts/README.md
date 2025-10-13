## Function selection
  a) Single sample screening and reporting \
  b) Batch samples screening and exploration \
  c) Family-based batch samples screening and inheritance identification

## Input options:
  a) VCF file (same as XOmics and VarXOmics) \
  b) PED file (required for family-based batch sample screening)
  c) SampleID list (required for single and batch sample screening)
  
## Filtrations on the home page:
### Sequencing quality filtrations
   a) Reference genome version: GRCH37 or GRCH38 \
   b) 'PASS only' or not
   c) Gender: Male, Female, Unknown

### Variant filtrations (Same as VarXOmics)

   **Clinical Pathogenicity (multi-Selections)** \
     a) ClinVar Pathogenicity (6 options): Pathogenic, Likely pathogenic, Uncertain Significance, Conflicting classifications, Benign, Likely benign \
     b) ACMG Predicted Classification (5 options):  Pathogenic, Likely pathogenic, Uncertain Significance, Benign, Likely benign
     c) ACMG Pathogenicity: default display: 0.564

  **Allele Frequency (float)** \
  a) ClinVar variant: default display: 1 \
  b) Predicted variant: default display: 0.05 
      
  **Variant Deleteriouseness** \
  a) PolyPhen Score: default display: 0.05 \
  b) BayesDel Score (addAF): default display: 0.0692655 \
  c) BayesDel Score (noAF): default display: -0.0570105
      
  **Missense Impact** \
  a) AlphaMissense Classification (multi-selection, 3 options): Likely pathogenic, Uncertain Significance, Likely benign \
  b) AlphaMissense Score: default display: default display: 0.564 \
  c) REVEL Score: default display: default display: 0.75 \
  d) Ada Score: default display: 0.6 \
  e) Rf score: default display: 0.6
  
  **Splicing Impact** \
  a) Acceptor Gain :default display: 0.5 \
  b) Acceptor Loss: default display: 0.5 \
  c) Donor Gain: default display: 0.5 \
  d) Donor Loss: default display: 0.5 

### Screening list selection
  a) BabySesq \
  b) BabyDetect \
  c) BabyScreen+ \
  d) GUARDIAN \
  e) NC NEXUS \
  f) NBScreening \
  g) Customized list (user should upload their own list in fixed format) 

Screening lists could be found in the local server `/path_to_NewbornRisk/Datasets/genelists/`
   
