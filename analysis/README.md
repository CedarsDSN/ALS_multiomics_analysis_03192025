**Analysis Strategy**
=====================

## Fast, intermediate, and slow progressors

### Categorizing patients into fast, intermediate, and slow progressors

#### I. Analysis Script:

#### 1. Use ALSFRS-R score and disease onset to derive the ALS progression and categorize patients

Required Packages/Software:

 - pandas package in python

**Scripts and Inputs:**
- **[calculate_progression.py](/scripts/calculate_progression.py)**
  - **Inputs:**
    - Clinical tables ALSFRS-R ([bladder_cancer_RNAseq_raw_counts.csv](analysis/categorization/inputs/ALSFRS-R.csv)), [Demographics.csv](analysis/categorization/inputs/Demographics.csv), [aals_dataportal_datatable.csv](analysis/categorization/inputs/aals_dataportal_datatable.csv)
  - **Description:** Generates a table [ALS_Progression_Rates_Full.csv](analysis/categorization/fast_slow_intermediate/ALS_Progression_Rates_Full.csv) which contains a column Progressor_Class which categorizes patients into fast, slow and intermediate progressors and a plot showing the distribution of patients into the three classes.
