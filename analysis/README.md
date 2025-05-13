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
    - Clinical tables ALSFRS-R [ALSFRS-R.csv](analysis/categorization/inputs/ALSFRS-R.csv)), [Demographics.csv](analysis/categorization/inputs/Demographics.csv), [aals_dataportal_datatable.csv](analysis/categorization/inputs/aals_dataportal_datatable.csv)
  - **Description:** Generates a table [ALS_Progression_Rates_Full.csv](analysis/categorization/fast_slow_intermediate/ALS_Progression_Rates_Full.csv) which contains a column Progressor_Class which categorizes patients into fast, slow and intermediate progressors and a plot showing the distribution of patients into the three classes.
 
#### II. Process
**Deriving the ALSFRS-R Progression Rate**
   - Merged visit-wise ALSFRS-R scores, demographic information (including date of birth), and baseline clinical metadata to reconstruct a full disease timeline. The script identifies each patient’s latest ALSFRS-R visit, computes disease duration (from symptom onset to latest visit), and calculates the monthly rate of functional decline.
   - Filtered out healthy controls and patients with unreliable follow-up (e.g., only baseline visit, negative or too-short disease durations).
   - Performed a consistency check between the recorded ALSFRS-R latest values and those computed from visit-level data.
  
#### III. Results
The table [ALS_Progression_Rates_Full.csv](analysis/categorization/fast_slow_intermediate/ALS_Progression_Rates_Full.csv) has the categorization of patients into slow, fast and intermediate progressors.




   
