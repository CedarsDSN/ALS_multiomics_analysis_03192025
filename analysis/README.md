**Analysis**
============

Data was requested from AnswerALS on March 27th, 2025. Data user agreement was received on 9th April, 2025. Access to data was granted on 28th April, 2025. On the website, the full clinical metadata was downloaded. This folder contains all clinical data that is available for the ALS patients, Non-ALS MND patients, carriers, and healthy controls. There is no readme available for the files. There is a file "Clinical_Dictionary_v3_2024.csv" which has the columns for each CSV file in the clinical folder. There are a few files outside the folder for which there is no readme available. The first step is to categorize the patients according to how the disease is progressing. The most natural way is to use the ALS progression rate if there is a column in the metadata files. The file "aals_dataportal_datatable.csv" has a column "ALSFRS_R_PROGRESSION_SLOPE", but this column has a lot of missing data. Another way to calculate the ALS progression rate. It is calculated between the maximum ALSFRS-R score (48) and the ALSFRS-R score at diagnosis, divided by disease duration (in months), as shown in two papers ([1](https://www.nature.com/articles/s41467-024-49954-y), [2](https://www.tandfonline.com/doi/full/10.1080/21678421.2023.2297795#d1e419)) among many others. 

The columns needed to perform the above calculations were in three different tables - "aals_dataportal_datatable.csv", "Demographics.csv", and "ALSFRS-R.csv". I wrote a script and uploaded it to the scripts folder, "calculate_progression.py," that can read all tables and add a column that adds the type of progression - fast, slow, or intermediate. It also adds more columns for sanity check.

It merges visit-wise ALSFRS-R scores, demographic information (including date of birth), and baseline clinical metadata to reconstruct a full disease timeline. The script identifies each patient’s latest ALSFRS-R visit, computes disease duration (from symptom onset to latest visit), and calculates the monthly rate of functional decline.

Patients are then classified into slow, intermediate, or fast progressors based on thresholds from the studies shared above.

Additional features:

* Filters out healthy controls and patients with unreliable follow-up (e.g., only baseline visit, negative or too-short disease durations).

* Performs a consistency check between the recorded ALSFRS-R latest values and those computed from visit-level data.

* Flags data quality issues with a QC_Flag column for transparency.

* Outputs a cleaned CSV file with progression metrics and classification.

* Includes optional plotting code to visualize class distributions.

The output file has been uploaded to analysis folder as "ALS_Progression_Rates_Full.csv". The plot showing distributions of the class is uploaded as "Progressor_Class_Counts.png". 

RNA-seq Differential Expression Summary Across ALS Progression Groups
Sample Distribution (Progressor Classes):

•	Slow: 163 samples

•	Intermediate: 141 samples

•	Fast: 116 samples

Gene Expression Overview:

•	Mean normalized expression: 1,080.7 counts per gene

•	Genes retained after filtering: 27,818 genes across 420 samples

1. Fast vs Slow
   
•	Min LFC: –3.33

•	Median LFC: –0.0007

•	Mean LFC: 0.0037

•	Max LFC: 2.20

•	Significant DEGs (padj < 0.05): 0

2. Intermediate vs Slow
   
•	Min LFC: –3.08

•	Median LFC: –0.0036

•	Mean LFC: –0.0160

•	Max LFC: 2.06

•	Significant DEGs (padj < 0.05): 2

3. Fast vs Intermediate
   
•	Min LFC: –2.70

•	Median LFC: 0.0036

•	Mean LFC: 0.0198

•	Max LFC: 2.67

•	Significant DEGs (padj < 0.05): 17

Key Insight:

The data shows very small fold changes across most genes, with the majority clustered tightly around 0. This suggests minimal transcriptomic differences between progression classes at the gene level.
I plan to categorize the data into fast and slow progressors only since the intermediate group might be masking clearer differences.

