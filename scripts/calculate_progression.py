import pandas as pd

# Load the three datasets (you'll replace these with actual file paths or DataFrame inputs)
alsfrs_df = pd.read_csv('/Users/bhandaryp/Documents/DSN/ALS/metadata/clinical/ALSFRS_R.csv')
demo_df = pd.read_csv('/Users/bhandaryp/Documents/DSN/ALS/metadata/clinical/Demographics.csv')
portal_df = pd.read_csv('/Users/bhandaryp/Documents/DSN/ALS/metadata/aals_dataportal_datatable.csv')

# Step 1: Extract relevant columns
alsfrs_df = alsfrs_df[['Participant_ID', 'Visit_Name', 'Visit_Date', 'alsfrst']]
demo_df = demo_df[['Participant_ID', 'dob']]  # Use real date of birth
portal_df = portal_df[['Participant_ID', 'ALSFRS_R_Baseline_Value', 'ALSFRS_R_Latest_Value', 'AGE_AT_SYMPTOM_ONSET', 'SUBJECT_GROUP']]

# Step 2: Find Screening Visit and Latest Visit
screening_df = alsfrs_df[alsfrs_df['Visit_Name'] == 'ANSWER-ALS Screening Visit'][['Participant_ID', 'Visit_Date']]
screening_df = screening_df.rename(columns={'Visit_Date': 'Screening_Day'})

latest_visits = alsfrs_df.loc[alsfrs_df.groupby('Participant_ID')['Visit_Date'].idxmax()]
latest_visits = latest_visits[['Participant_ID', 'Visit_Date', 'alsfrst']]
latest_visits = latest_visits.rename(columns={'Visit_Date': 'Latest_Visit_Day', 'alsfrst': 'ALSFRS_Latest_Score'})

# Step 3: Merge everything
merged = portal_df.merge(screening_df, on='Participant_ID', how='left') \
                  .merge(latest_visits, on='Participant_ID', how='left') \
                  .merge(demo_df, on='Participant_ID', how='left')

# Step 4: Calculate symptom onset, latest visit, disease duration
merged['Symptom_Onset_Day'] = merged['dob'] + merged['AGE_AT_SYMPTOM_ONSET'] * 365.25
merged['Latest_Visit_Absolute_Day'] = merged['Screening_Day'] + merged['Latest_Visit_Day']
merged['Disease_Duration_days'] = merged['Latest_Visit_Absolute_Day'] - merged['Symptom_Onset_Day']
merged['Disease_Duration_months'] = merged['Disease_Duration_days'] / 30.44

# Step 5: QC flag based on data integrity
def qc_flag(row):
    if row['SUBJECT_GROUP'] == 'Healthy Control':
        return 'Healthy Control'
    elif pd.isnull(row['ALSFRS_Latest_Score']) or pd.isnull(row['Disease_Duration_months']):
        return 'Missing Data'
    elif row['Disease_Duration_months'] <= 0:
        return 'Negative Disease Duration'
    elif row['Disease_Duration_months'] < 1:
        return 'Too Short Disease Duration (<1 month)'
    else:
        return 'OK'

merged['QC_Flag'] = merged.apply(qc_flag, axis=1)

# Step 6: Progression Rate (only if QC passes)
merged['Progression_Rate'] = None
merged.loc[merged['QC_Flag'] == 'OK', 'Progression_Rate'] = (
    (48 - merged.loc[merged['QC_Flag'] == 'OK', 'ALSFRS_Latest_Score']) /
    merged.loc[merged['QC_Flag'] == 'OK', 'Disease_Duration_months']
)

# Step 7: Progressor Class (only if QC passes)
def classify(rate):
    if pd.isnull(rate):
        return 'Unknown'
    elif rate <= 0.36:
        return 'Slow'
    elif rate < 0.77:
        return 'Intermediate'
    else:
        return 'Fast'

merged['Progressor_Class'] = merged['Progression_Rate'].apply(classify)

# Step 8: ALSFRS-R Consistency Check
def check_mismatch(row):
    if pd.isnull(row['ALSFRS_Latest_Score']) or pd.isnull(row['ALSFRS_R_Latest_Value']):
        return 'Missing'
    elif row['ALSFRS_Latest_Score'] == row['ALSFRS_R_Latest_Value']:
        return 'Match'
    else:
        return 'Mismatch'

merged['ALSFRS_Consistency_Check'] = merged.apply(check_mismatch, axis=1)

# Step 9: Save full output
output_path = '/Users/bhandaryp/Documents/DSN/ALS/ALS_Progression_Rates_Full.csv'
merged[['Participant_ID',
        'ALSFRS_R_Latest_Value',
        'ALSFRS_Latest_Score',
        'ALSFRS_Consistency_Check',
        'Disease_Duration_months',
        'Progression_Rate',
        'Progressor_Class',
        'QC_Flag']].to_csv(output_path, index=False)

print(f"✅ Finished! Full results (with ALSFRS consistency check) saved to: {output_path}")

import matplotlib.pyplot as plt
import seaborn as sns

# Filter only valid rows (where progression rate is not null and QC passed)
plot_df = merged[merged['QC_Flag'] == 'OK'].copy()

# --- Bar plot of Progressor Class counts ---
plt.figure(figsize=(6, 4))
sns.countplot(data=plot_df, x='Progressor_Class', order=['Slow', 'Intermediate', 'Fast'])
plt.title('Number of Patients per Progression Class')
plt.ylabel('Count')
plt.xlabel('Progressor Class')
plt.tight_layout()
plt.savefig('/Users/bhandaryp/Documents/DSN/ALS/Progressor_Class_Counts.png')
plt.show()

