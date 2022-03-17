####################################################################################
# Boxplot : Mean methylation (Tumor vs Normal) ==> deprecated
import pandas as pd
import seaborn as sns

df = dmr_t.obs[['mean_met_normal', 'mean_met_tumor']]
df.rename(columns={'mean_met_normal': 'Normal', 'mean_met_tumor': 'Tumor'}, inplace=True)

sns.boxplot(data=df, palette={'Normal': '#00008b','Tumor': '#8b0000'})
sns.swarmplot(data=df, color="k").set(ylabel='Mean Methylation Level')

new_df = df['Tumor'] - df['Normal']
sns.boxplot(data=new_df, color='lightsalmon')
sns.swarmplot(data=new_df, color='k').set(ylabel='(Tumor - Normal) Mean methylation')

stats.ttest_rel(df['Normal'], df['Tumor'])

# Violinplot : Mean methylation (Tumor vs Normal) ==> recent version
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')
p = sns.violinplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], palette={'PercentMet_COV5_Normal':'navy', 'PercentMet_COV5_Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], color="black")
p.set_xticklabels(['Normal', 'Tumor'])
p.set_title('Average Methylation')

####################################################################################
# Violinplot : PMD counts (Tumor vs Normal)
df = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/PMD/PMD_counts.txt", sep='\t')
sns.violinplot(data=df, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_xticklabels(['Normal', 'Tumor'])
sns.stripplot(data=df, color="k").set_title('PMD counts')

####################################################################################
# Violinplot : DNMT (Tumor vs Normal)

dnmt1 = pd.concat([pd.DataFrame(rna.loc[['DNMT1']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(rna.loc[['DNMT1']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(rna.columns[84:])
sns.violinplot(data=dnmt1, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_title('DNMT1')

dnmt3a = pd.concat([pd.DataFrame(rna.loc[['DNMT3A']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(rna.loc[['DNMT3A']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(rna.columns[84:])
sns.violinplot(data=dnmt3a, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_title('DNMT3A')

dnmt3b = pd.concat([pd.DataFrame(rna.loc[['DNMT3B']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(rna.loc[['DNMT3B']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(rna.columns[84:])
sns.violinplot(data=dnmt3b, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_title('DNMT3B')


