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

# Multiple DNMT
dnmt1 = pd.concat([pd.DataFrame(rna.loc['DNMT1'].values, columns=['Expression']), pd.DataFrame(['DNMT1']*168, columns=['DNMT']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
dnmt3a = pd.concat([pd.DataFrame(rna.loc['DNMT3A'].values, columns=['Expression']), pd.DataFrame(['DNMT3A']*168, columns=['DNMT']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
dnmt3b = pd.concat([pd.DataFrame(rna.loc['DNMT3B'].values, columns=['Expression']), pd.DataFrame(['DNMT3B']*168, columns=['DNMT']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
dnmt = pd.concat([dnmt1, dnmt3a, dnmt3b], ignore_index=True)
sns.violinplot(data=dnmt, x='DNMT', y='Expression', hue='TN', palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
sns.stripplot(data=dnmt, x='DNMT', y='Expression', hue='TN', color=".3", size=3, dodge=True)

# Multiple GATA
gata5 = pd.concat([pd.DataFrame(rna.loc['GATA5'].values, columns=['Expression']), pd.DataFrame(['GATA5']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
gata6 = pd.concat([pd.DataFrame(rna.loc['GATA6'].values, columns=['Expression']), pd.DataFrame(['GATA6']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
gata4 = pd.concat([pd.DataFrame(rna.loc['GATA4'].values, columns=['Expression']), pd.DataFrame(['GATA4']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
gata2 = pd.concat([pd.DataFrame(rna.loc['GATA2'].values, columns=['Expression']), pd.DataFrame(['GATA2']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)

gata = pd.concat([gata5, gata6, gata4, gata2], ignore_index=True)
sns.violinplot(data=gata, x='GATA', y='Expression', hue='TN', palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
sns.stripplot(data=gata, x='GATA', y='Expression', hue='TN', color=".3", size=3, dodge=True)

# 5 core EMT TFs
snai1 = pd.concat([pd.DataFrame(rna.loc['SNAI1'].values, columns=['Expression']), pd.DataFrame(['SNAI1']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
twist1 = pd.concat([pd.DataFrame(rna.loc['TWIST1'].values, columns=['Expression']), pd.DataFrame(['TWIST1']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
snai2 = pd.concat([pd.DataFrame(rna.loc['SNAI2'].values, columns=['Expression']), pd.DataFrame(['SNAI2']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
zeb1 = pd.concat([pd.DataFrame(rna.loc['ZEB1'].values, columns=['Expression']), pd.DataFrame(['ZEB1']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
zeb2 = pd.concat([pd.DataFrame(rna.loc['ZEB2'].values, columns=['Expression']), pd.DataFrame(['ZEB2']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)

emt_tfs = pd.concat([snai1, twist1, snai2, zeb1, zeb2], ignore_index=True)
sns.violinplot(data=emt_tfs, x='EMT_TFs', y='Expression', hue='TN', palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
sns.stripplot(data=emt_tfs, x='EMT_TFs', y='Expression', hue='TN', color=".3", size=3, dodge=True)



