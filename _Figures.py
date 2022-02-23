####################################################################################
# Boxplot : Mean methylation (Tumor vs Normal)
import seaborn as sns

df = dmr_t.obs[['mean_met_normal', 'mean_met_tumor']]
df.rename(columns={'mean_met_normal': 'Normal', 'mean_met_tumor': 'Tumor'}, inplace=True)

sns.boxplot(data=df, palette={'Normal': '#00008b','Tumor': '#8b0000'})
sns.swarmplot(data=df, color="k").set(ylabel='Mean Methylation Level')

new_df = df['Tumor'] - df['Normal']
sns.boxplot(data=new_df, color='lightsalmon')
sns.swarmplot(data=new_df, color='k').set(ylabel='(Tumor - Normal) Mean methylation')

stats.ttest_rel(df['Normal'], df['Tumor'])

####################################################################################
