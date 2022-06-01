# source activate wc300
# ipython --profile=wc300

import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
sns.set(font="Arial", font_scale=1.2, style='ticks')
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap4 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#191970", "#ffdab9", "#8B0000"])
%matplotlib
%autoindent

# Clinical information
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')
#clinic_info = pd.read_csv('/home/mhryan/Workspace/02.Projects/02.WC300/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')


