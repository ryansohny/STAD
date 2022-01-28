#export PATH=/ruby/Users/uugi0620/phenomata/Anaconda3/bin:$PATH
#source activate cemitool

# source activate cemitool
library(CEMiTool)
setwd("/mnt/mone/Project/WC300/02.RNA-seq/02.Kallisto_Analysis/CEMiTool_tumor")
expr <- read.csv("STAD_SNUH_Tumor_leiden_vst.csv", row.names=1, header=T)

## Make annotation files (using python)
#import pandas as pd
#a = ['leiden_0', 'leiden_2', 'leiden_2', 'leiden_0', 'leiden_1', 'leiden_3', 'leiden_0', 'leiden_0', 'leiden_2', 'leiden_4', 'leiden_0', 'leiden_1', 'leiden_2', 'leiden_0', 'leiden_0', 'leiden_2', 'leiden_0', 'leiden_2', 'leiden_0', 'leiden_0', 'leiden_4', 'leiden_3', 'leiden_3', 'leiden_1', 'leiden_1', 'leiden_4', 'leiden_1', 'leiden_0', 'leiden_3', 'leiden_4', 'leiden_4', 'leiden_2', 'leiden_0', 'leiden_4', 'leiden_1', 'leiden_0', 'leiden_1', 'leiden_1', 'leiden_2', 'leiden_4', 'leiden_0', 'leiden_2', 'leiden_1', 'leiden_3', 'leiden_3', 'leiden_2', 'leiden_0', 'leiden_2', 'leiden_1', 'leiden_1', 'leiden_0', 'leiden_2', 'leiden_2', 'leiden_1', 'leiden_1', 'leiden_0', 'leiden_0', 'leiden_4', 'leiden_3', 'leiden_3', 'leiden_1', 'leiden_2', 'leiden_0', 'leiden_0', 'leiden_1', 'leiden_2', 'leiden_1', 'leiden_1', 'leiden_4', 'leiden_4', 'leiden_0', 'leiden_1', 'leiden_3', 'leiden_1', 'leiden_0', 'leiden_3', 'leiden_2', 'leiden_0', 'leiden_1', 'leiden_0', 'leiden_1', 'leiden_2', 'leiden_1', 'leiden_0']
#leiden_dict = {'leiden_4':'leiden_A', 'leiden_1': 'leiden_B', 'leiden_0': 'leiden_C', 'leiden_3': 'leiden_D', 'leiden_2': 'leiden_E'}
#new_a = list()
#for i in a:
#  new_a.append(leiden_dict[i])
#b = ['14002_T', '14006_T', '14007_T', '14009_T', '14020_T', '14021_T', '15001_T', '15002_T', '15005_T', '15006_T', '16002_T', '16010_T', '16015_T', '16018_T', '16021_T', '16022_T', '16025_T', '16031_T', '17014_T', '17017_T', '17024_T', '17027_T', '17032_T', '17039_T', '17040_T', '17041_T', '17043_T', '17044_T', '17049_T', '18001_T', '18002_T', '18003_T', '18004_T', '18005_T', '18006_T', '18007_T', '18008_T', '18009_T', '18010_T', '18011_T', '18013_T', '18014_T', '18015_T', '18016_T', '18017_T', '19001_T', '19004_T', '19005_T', '19006_T', '19007_T', '19008_T', '19010_T', '19014_T', '19021_T', '19025_T', '19028_T', '19032_T', '19033_T', '19035_T', '19036_T', '19038_T', '19040_T', '19041_T', '19055_T', '19081_T', '19082_T', '20001_T', '20002_T', '20004_T', '20018_T', '20037_T', '20054_T', '20068_T', '20080_T', '20082_T', '20084_T', '20087_T', '20089_T', '20092_T', '20094_T', '20095_T', '20102_T', '20114_T', '20115_T']
#c = pd.read_csv("/mnt/mone/Project/WC300/02.RNA-seq/02.Kallisto_Analysis/CEMiTool_tumor/2022_WC300_clinical_information_Xadded.csv", index_col=0)

#db = dict()
#for i in range(len(new_a)):
#  db["X" + b[i]] = new_a[i]

#c["DMR_leiden"] = pd.DataFrame.from_dict(db, orient="index", columns=["DMR_leiden"])
#c.to_csv("2022_WC300_clinical_information_Xadded_DMRleidenadded2.csv", sep=",")

#cem <- cemitool(expr)
#diagnostic_report(cem, force=TRUE) # 아래 save_plots 안해도 이걸로 해도 됨
#save_plots(cem, "all", force=TRUE, directory="./Plots") # force=TRUE 하면 이미 생성된 Plots 디렉토리 안에 덮어쓰기 함.
# 위처럼 하면, 아래와 같이, 일곱개의 PDF파일이 Plots 디렉토리와 함께 저장된다.
# sample_tree.pdf ==> 샘플 간 expression에서의 Hierarchical Clustering
# profile.pdf ==> *중요: 밝혀진 각각의 module에 속하는 유전자의 expression을 샘플별로 보여주는 plot
# mean_k.pdf ==> Mean connectivity plot. 아래의 beta R squared plot 과 비슷 (이건 나중에 다시 이해할 것)
# beta_r2.pdf ==> Soft threshold beta vs Scale-free topology model fit (R squared). automatically 설정된, B value가 어떻게 선택되었는 지 확인할 수 있다.
# mean_var.pdf ==> Expression file에서 각 gene들의 평균과 분산을 log scale로 보여준 plot. 만약 여기의 regression R squared value가 meaningful하다면 즉 크다면, cemitool() function에서 apply_vst argument를 TRUE로 설정하는 것을 통해 이 관계를 없애야 함. ==> 우리는 이미 vst 된 것을 사용하므로 상관X
# hist.pdf ==> Gene expression의 histogram을 보여준 것 같다. 보면, negative bionomial distribution처럼 생김 ==> RNA-seq 특성
# qq.pdf ==> Theoretical value와 real value를 비교하는 것을 통해서 우리가 가지고 있는 real value가 normal distribution화 되어있는지, 즉 어떤 distribution을 나타내고 있는지를 보여주기 위함 ==> RNA-seq의 일반적 특성, 그리고 우리의 real data가 NB distribution을 보이므로, 이 Q-Q plot은 y=x graph에서 멀어져있는 것이 확인됨.

# Adding sample annotation ==> data.frame

#annot <- read.table("/clinix1/Analysis/mongol/phenomata/04.GC_CC/Whole_new_2020/02.RNA-seq/04.HTSeq/Analysis/CEMiTool/Annot_STAD_n174_vst_modified.txt", sep="\t", header=T)

sample_annot <- read.csv("/mnt/mone/Project/WC300/02.RNA-seq/02.Kallisto_Analysis/CEMiTool_tumor/2022_WC300_clinical_information_Xadded_DMRleidenadded2.csv")
# head(cem@module) ## Module에 속하는 gene을 알아내는 방법
#go_gmt <- read_gmt("/clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/RNA-seq/Analysis_new/CEMiTool/c5.go.bp.v7.4.symbols.gmt")
hallmark_gmt <- read_gmt("/mnt/mone/Project/WC300/02.RNA-seq/02.Kallisto_Analysis/CEMiTool_tumor/h.all.v7.2.symbols.gmt")
#interactome_gmt <- read_gmt("/clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/RNA-seq/Analysis_new/CEMiTool/c2.cp.v7.4.symbols.gmt")
#int_string <- read.csv("/clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/RNA-seq/Analysis_new/CEMiTool/STRING_interaction.csv")
#int_string <- read.csv("/clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/RNA-seq/Analysis_new/CEMiTool/STRING_interaction_scorecutoff217.csv")
#int_string <- read.csv("/clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/RNA-seq/Analysis_new/CEMiTool/ARCHS4_Coexpression_interaction.csv")

########################## (i) ##########################
# network type ==> unsigned
# annotation ==> hallmark
#cem <- cemitool(expr=expr, annot=sample_annot, gmt=go_gmt, cor_method="spearman", network_type="unsigned", tom_type="signed", rank_method="mean", class_column="DMR_leiden", gsea_max_siz=1200, verbose=TRUE)
cem <- cemitool(expr=expr, annot=sample_annot, sample_name_column="Sample", gmt=hallmark_gmt, cor_method="spearman", network_type="unsigned", tom_type="signed", rank_method="mean",class_column="DMR_leiden", gsea_max_siz=1200, verbose=TRUE)

########################## (ii) ##########################
# network type ==> signed
cem <- cemitool(expr=expr, annot=sample_annot, gmt=hallmark_gmt, interactions=int_string, cor_method="spearman", network_type="signed", tom_type="signed", rank_method="mean",class_column="cluster", gsea_max_siz=1200, verbose=TRUE)

########################## (iii) ##########################
# network type ==> signed
# merge similar ==> FALSE
cem <- cemitool(expr=expr, annot=sample_annot, gmt=hallmark_gmt, cor_method="spearman", interactions=int_string, merge_similar=FALSE, network_type="signed", tom_type="signed", rank_method="mean",class_column="cluster", gsea_max_siz=2000, verbose=TRUE)

########################## (iv) ##########################
# network type ==> signed
# merge similar ==> TRUE
cem <- cemitool(expr=expr, annot=sample_annot, gmt=hallmark_gmt, cor_method="spearman", interactions=int_string, merge_similar=TRUE, network_type="signed", tom_type="signed", rank_method="mean",class_column="cluster", gsea_max_siz=2000, verbose=TRUE)

write_files(cem, directory="./Tables", force=TRUE)
# Create report as html file ################################## Recommended ##################################
library(ggrepel)
options(ggrepel.max.overlaps=Inf)
generate_report(cem, directory = "./Report", output_format="html_document", force=TRUE)

cem <- mod_gsea(cem)
show_plot(cem, "gsea")
dev.off()

cem <- mod_ora(cem, go_gmt)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
pdf("ORA.pdf")
plots
dev.off()


library(ggplot2)
interactions_data(cem) <- int_string
cem <- plot_interactions(cem)
plots <- show_plot(cem, "interaction")
pdf("interaction_ARCHS4.pdf")
plots
dev.off()



##################################################################################################################################################################
################################################### Raw counts ==> CEMiTool vst (2021-06-16) #####################################################################
##################################################################################################################################################################

library(CEMiTool)
setwd("/clinix1/Analysis/mongol/phenomata/04.GC_CC/Whole_new_2020/02.RNA-seq/04.HTSeq/Analysis/CEMiTool/RawCounts")
expr <- read.table("/clinix1/Analysis/mongol/phenomata/04.GC_CC/Whole_new_2020/02.RNA-seq/04.HTSeq/Analysis/CEMiTool/RawCounts/STAD_n174_STAR_counts_modified_onlyTumor.txt", row.names = 1, sep="\t", header=T)

sample_annot <- read.csv("/clinix1/Analysis/mongol/phenomata/04.GC_CC/Whole_new_2020/02.RNA-seq/04.HTSeq/Analysis/CEMiTool/DMR_cluster.csv")
hallmark_gmt <- read_gmt("/clinix1/Analysis/mongol/phenomata/04.GC_CC/Whole_new_2020/02.RNA-seq/04.HTSeq/Analysis/CEMiTool/GMT/h.all.v7.2.symbols.ensemblID.gmt")

cem <- cemitool(expr=expr, annot=sample_annot, gmt=hallmark_gmt, cor_method="pearson", apply_vst=TRUE, merge_similar=TRUE, network_type="signed", tom_type="signed", rank_method="mean",class_column="Methylation_subcluster_DMRabs15", gsea_max_size=2000, verbose=TRUE)

diagnostic_report(cem, force=TRUE)

write_files(cem, directory="./Tables", force=TRUE)
# Create report as html file ################################## Recommended ##################################
library(ggrepel)
options(ggrepel.max.overlaps=Inf)
generate_report(cem, directory = "./Report", output_format="html_document", force=TRUE)
