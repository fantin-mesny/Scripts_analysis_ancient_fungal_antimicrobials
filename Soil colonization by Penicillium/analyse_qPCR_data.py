import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Parse qPCR deltaCt
df=pd.read_csv('data/deltaCtvalues.csv',sep=';',na_values='n.a',decimal=',').fillna(0)

# Plot figure
fig,ax=plt.subplots(1,2,figsize=(8,5))
sns.boxplot(x='Timepoint',y='Values',hue='Fungus',ax=ax[0],data=df[df['Soil']=='Sterile soil'],fliersize=0,palette={'WT':'grey','Mutant':'#ff5555ff'})
sns.swarmplot(x='Timepoint',y='Values',hue='Fungus',ax=ax[0],data=df[df['Soil']=='Sterile soil'],dodge=True,size=2,color='black',legend=None)
ax[0].set_title('Sterile soil')
ax[0].set_ylabel('Soil colonization')
sns.boxplot(x='Timepoint',y='Values',hue='Fungus',ax=ax[1],data=df[df['Soil']=='Recolonized soil'],fliersize=0,legend=None,palette={'WT':'grey','Mutant':'#ff5555ff'})
sns.swarmplot(x='Timepoint',y='Values',hue='Fungus',ax=ax[1],data=df[df['Soil']=='Recolonized soil'],dodge=True,size=2,color='black',legend=None)
ax[1].set_title('Recolonized soil')
ax[1].set_ylabel(' ')
plt.savefig('figure.pdf')
plt.close()

## Export data for statistical analysis in R
#df['Category']=df['Timepoint'].astype(str)+'_'+df['Fungus']
#df[df['Soil']=='Recolonized soil'].to_csv('data_recol.csv')
#df[df['Soil']=='Sterile soil'].to_csv('data_sterile.csv')

###Stats calculated in R with the following commands:
##data<-read.csv('data4R_recol.csv')
##K<-kruskal.test(Values~Category,data=data)
##print(K)
##library(DescTools)
##library(multcompView)
##D<-DunnTest(Values~Category,data=data)
##Ddf<-data.frame(D[[1]])
##p<-c(Ddf$pval)
##names(p)<-rownames(Ddf)
##print(multcompLetters(p)$Letters)