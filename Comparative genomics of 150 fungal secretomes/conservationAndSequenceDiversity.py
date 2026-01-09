import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

## parse orthology prediction data
o=pd.read_csv('data/Orthogroups.tsv',sep='\t').set_index('Orthogroup').fillna('') #orthogroup composition in proteins

opa=pd.read_csv('data/Orthogroups.GeneCount.tsv',sep='\t').set_index('Orthogroup') 
opa[opa>1]=1  #orthogroup presence/absence across 150 secretomes
opa['Total']=opa.sum(axis=1)
genomeCountPerOG=opa['Total'].to_dict()
opa=opa.sort_values(by='Total',ascending=False) #sorting by decreasing conservation
opa=opa.drop(columns=['Total'])
order=list(pd.read_csv('data/orderTreeFigure.txt',header=None)[0])[::-1] # order of fungal strains in the phylogenetic tree
opa=opa[order]

ogc=pd.read_csv('data/Orthogroups.GeneCount.tsv',sep='\t').set_index('Orthogroup') #orthogroup gene counts per secretome
ogc['Total']=ogc.sum(axis=1)

## parse effector annotation data
effectors_top600=pd.read_csv('data/top600_family_annotation.csv').set_index('Unnamed: 0')
effectors_top600=effectors_top600.sort_values(by='rank',ascending=True)
effectors_top600=effectors_top600[(effectors_top600['isTM']==0)] #exclude transmembrane proteins

## parse sequence lengths
seqlen=pd.read_csv('data/allSecretomes.length.tsv',sep='\t')
seqlen['#name']=seqlen['#name'].str.split(' ').str[0]
seqlen1=seqlen[(seqlen['#name'].str.contains('|',regex=False))]
seqlen2=seqlen[~(seqlen['#name'].str.contains('|',regex=False))]
seqlen1['#name']=seqlen1['#name'].str.split('|').str[0]+'|'+seqlen1['#name'].str.split('|').str[1]+'|'+seqlen1['#name'].str.split('|').str[2]
seqlen=pd.concat([seqlen1,seqlen2])
seqlen=seqlen.set_index('#name')
lenPerOG={}
for og in o.index:
	lenPerOG[og]=[]
	for org in o.columns:
		for prot in o.loc[og,org].split(', '):
			if prot!='':
				if '|' in prot:
					prot2='|'.join([prot.split('|')[0],prot.split('|')[1],prot.split('|')[2]])
				else:
					prot2=prot
				try:
					lenPerOG[og].append((prot2,seqlen.loc[prot2,'length']))
				except:
					print(prot2)
			else:
				pass
			
length=[]
for og in opa.index:
	for value in lenPerOG[og]:
		length.append({'og':og,
				'length':value[1],
				'prot':value[0],
				})
length=pd.DataFrame(length)



## plot presence/absence heatmap

top150=list(effectors_top600.index[:150])
effectors_top150=effectors_top600[effectors_top600.index.isin(top150)]

fig,ax=plt.subplots(4,1, gridspec_kw={'height_ratios': [4, 30, 1,1]})
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.01)

df2plot=opa[opa.index.isin(top150)].T
print(df2plot)
sns.heatmap(df2plot,cmap='Greys',ax=ax[1],cbar=False)
ax[1].set_xticks([])

sns.heatmap(effectors_top600[['isCAZyme']].T[list(df2plot.columns)],cmap='Greens',ax=ax[2],cbar=False)
ax[2].set_xticks([])
ax[2].set_xlabel('')

pred_top600=effectors_top600[['Prediction']]
pred_top600['Prediction']=pred_top600['Prediction'].map({'Non-antimicrobial':0,'Antimicrobial':1})
sns.heatmap(pred_top600.T[list(df2plot.columns)],cmap='Reds',ax=ax[3],cbar=False)

ax[3].set_xticks([])
ax[3].set_xlabel('')
length2plot=length[length['og'].isin(df2plot.columns)]
length2plot['rank']=length2plot['og'].map(effectors_top150['rank'].to_dict())
sns.barplot(data=length2plot.sort_values(by='rank',ascending=True),x='og',y='length',color='grey',ax=ax[0],errwidth=0.2,errorbar='sd')
ax[0].set_xticks([])
plt.savefig('presenceAbsenceHeatmap.pdf')
plt.close()



## plot diversity
effectors_top600=effectors_top600[(effectors_top600['isCAZyme']==0) & (effectors_top600['isTM']==0)] #exclude both transmembrane proteins and CAZymes
effectors_top600['genomeCount']=effectors_top600.index.map(genomeCountPerOG)

#parse blast
blast=pd.read_csv('data/publishedAM_on_dataset.blast.tsv',sep='\t',header=None)
blast=blast[blast[10]<0.001]
hitToOG={}
for ind in blast.index:
	for col in o.columns:
		if len(o[o[col].str.contains(blast.loc[ind,1],regex=False)])>0:
			hitToOG[blast.loc[ind,1]]=o[o[col].str.contains(blast.loc[ind,1],regex=False)].index[0]
blast['OG']=blast[1].map(hitToOG)
effectors_top600['knownHomologs']=effectors_top600.index.isin(blast['OG'])
#print(effectors_top600[effectors_top600['knownHomologs']]) #this shows the 4 orthogroups that include antimicrobials included in the AMAPEC dataset
containKnownAM=['OG0000020','OG0000128'] #these 2 include proteins which are in the AMAPEC dataset
containHomologsOfKnownAM=['OG0000375','OG0000673'] #these 2 include homologs of the some antimicrobials in the AMAPEC dataset

#parse mercat
shannon=pd.read_csv('data/Mercat_ShannonIndexes.csv').set_index('Orthogroup')
effectors_top600=effectors_top600.merge(shannon,left_index=True,right_index=True,how='left')
effectors_top600=effectors_top600.dropna()

fit=np.polyfit(np.log(effectors_top600['genomeCount']),effectors_top600['shannon'],1)
print(fit)
fit_am=np.polyfit(np.log(effectors_top600[effectors_top600['Prediction']=='Antimicrobial']['genomeCount']),effectors_top600[effectors_top600['Prediction']=='Antimicrobial']['shannon'],1)
fit_nam=np.polyfit(np.log(effectors_top600[effectors_top600['Prediction']=='Non-antimicrobial']['genomeCount']),effectors_top600[effectors_top600['Prediction']=='Non-antimicrobial']['shannon'],1)
effectors_top600['fit']=fit[0]*np.log(effectors_top600['genomeCount'])+fit[1]
effectors_top600['fit_am']=fit_am[0]*np.log(effectors_top600['genomeCount'])+fit_am[1]
effectors_top600['fit_nam']=fit_nam[0]*np.log(effectors_top600['genomeCount'])+fit_nam[1]


fig,ax=plt.subplots(1,1)
sns.scatterplot(ax=ax,x='genomeCount',y='shannon',data=effectors_top600, hue='Prediction',palette={'Antimicrobial':'#ff5555', 'Non-antimicrobial':'#999999'})
sns.scatterplot(ax=ax,x='genomeCount',y='shannon',data=effectors_top600[effectors_top600.index.isin(containKnownAM)], hue='Prediction',palette={'Antimicrobial':'#ff5555', 'Non-antimicrobial':'#999999'},edgecolor='black',linewidth=1.5,legend=None)
sns.lineplot(data=effectors_top600,y='fit',x='genomeCount',ax=ax,color='black')
ax.set_xlabel('Number of genomes containing gene family')
ax.set_ylabel('Gene family sequence diversity (Shannon index)')
ax.invert_xaxis()
plt.savefig('sequenceDiversityScatterplot.pdf')
plt.close()