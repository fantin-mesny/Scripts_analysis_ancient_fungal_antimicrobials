import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm

##parse taxonomy
metadata=pd.read_csv('data/genomeTaxonomy.csv').set_index('JGI ID')

##parse OG gene counts
og=pd.read_csv('data/Orthogroups.GeneCount.tsv',sep='\t').set_index('Orthogroup').drop(columns='Total').T


##classify each orthogroup by conservation level
output=[]
for o in og.columns:
	haveOG=list(og[og[o]>0].index)
	if len(haveOG)>1:
		for c in metadata.columns:
			occurrence=pd.DataFrame(metadata[metadata.index.isin(haveOG)][c].value_counts())
			if len(occurrence[occurrence['count']>0])==1:
				group=occurrence[occurrence['count']>0].index[0]
				output.append({'Orthogroup':o,'Specificity':c+'-specific','Taxon':group,'Organisms':';'.join(haveOG)})
				break
			elif c=='Phylum' and len(occurrence[occurrence['count']>0])>1:
				group=occurrence[occurrence['count']>0].index[0]
				output.append({'Orthogroup':o,'Specificity':'Multi-phyla','Taxon':group,'Organisms':';'.join(haveOG)})
output=pd.DataFrame(output)
OGspec=output.set_index('Orthogroup')['Specificity'].to_dict()

##parse prediction of antimicrobial activity on representative sequences of all orthogroups
annot=pd.read_csv('data/effectors_all_annotation.csv').set_index('index').drop(columns='Unnamed: 0')
annot['Specificity']=annot.index.map(OGspec)
annot=annot[(annot['isCAZyme']==0) & (annot['isTM']==0)].dropna() ## exclude transmembrane and CAZymes

##link AMAPEC prediction to conservation level
groups=list(set(output['Specificity']))
outp=[]
for spec in groups:
	annotspec=annot[(annot['Specificity']==spec)]
	outp.append({'Specificity':spec,'Number of antimicrobials':len(annotspec[annotspec['Prediction']=='Antimicrobial']),'Number of non-antimicrobials':len(annotspec[annotspec['Prediction']=='Non-antimicrobial'])})

outp=pd.DataFrame(outp).set_index('Specificity')
outp['Sum']=outp.sum(axis=1)
outp['pcAM']=(outp['Number of antimicrobials']/outp['Sum'])*100
order=['Multi-phyla','Phylum-specific','Class-specific','Order-specific','Family-specific','Genus-specific']
outp['order']=outp.index.map({val:order.index(val) for val in order})
outp=outp.sort_values(by='order').drop(columns='order')
print(outp)

fig,ax=plt.subplots(1,2,figsize=(8,2),sharey=True, gridspec_kw={'width_ratios': [3, 1]})
sns.barplot(ax=ax[0],data=outp.reset_index(drop=False),y='Specificity',x='Sum',color='#999999ff',order=order)
sns.barplot(ax=ax[0],data=outp.reset_index(drop=False),y='Specificity',x='Number of antimicrobials',color='#ff5555ff',order=order)
sns.scatterplot(ax=ax[1],data=outp.reset_index(drop=False),y='Specificity',x='pcAM',color='black')
plt.savefig('proportionOfAntimicrobials.pdf')
plt.close()




## Cochrane-Armitage test

outp2=outp.drop(columns=['Sum','pcAM'])
contingency_table = sm.stats.Table(outp2)
result = contingency_table.test_ordinal_association(
    col_scores=np.array([1,0]),
    row_scores=np.array([5,4,3,2,1,0])
)
print(f"Cochran-Armitage test statistic: {result.statistic:.4f}")
print(f"P-value: {result.pvalue}")

