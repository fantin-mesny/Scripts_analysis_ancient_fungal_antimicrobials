#This script uses orthology prediction data (Orthogroups.tsv, Orthogroups.GeneCount.tsv, orderTreeFigure.txt) that can be found in the folder "Comparative genomics of 150 fungal secretomes/data/"

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def getFungusID(seqName):
	if seqName.startswith('jgi'):
		return seqName.split('|')[1]
	elif seqName.startswith('MG'):
		return 'MagnaMG8'
	elif seqName.startswith('VDAG'):
		return 'VerdaJR2'
	elif seqName.startswith('AB'):
		return 'Altbr1'
	else:
		print('ERROR:',seqName)
		return ''
		


df=pd.read_csv('data/prediction_on_phiBaseEffectors.csv')

# Exclude oomycetes in the dataset
df=df[~(df['Protein ID'].str.contains('Hyaloperonospora'))]
df=df[~(df['Protein ID'].str.contains('Phytophthora'))]
df=df[~(df['Protein ID'].str.contains('Pythium'))]
print(df)

# Link effector in Phi-base to orthogroup in 150-genome dataset 
blast=pd.read_csv('data/blastpResults_bestHit.tsv',sep='\t',header=None)
og=pd.read_csv('Orthogroups.tsv',sep='\t').set_index('Orthogroup').fillna('') # file can be found in the folder "Comparative genomics of 150 fungal secretomes/data/"
og_gc=pd.read_csv('Orthogroups.GeneCount.tsv',sep='\t').set_index('Orthogroup').drop(columns='Total') # file can be found in the folder "Comparative genomics of 150 fungal secretomes/data/"
og_gc[og_gc>1]=1
og_gc['Ngenomes']=og_gc.sum(axis=1)
Ngenomes=og_gc['Ngenomes'].to_dict()
for ind in blast.index:
	hit=blast.loc[ind,1]
	hit_sp=getFungusID(hit)
	hit_og=og[og[hit_sp].str.contains(hit.replace('|','\|'),regex='False')].index
	blast.loc[ind,'OG']=hit_og[0]
blast['Ngenomes']=blast['OG'].map(Ngenomes)
blast=blast.rename(columns={0:'Protein ID',1:'Best hit'})
blast=blast.merge(df,on='Protein ID',how='right').sort_values(by='Ngenomes',ascending=False)

# Prepare data for supplementary table and plot
df4scatterplot=blast[['Protein ID','Best hit',2,'OG','Ngenomes','Prediction','pLDDT','Probability of antimicrobial activity']]
df4scatterplot['LysM effector']=(df4scatterplot['Protein ID'].str.contains('lysm')) | (df4scatterplot['Protein ID'].str.contains('Lysm')) | (df4scatterplot['Protein ID'].str.contains('Blys')) | (df4scatterplot['Protein ID'].str.contains('Ecp6'))
df4scatterplot['Ngenomes']=df4scatterplot['Ngenomes'].fillna(0)
#df4scatterplot.to_csv('suptable_phibase.tsv',sep='\t') 

# Plot prediction/conservation for phi-base fungal effectors
g=sns.JointGrid(y='Ngenomes',x='Probability of antimicrobial activity',data=df4scatterplot[~(df4scatterplot['LysM effector'])],hue='Prediction',xlim=(-0.01,1.01),ylim=(-2,152))
sns.scatterplot(y='Ngenomes',x='Probability of antimicrobial activity', style='LysM effector',hue='Prediction',data=df4scatterplot.sort_values(by='LysM effector'),palette={'Antimicrobial':'#ff5555','Non-antimicrobial':'#b3b3b3'},ax=g.ax_joint,s=75)
g.plot_marginals(sns.kdeplot, fill=True,palette={'Antimicrobial':'#ff5555','Non-antimicrobial':'#b3b3b3'})
plt.savefig('jointplot_github.pdf')
plt.close()

# Plot conservation heatmap
og_gc=og_gc[og_gc.index.isin(blast['OG'])].sort_values(by='Ngenomes',ascending=False).drop(columns='Ngenomes')
order=list(pd.read_csv('orderTreeFigure.txt',header=None)[0])[::-1] # file can be found in the folder "Comparative genomics of 150 fungal secretomes/data/"
og_gc=og_gc[order]

fig,ax=plt.subplots(1,1,figsize=(2.5,10))
sns.heatmap(og_gc.T,cmap='Greys',ax=ax,cbar=False)
plt.savefig('heatmap_github.pdf')
plt.close()

