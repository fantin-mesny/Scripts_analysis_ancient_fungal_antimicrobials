import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
from Bio import SeqIO
from scipy import stats
from copy import deepcopy

def parseFasta(fa):
    fa=SeqIO.to_dict(SeqIO.parse(fa,'fasta'))
    return fa
    
def parseTMbed(filename):
    pred=pd.DataFrame(columns=['Number of TM domains'])
    with open(filename,'r') as inp:
        fileLines=inp.readlines()
    for lineN in range(len(fileLines)):
        if fileLines[lineN].startswith('>'):
            pred.loc[fileLines[lineN][1:-1],'Number of TM domains']=len([domain for domain in fileLines[lineN+2][:-1].split('.') if domain!='' and list(set(domain))[0]!='S'])
    return pred[pred['Number of TM domains']>0]

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format


org=['Verticillium','Coprinopsis','Rhizophagus']
dic={}
n=0
plddts=[]
dfs=[]
stacked_annot=[]
stacked_pred=[]
am={}
nam={}

fig,ax=plt.subplots(5,3,figsize=(8,15),height_ratios=[3, 3, 2,3,3])
for o in org:
	print(o)
    # parse secretome fasta
	fa=parseFasta('data/'+o+'/secretome.fa')
	# parse TMbed outputs and filter out transmembrane proteins
	tm=parseTMbed('data/'+o+'/tmbed.pred')
	fa_filtered=[prot for prot in fa if prot not in tm.index]
	
	# parse CAZyme annotation from dbcan
	caz=pd.read_csv('data/'+o+'/dbcan.tsv',sep='\t').set_index('Gene ID')
	caz['allAnnot']=caz['HMMER']+'/'+caz['dbCAN_sub']+'/'+caz['DIAMOND']
	caz=caz[(caz['allAnnot'].str.contains('GH')) | (caz['allAnnot'].str.contains('GT')) | (caz['allAnnot'].str.contains('PL')) | (caz['allAnnot'].str.contains('CE')) | (caz['allAnnot'].str.contains('AA'))]
	caz_rename={c:'dbcan:'+c for c in caz.columns}

	# parse annotations from EggNog
	columnsEN=['COG_category','Description','Preferred_name','GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module','KEGG_Reaction','KEGG_rclass','BRITE','KEGG_TC','CAZy','BiGG_Reaction','PFAMs']
	columnsEN_rename={c:'emapper:'+c for c in columnsEN}
	egg=pd.read_csv('data/'+o+'/emapper.annotations.tsv',sep='\t',skiprows=4,skipfooter=3,engine='python').set_index('#query')[columnsEN]
	for ind in egg.index:
		for c in columnsEN:
			if egg.loc[ind,c]!='-':
				egg.loc[ind,'non-']=egg.loc[ind,c]
	egg=egg[columnsEN+['non-']].dropna().drop(columns=['non-'])


	# assemble all annotations in dataframes 'annot' and 'df'
	annot=egg.rename(columns=columnsEN_rename).merge(caz.rename(columns=caz_rename),left_index=True,right_index=True, how='outer')
	annot['emapper:GOs']=annot['emapper:GOs'].str.replace(',',';')
	df=pd.DataFrame(index=fa_filtered,columns=['TM protein','CAZyme','Other annotation (eggnog)','Not annotated'])
	seqs={}
	TMs={}
	for prot in fa:
		if prot in tm.index:
			df.loc[prot,'TM protein']=1
			TMs[prot]='YES'
		elif prot in caz.index:
			df.loc[prot,'CAZyme']=1
		elif prot in egg.index:
			df.loc[prot,'Other annotation (eggnog)']=1
		else:
			df.loc[prot,'Not annotated']=1
			annot.loc[prot,'NA']=True
		seqs[prot]=str(fa[prot].seq)
	df=df.fillna(0)
	annot['Sequence']=annot.index.map(seqs) # add sequences
	col=list(annot.columns)
	annot['tmpred:Transmembrane']=annot.index.map(TMs)
	annot['tmpred:Transmembrane']=annot['tmpred:Transmembrane'].fillna('NO')
	col=['tmpred:Transmembrane']+col

	# export supplementary tables
	annot[col].drop(columns='NA').to_csv(o+'.suptable.tsv',sep='\t')

	# parse antimicrobial activity predictions from AMAPEC
	pred=pd.read_csv('data/'+o+'/AMAPEC_prediction.csv').set_index('Protein ID')
	pred=pred[pred.index.isin(df[(df['CAZyme']==0) & (df['TM protein']==0)].index)]
	pred_rename={c:'amapec:'+c for c in pred.columns}
	annot=annot.merge(pred.rename(columns=pred_rename),left_index=True,right_index=True,how='outer')
	plddts.append(annot[annot['dbcan:#ofTools'].isna()][['amapec:pLDDT']])
	plddts[-1]['Org']=o
	am[o]=list(pred[pred['Prediction']=='Antimicrobial'].index)
	nam[o]=list(pred[pred['Prediction']=='Non-antimicrobial'].index)

	# prepare data for stacked barplot figure
	stacked_pred.append({'Organism':o, 'Antimicrobial': len(pred[pred['Prediction']=='Antimicrobial']), 'Non-antimicrobial': len(pred[pred['Prediction']=='Non-antimicrobial'])})
	stacked_annot.append({'Organism':o,'TM protein': sum(df['TM protein']), 'CAZyme': sum(df['CAZyme']), 'Other annot': sum(df['Other annotation (eggnog)']), 'Not annot': sum(df['Not annotated'])})

	dfs.append(df)


# export figure with barcharts
stacked_annot=pd.DataFrame(stacked_annot).set_index('Organism')
stacked_pred=pd.DataFrame(stacked_pred).set_index('Organism')
stacked_pred_pc=stacked_pred.div(stacked_pred.sum(axis=1), axis=0)
fig,ax=plt.subplots(1,2,figsize=(9,4))
stacked_annot.plot.barh(ax=ax[0],stacked=True,color=['#ebc934','#91cf60', '#91bfdb','lightgrey'])
stacked_pred_pc.plot.barh(ax=ax[1],stacked=True,color=['#ff5555', '#999999'])
plt.savefig('barchart_figure.pdf')
plt.close()



plddts=pd.concat(plddts)
dfs=pd.concat(dfs)


# Analysing orthology prediction
ortho=pd.concat([pd.read_csv('data/Orthogroups.tsv',sep='\t').set_index('Orthogroup'),pd.read_csv('data/Orthogroups_UnassignedGenes.tsv',sep='\t').set_index('Orthogroup')]).fillna('')
prot2exclude=list(dfs[(dfs['CAZyme']==1) | (dfs['TM protein']==1)].index)
for og in ortho.index:
	for fung in ortho.columns:
		prots=ortho.loc[og,fung].split(', ')
		for prot in prots:
			if prot in prot2exclude:
				prots.remove(prot)
			ortho.loc[og,fung]=', '.join(prots)
ortho_am=deepcopy(ortho)
for og in ortho_am.index:
	for fung in ortho_am.columns:
		for prot in ortho_am.loc[og,fung].split(', '):
			if prot not in am[fung]:
				if ortho_am.loc[og,fung].endswith(prot):
					ortho_am.loc[og,fung]=ortho_am.loc[og,fung].replace(prot,'')
				else:
					ortho_am.loc[og,fung]=ortho_am.loc[og,fung].replace(prot+', ', '')
ortho_am[ortho_am!='']=1
ortho_am[ortho_am=='']=0
ortho_am['nbGenomes']=ortho_am.sum(axis=1)
am_vc=pd.DataFrame(ortho_am['nbGenomes'].value_counts())
am_vc.columns=['AM']

# Venn diagrams showing conserved antimicrobials and non-antimicrobials
plt.figure()
Vd = set(ortho_am[ortho_am['Verticillium']==1].index)
Ri = set(ortho_am[ortho_am['Rhizophagus']==1].index)
Cc = set(ortho_am[ortho_am['Coprinopsis']==1].index)
venn3([Vd, Ri, Cc], ('Verticillium', 'Rhizophagus', 'Coprinopsis'))
plt.savefig('venn_am.pdf')
plt.close()
ortho_nam=deepcopy(ortho)
for og in ortho_nam.index:
	for fung in ortho_nam.columns:
		for prot in ortho_nam.loc[og,fung].split(', '):
			if prot not in nam[fung]:
				if ortho_nam.loc[og,fung].endswith(prot):
					ortho_nam.loc[og,fung]=ortho_nam.loc[og,fung].replace(prot,'')
				else:
					ortho_nam.loc[og,fung]=ortho_nam.loc[og,fung].replace(prot+', ', '')
ortho_nam[ortho_nam!='']=1
ortho_nam[ortho_nam=='']=0
ortho_nam['nbGenomes']=ortho_nam.sum(axis=1)
nam_vc=pd.DataFrame(ortho_nam['nbGenomes'].value_counts())
nam_vc.columns=['NAM']
plt.figure()
Vd = set(ortho_nam[ortho_nam['Verticillium']==1].index)
Ri = set(ortho_nam[ortho_nam['Rhizophagus']==1].index)
Cc = set(ortho_nam[ortho_nam['Coprinopsis']==1].index)
venn3([Vd, Ri, Cc], ('Verticillium', 'Rhizophagus', 'Coprinopsis'))
plt.savefig('venn_nam.pdf')
plt.close()

# Barplots showing conservation of antimicrobials
bar_df=pd.concat([am_vc,nam_vc],axis=1).T.drop(columns=0.0)
bar_df['3+2']=bar_df[3.0]+bar_df[2.0]
bar_df['3+2+1']=bar_df[3.0]+bar_df[2.0]+bar_df[1.0]
bar_df=pd.DataFrame(bar_df.stack()).reset_index(drop=False).rename(columns={'level_0':'Group','level_1':'Sharedness',0:'Number of proteins'})
fig,ax=plt.subplots(1,1,figsize=(1,1.75))
sns.barplot(ax=ax,x='Group',y='Number of proteins',data=bar_df[bar_df['Sharedness']=='3+2+1'],color='lightgrey')
sns.barplot(ax=ax,x='Group',y='Number of proteins',data=bar_df[bar_df['Sharedness']=='3+2'],color='grey')
sns.barplot(ax=ax,x='Group',y='Number of proteins',data=bar_df[bar_df['Sharedness']==3.0],color='black')
plt.savefig('barplot_conservation.pdf')
plt.close()


#  Group Sharedness  Number of proteins
#0    AM        1.0                 504
#1    AM        2.0                  45
#2    AM        3.0                  13
#3    AM        3+2                  58
#4    AM      3+2+1                 562
#5   NAM        1.0                 583
#6   NAM        2.0                  25
#7   NAM        3.0                   6
#8   NAM        3+2                  31
#9   NAM      3+2+1                 614


# Fisher's test for enrichment
fet=stats.fisher_exact([[bar_df.loc[3,'Number of proteins'],bar_df.loc[8,'Number of proteins']],[bar_df.loc[0,'Number of proteins'],bar_df.loc[5,'Number of proteins']]], alternative='greater')
print(fet)
