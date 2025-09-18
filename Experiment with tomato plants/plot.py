import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df=pd.read_csv('data/WPFMAllReps.csv')
df['Treatment']=df['Condition']+'_'+df['Genotype']
order=['Sterile_Mock','Sterile_JR2','Sterile_dAve1','Sterile_424Y','Recol_Mock','Recol_JR2','Recol_dAve1','Recol_424Y']

fig,ax=plt.subplots(1,1,figsize=(10,5))
sns.boxplot(ax=ax,data=df,x='Treatment',y='Freshweight',order=order,fliersize=0)
sns.swarmplot(ax=ax,data=df,x='Treatment',y='Freshweight',order=order,size=3.8,color='black')
ax.set_ylim(0,1400)
plt.savefig('plot4.pdf')
plt.close()
