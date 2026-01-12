import pandas as pd
import os
from scipy.spatial.distance import braycurtis, pdist, squareform
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def pcoa_from_braycurtis(df):
    """
    df: pandas DataFrame, samples × features
    Returns:
        coords: DataFrame of PCoA coordinates
        var_exp: Series of variance explained per axis
    """
    # --- Step 1: Compute Bray–Curtis distances ---
    D = squareform(pdist(df.values, metric="braycurtis"))
    # --- Step 2: Gower’s double-centering to convert distances -> inner product matrix ---
    # B = -0.5 * J * D^2 * J
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    # --- Step 3: Eigendecomposition ---
    eigvals, eigvecs = np.linalg.eigh(B)
    # Sort eigenvalues & eigenvectors descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    # Keep only positive eigenvalues
    pos = eigvals > 0
    eigvals = eigvals[pos]
    eigvecs = eigvecs[:, pos]
    # --- Step 4: Compute coordinates ---
    coords = eigvecs * np.sqrt(eigvals)
    # --- Step 5: Variance explained ---
    var_exp = eigvals / eigvals.sum()
    # Return top 2 axes
    coords_df = pd.DataFrame(coords[:, :2], index=df.index, columns=["PC1", "PC2"])
    var_df = pd.Series(var_exp[:2], index=["PC1_var", "PC2_var"]) * 100
    return coords_df, var_df, D



# combine number of reads per taxon identified by EMU in a single dataframe
ranks=['superkingdom','kingdom','phylum','class','order','family','genus','species']
rankOI='genus'
dfs=[]
dfs_notAggregated=[]
for i in os.listdir('./data'):
    if i.endswith('_rel-abundance.tsv'):
        df=pd.read_csv('./data/'+i,sep='\t')
        for rank in ranks:
            df[rank]=df['lineage'].str.split(';').str[:ranks.index(rank)]
            df[rank]=df[rank].astype(str)
            df[rank]=df[rank].str[2:-2].str.replace("', '",'_')
        genus_level=pd.DataFrame(df[[rankOI,'estimated counts']].groupby(rankOI).sum())
        genus_level=genus_level.rename(columns={'estimated counts':i.split('.')[1]})
        dfs_notAggregated.append(df.set_index('tax_id')[['estimated counts']].rename(columns={'estimated counts':i.split('.')[1]}))
        dfs.append(genus_level)
df=pd.concat(dfs,axis=1).fillna(0)
df=df[sorted(list(df.columns))]
df_notAggregated=pd.concat(dfs_notAggregated,axis=1).fillna(0)
df_notAggregated=df_notAggregated[sorted(list(df.columns))]
df_notAggregated=df_notAggregated.T

# calculate PCoA
pcoa_res = pcoa_from_braycurtis(df_notAggregated)
pcoa_coord=pcoa_res[0]
meta=pd.read_csv('data/metadata.csv').set_index('Sample')
pcoa_coord=pcoa_coord.merge(meta,left_index=True,right_index=True)
pcoa_coord['Barcode']=pcoa_coord.index.str[-2:]
vars=pcoa_res[1].to_dict()
distMatrix=pd.DataFrame(pcoa_res[2])
distMatrix.index=df_notAggregated.index
distMatrix.columns=df_notAggregated.index

#distMatrix.to_csv('distMatrix_brayCurtis.csv') # Used to calculate PERMANOVA in R (with vegan adonis2). Results:

### Permutation test for adonis under reduced model
### Terms added sequentially (first to last)
### Permutation: free
### Number of permutations: 999

### adonis2(formula = distmat ~ Condition, data = design, permutations = 999)
###           Df  SumOfSqs      R2      F Pr(>F)  
### Condition  1 0.0099469 0.37724 3.6346  0.046 *
### Residual   6 0.0164203 0.62276                
### Total      7 0.0263672 1.00000                
### ---
### Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# parse results of DESeq
deseq=pd.read_csv('data/DESeq_results.csv').set_index('Unnamed: 0')
deseq['Genus']=deseq.index.str.split('_').str[-2]+' '+deseq.index.str.split('_').str[-1]
deseq_top10=deseq.sort_values(by='log2FoldChange',ascending=True)[:10].sort_values(by='log2FoldChange',ascending=False)

#plot figures
fig,ax=plt.subplots(1,2,figsize=(10,4))
sns.scatterplot(x='PC1',y='PC2',hue='Condition',data=pcoa_coord,ax=ax[0],palette={'JR2':'black','d424y':'grey'})
ax[0].set_xlabel('PC1 ('+str(round(vars['PC1_var'],2))+'%)')
ax[0].set_ylabel('PC2 ('+str(round(vars['PC2_var'],2))+'%)')
ax[1].barh(deseq_top10['Genus'],deseq_top10['log2FoldChange'], xerr=deseq_top10['lfcSE'], capsize=1, color='grey')
ax[1].set_xlim(-32,0)
plt.savefig('Figure.pdf')
plt.close()