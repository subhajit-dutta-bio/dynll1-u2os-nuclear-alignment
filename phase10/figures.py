from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results';F=P/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.dpi':220})
BLUE='#286894';ORANGE='#BA6636';GRAY='#78838A'
def save(fig,name):
    fig.savefig(F/(name+'.png'),bbox_inches='tight',facecolor='white');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def main():
    a=pd.read_csv(ROOT/'phase6/results/adaptive_conditional_contrasts.csv');order=['DYNLL1','TP53BP1','H2AX','RNF8','RNF168','MRE11','RAD50','CYP7B1','MAT1A','SLC35A2','SMARCE1'];fig,ax=plt.subplots(figsize=(7,4.5))
    for y,g in enumerate(order):
        q=a[a.gene.eq(g)&a.evaluable];v=q.effect.to_numpy();color=BLUE if g=='DYNLL1' else GRAY
        ax.scatter(v,y+np.linspace(-.14,.14,len(v)),s=20,color=color,alpha=.65);ax.scatter(np.mean(v),y,marker='D',s=40,color=color,zorder=3)
    ax.axvline(0,color='black',lw=.8);ax.set_yticks(range(len(order)),order);ax.invert_yaxis();ax.set_xlabel('Adjusted nuclear alignment contrast');ax.set_title('Date-level contrasts in U2OS',loc='left',fontweight='bold');ax.text(.99,.02,'Dots are dates   Diamonds are date means\nOne pooled reagent per gene',ha='right',va='bottom',transform=ax.transAxes,fontsize=8);ax.grid(axis='x',alpha=.18);fig.tight_layout();save(fig,'figure1_spatial_effects');a.to_csv(F/'figure1_source.csv',index=False)
    audit=pd.read_csv(R/'hTERT_DYNLL1_DYNLL2_reagent_and_interaction_audit.csv');dyn=audit[audit.gene.eq('DYNLL1')].sort_values('cut_position');allgenes=pd.read_csv(R/'hTERT_all_gene_interactions.csv.gz');summary=json.loads((R/'hTERT_candidate_audit.json').read_text());fig,axes=plt.subplots(1,3,figsize=(8.3,4.4),gridspec_kw={'width_ratios':[1.25,1,1.3]});ax=axes[0]
    colors=[BLUE,ORANGE,'#517C4E','#896D9C']
    for i,(_,q) in enumerate(dyn.iterrows()):
        yy=[np.mean([q.LFC_r1_WT,q.LFC_r2_WT]),np.mean([q.LFC_r1_KO,q.LFC_r2_KO])];ax.plot([0,1],yy,color=colors[i],lw=1.3,marker='o',label=str(int(q.cut_position)))
        for j,context in enumerate(['WT','KO']):ax.scatter([j-.035,j+.035],[q[f'LFC_r1_{context}'],q[f'LFC_r2_{context}']],s=10,color=colors[i],alpha=.55)
    ax.axhline(0,color='black',lw=.6);ax.set_xticks([0,1],['Parental','hTERT\nengineered']);ax.set_ylabel('DYNLL1 guide log2 abundance ratio');ax.set_title('A  Four guides',loc='left',fontweight='bold');ax.legend(title='Cut coordinate',fontsize=8,title_fontsize=8,frameon=False,loc='upper left',bbox_to_anchor=(-.05,-.24),ncol=2,columnspacing=.6,handlelength=1,borderaxespad=0);ax=axes[1]
    for i,(_,q) in enumerate(dyn.iterrows()):
        v=[q.delta_r1,q.delta_r2];ax.scatter(v,[i-.07,i+.07],color=colors[i],s=15);ax.scatter(np.mean(v),i,marker='D',color=colors[i],s=24)
    ax.axvline(0,color='black',lw=.6);ax.set_yticks(range(4),[str(int(x)) for x in dyn.cut_position],fontsize=9);ax.set_xlabel('Mutant minus parental\nlog2 ratio');ax.set_title('B  Guide interactions',loc='left',fontweight='bold');ax.invert_yaxis();ax=axes[2]
    values=allgenes.mean_delta.dropna();ax.hist(values,bins=np.linspace(-2,2,65),color='#CBD4DA',edgecolor='white');ax.axvline(summary['DYNLL1_paired_median_interaction_mean'],color=BLUE,lw=2,label='DYNLL1');d2=float(allgenes.loc[allgenes.gene.eq('DYNLL2'),'mean_delta'].iloc[0]);ax.axvline(d2,color=ORANGE,lw=1.4,label='DYNLL2');ax.set_xlabel('Gene median interaction');ax.set_ylabel('Genes');ax.set_title('C  Genome context',loc='left',fontweight='bold');ax.legend(frameon=False,fontsize=8);fig.tight_layout(w_pad=1.2);save(fig,'figure2_hTERT_context');audit.to_csv(F/'figure2_source_guides.csv',index=False)
    net=pd.read_csv(R/'APMS_paralog_neighbors.csv');genes=['ATMIN','TP53BP1','FAM83D','HMMR','TLK1','TLK2','ZMYM2','ZMYM3','ZMYM4','ZMYM6','RCOR1','RCOR2'];matrix=net.set_index('neighbor').reindex(genes)[['DYNLL1_as_bait','DYNLL2_as_bait']].astype(float);fig,(ax,bx)=plt.subplots(1,2,figsize=(7,4.6),gridspec_kw={'width_ratios':[1,1.7]});ax.imshow(matrix,cmap=matplotlib.colors.ListedColormap(['#F0F0F0',BLUE]),vmin=0,vmax=1,aspect='auto');ax.set_yticks(range(len(genes)),genes);ax.set_xticks([0,1],['DYNLL1 bait','DYNLL2 bait']);ax.tick_params(length=0);ax.set_title('A  Selected associations',loc='left',fontweight='bold');counts=[int(net['class'].eq(k).sum()) for k in ['DYNLL1_only','shared','DYNLL2_only']];bx.bar(['DYNLL1\nonly','Shared','DYNLL2\nonly'],counts,color=[BLUE,GRAY,ORANGE],width=.65)
    for i,v in enumerate(counts):bx.text(i,v+1,str(v),ha='center')
    bx.set_ylim(0,94);bx.set_ylabel('Unique neighboring proteins');bx.set_title('B  Complete neighbor overlap',loc='left',fontweight='bold');bx.text(.5,.93,'Both paralogs were affinity baits',ha='center',transform=bx.transAxes,fontsize=8);fig.tight_layout(w_pad=2);save(fig,'figure3_APMS_context');matrix.to_csv(F/'figure3_source_selected_associations.csv')
    effects=pd.read_csv(R/'DYNLL1_DYNLL2_all_contrasts.csv');names=list(dict.fromkeys(effects.contrast));mat=effects.pivot(index='contrast',columns='gene',values='mean_LFC').reindex(names);fig,(ax,bx)=plt.subplots(1,2,figsize=(10,6),gridspec_kw={'width_ratios':[1,1.2]});im=ax.imshow(mat,cmap='RdBu_r',vmin=-4,vmax=4,aspect='auto');ax.set_yticks(range(len(mat)),[n.replace('GSE','').replace('_',' ') for n in names],fontsize=10);ax.set_xticks([0,1],mat.columns,fontsize=11);ax.set_title('A  All screen contrasts',loc='left',fontweight='bold',fontsize=14)
    for i in range(len(mat)):
        for j in range(2):ax.text(j,i,f'{mat.iloc[i,j]:.2f}',ha='center',va='center',fontsize=11,color='white' if abs(mat.iloc[i,j])>2 else 'black')
    h=pd.read_csv(R/'GSE290538_HIF_high_vs_low_DYNLL1_guide_audit.csv');vals=h[['LFC_r1','LFC_r2','LFC_r3']].to_numpy();bx.imshow(vals,cmap='RdBu_r',vmin=-10,vmax=10,aspect='auto');bx.set_yticks(range(len(h)),[s.split('_',1)[1][:8]+(' *' if not e else '') for s,e in zip(h.guide,h.eligible)],fontsize=11);bx.set_xticks([0,1,2],['Rep 1','Rep 2','Rep 3'],fontsize=11);bx.set_title('B  HIF guide heterogeneity',loc='left',fontweight='bold',fontsize=14)
    for i in range(len(h)):
        for j in range(3):bx.text(j,i,f'{vals[i,j]:.1f}',ha='center',va='center',fontsize=11,color='white' if abs(vals[i,j])>5 else 'black')
    bx.set_xlabel('* Fails denominator coverage rule\nShown but excluded from primary median',fontsize=10);fig.tight_layout(w_pad=2);save(fig,'supplement_figure1_complete_context')
    print('Four PNG and SVG figure pairs created')
if __name__=='__main__':main()
