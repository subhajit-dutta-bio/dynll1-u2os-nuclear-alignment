"""Shared plotting and panel-source conventions for the revised manuscript."""
from pathlib import Path
import json,textwrap,re
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import AutoLocator,MaxNLocator
from matplotlib.lines import Line2D
import tifffile
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=ROOT/'phase14';R=P/'results';F=OUT/'figures';D=OUT/'panel_data';A=OUT/'results'
COLORS={'anchor':'#236A98','motor_context':'#C56B36','proximal_12q':'#9467A8','distal_12q':'#41957A','matched_other_chromosome':'#84939D'}
NAMES={'anchor':'DYNLL1','motor_context':'Motor context','proximal_12q':'Within 2 Mb','distal_12q':'Elsewhere 12q','matched_other_chromosome':'Matched other chr.'}
DATES=['20220914_Run1','20221009_Run2','20221017_Run3','20221024_Run4','20221109_Run5'];DATECOL=['#327DA4','#CC783E','#51946E','#9475AA','#B2A14F'];REG=[]
SCALES=json.loads((ROOT/'phase5/jump_raw_policy.json').read_text())['source_control_scales']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':9,'axes.labelsize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
def panel(ax,letter,title):ax.set_title(letter+'  '+textwrap.fill(title,43),loc='left',fontweight='normal',pad=7)
def zero(ax):ax.axhline(0,color='#555555',lw=.7,zorder=0)
def diagonal(ax,x,y):
    low=float(np.nanmin(np.r_[x,y]));high=float(np.nanmax(np.r_[x,y]));ax.plot([low,high],[low,high],color='#777777',lw=.8,zorder=0)
def register(key,letter,title,data,unit,interpretation,source):
    D.mkdir(exist_ok=True);path=D/f'{key}_{letter}.csv'
    if not isinstance(data,pd.DataFrame):data=pd.DataFrame(data)
    assert len(data)>0,(key,letter);data.to_csv(path,index=False);REG.append({'figure':key,'panel':letter,'title':title,'source_data':path.relative_to(ROOT).as_posix(),'input_source':source,'sample_unit':unit,'interpretation':interpretation})
def grid():
    fig,axes=plt.subplots(4,2,figsize=(7.4,9.5));fig.subplots_adjust(left=.14,right=.97,top=.965,bottom=.06,wspace=.52,hspace=.78);return fig,list(axes.flat)

def group_legend(fig,title='Gene groups'):
    handles=[Line2D([],[],marker='o',linestyle='none',color=color,markersize=4,label=NAMES[group]) for group,color in COLORS.items()]
    fig.legend(handles=handles,ncol=5,frameon=False,loc='upper center',bbox_to_anchor=(.5,-.012),fontsize=8,title=title,title_fontsize=8,handletextpad=.3,columnspacing=.9)
def finish(fig,key):
    F.mkdir(exist_ok=True)
    from matplotlib.text import Text
    from matplotlib.image import AxesImage
    for ax in fig.axes:
        if isinstance(ax.xaxis.get_major_locator(),AutoLocator):ax.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3))
    for artist in fig.findobj(match=AxesImage):
        cmap=artist.get_cmap().copy();cmap.set_bad('#E0E0E0');artist.set_cmap(cmap)
    for obj in fig.findobj(match=Text):
        if obj.get_text() and obj.get_fontsize()<7.5:obj.set_fontsize(7.5)
        label=obj.get_text();label=re.sub(r'(?i)proximal\s*12q','Within 2 Mb',label);label=re.sub(r'(?i)distal\s*12q','Elsewhere 12q',label);label=re.sub(r'\blog2\b','log₂',label);label=re.sub(r'\blog10\b','log₁₀',label)
        for word in ['chromosome','Day','Replicate','Rep','site','proximal','distal','All','across','within','Manual','Watershed','SAM']:label=re.sub(r'\b('+word+r')(?=\d)',r'\1 ',label,flags=re.I)
        label=re.sub(r'(\d)(µg|µm)\b',r'\1 \2',label)
        label=re.sub(r'\brho(?=\d)','ρ = ',label);label=re.sub(r'\bq(?=\d)','q = ',label);label=label.replace('Intended chromosome 12 cut coordinate, Mb','Brunello gene-position proxy, Mb').replace('Intended cut coordinate on chromosome 12, Mb','Brunello median cut coordinate, Mb');obj.set_text(label)
    if key=='main_figure2':
        for ax in fig.axes:
            if len(ax.get_yticklabels())>=4:
                for tick in ax.get_yticklabels():tick.set_fontsize(9.5)
    for ext in ['png','svg','pdf']:fig.savefig(F/f'{key}.{ext}',dpi=320,bbox_inches='tight',pad_inches=.05)
    plt.close(fig)
def image_panel(ax,row,label,control=False):
    field=row.field;raw=tifffile.imread(P/'data/images'/f'{field}.tif');source=f'{row.plate}__A02__s0';reference=tifffile.imread(P/'data/images'/f'{source}.tif');lo,hi=np.quantile(reference,[.01,.995]);crop=raw[284:796,384:896];ax.imshow(crop,cmap='gray',vmin=lo,vmax=hi,interpolation='none');g=pd.read_parquet(R/'geometry'/f'{field}.parquet');scale=SCALES[row.batch]['area'];valid=g.area.between(.1*scale,10*scale)&np.isfinite(g[['x','y','theta','ar']]).all(axis=1)&g.minor.gt(0)&g.x0.gt(0)&g.y0.gt(0)&g.x1.lt(1280)&g.y1.lt(1080);g=g[valid&g.ar.ge(1.2)&g.x.between(384,896)&g.y.between(284,796)].copy();theta=np.pi/2-g.theta.to_numpy();dx=.5*g.major.to_numpy()*np.cos(theta);dy=.5*g.major.to_numpy()*np.sin(theta);xy=g[['x','y']].to_numpy()-[384,284]
    if len(g):ax.add_collection(LineCollection(np.stack([xy-np.column_stack([dx,dy]),xy+np.column_stack([dx,dy])],axis=1),colors='#62C5DB',linewidths=.45,alpha=.9))
    ax.plot([20,120],[490,490],color='white',lw=2);ax.text(20,477,'65 µm',color='white',fontsize=8);ax.set(xlim=(0,512),ylim=(512,0));ax.set_xticks([]);ax.set_yticks([])
    source=g[['nucleus_id','x','y','major','theta','ar']].copy()
    if not len(source):source=pd.DataFrame([{'nucleus_id':np.nan}])
    source['field']=field;source['gene']=label;source['crop_x0']=384;source['crop_y0']=284;source['crop_width']=512;source['display_low']=lo;source['display_high']=hi;source['pixel_size_um']=.65;source['raw_path']='phase14/data/images/'+field+'.tif';return source
def paired_dates(ax,data,x='model',y='effect',order=None):
    order=order or data[x].drop_duplicates().tolist()
    for i,date in enumerate(DATES):
        q=data[data.batch.eq(date)].set_index(x).reindex(order);ax.plot(range(len(order)),q[y],'-o',color=DATECOL[i],ms=3,lw=1,label=f'R{i+1}')
    ax.set_xticks(range(len(order)),order);zero(ax)
class Data:
    def __init__(self):
        self.genes=pd.read_csv(R/'gene_context_and_covariates.csv');assert len(self.genes)==67
        self.alignment=pd.read_csv(R/'alignment_contrasts.csv');self.primary=self.alignment[self.alignment.model.eq('adaptive_conditional')];assert len(self.primary)==335
        self.fields=pd.read_csv(R/'field_descriptors.csv').fillna({'gene':'','panel':''});assert len(self.fields)==7875
        self.descriptors=pd.read_csv(R/'descriptor_contrasts.csv');self.design=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''});self.groups=pd.read_csv(ROOT/'phase13/results/proposed_comparator_gene_sets.csv');self.dates=DATES;self.summary=json.loads((R/'comparison_summary.json').read_text());self.support=pd.read_csv(R/'density_area_support.csv');self.pred=pd.read_csv(R/'prediction/summary.csv');self.oof=pd.read_csv(R/'prediction/OOF_predictions.csv')
    def fixed(self,gene,control=False):
        q=self.design[self.design.gene.eq(gene)&self.design.batch.eq(DATES[0])&self.design.site.eq(0)];assert len(q)==1;row=q.iloc[0]
        if control:row=self.design[self.design.plate.eq(row.plate)&self.design.well.eq('A02')&self.design.site.eq(0)].iloc[0]
        return row
    def gene_scatter(self,ax,x,y):
        for group,color in COLORS.items():
            q=self.genes[self.genes.panel.eq(group)];ax.scatter(q[x],q[y],s=42 if group=='anchor' else 19,color=color,alpha=.95 if group=='anchor' else .65,zorder=4 if group=='anchor' else 2,label=NAMES[group])
        q=self.genes[self.genes.gene.eq('DYNLL1')].iloc[0];low=self.genes[y].min();high=self.genes[y].max();top=q[y]>low+.85*(high-low);ax.annotate('DYNLL1',(q[x],q[y]),xytext=(5,-6 if top else 5),va='top' if top else 'bottom',textcoords='offset points',fontsize=7.5,color=COLORS['anchor'])
