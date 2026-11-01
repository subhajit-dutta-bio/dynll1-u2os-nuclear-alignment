"""Inspect only the prespecified first-date raw examples, without gene ranking."""
from figure_common import *
def main():
    d=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''});fig,axs=plt.subplots(2,3,figsize=(9,6));rows=[]
    for i,gene in enumerate(['DYNLL1','DYNLL1','DYNC1H1','DYNC1I2','PAFAH1B1']):
        row=d[d.gene.eq(gene)&d.batch.eq(DATES[0])&d.site.eq(0)].iloc[0]
        if i==0:row=d[d.plate.eq(row.plate)&d.well.eq('A02')&d.site.eq(0)].iloc[0]
        label='Source nontargeting' if i==0 else gene+' targeted';q=image_panel(axs.flat[i],row,label);axs.flat[i].set_title(label,fontsize=10);rows.append(q)
    axs.flat[5].axis('off');axs.flat[5].text(0,.9,'Fixed image examples\n\nFirst date and site0\nCentered512-pixel native crops\nSource-control intensity limits\nCyan: eligible nuclear axes\nWhite bar:65µm\n\nNo spatial effect ranking\nis performed for this proof.',transform=axs.flat[5].transAxes,va='top',fontsize=10);fig.tight_layout();fig.savefig(P/'qa/fixed_image_examples.png',dpi=220);plt.close(fig);pd.concat(rows,ignore_index=True).to_csv(P/'qa/fixed_image_examples_source.csv',index=False)
if __name__=='__main__':main()
