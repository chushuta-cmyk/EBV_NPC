"""Restyle saved EBV results. No model fitting, enrichment or DEG reruns."""
from pathlib import Path
import json, shutil, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
import networkx as nx
from adjustText import adjust_text

ROOT = Path(__file__).resolve().parent
OUT = ROOT
DATA = OUT/'source_data'
BASE = Path('.')  # Inputs are resolved from packaged source_data by filename.
EXP = Path('.')  # Inputs are resolved from packaged source_data by filename.
SC = Path('.')  # Inputs are resolved from packaged source_data by filename.
for d in [OUT, DATA, OUT/'panels']: d.mkdir(parents=True, exist_ok=True)
for font in ['Arial.ttf','Arial Bold.ttf']:
    fp=Path('/System/Library/Fonts/Supplemental')/font
    if fp.exists(): font_manager.fontManager.addfont(str(fp))
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,
 'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.5,
 'axes.linewidth':.6,'lines.linewidth':.8,'pdf.fonttype':42,'ps.fonttype':42,
 'savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
COL={'Low_CP':'#377EB8','High_CP':'#E15A63'}
MM=1/25.4
manifest=[]
def source(path, name=None):
    path=Path(path); target=DATA/(name or path.name)
    if not target.exists(): raise FileNotFoundError(target)
    manifest.append({'source':str(path),'copy':str(target.relative_to(OUT))})
    return target
def table(name): return pd.read_csv(source(BASE/'output/tables'/f'{name}.csv'))
def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(width=.5,length=2,pad=2)
def label(fig,letter,x,y): fig.text(x,y,letter,weight='bold',fontsize=11,va='top')
def save(fig,num):
    fig.savefig(OUT/f'Figure_{num}.pdf',dpi=600)
    fig.savefig(OUT/f'Figure_{num}.png',dpi=600)

    plt.close(fig)
def short(s):
    return s.replace('T cells ','T ').replace('B cells ','B ').replace('NK cells ','NK ').replace('Dendritic cells ','DCs ').replace('Mast cells ','Mast ').replace('T follicular helper','T follicular helper').replace('T regulatory (Tregs)','Tregs')

def fig3():
    fig=plt.figure(figsize=(170*MM,185*MM))
    specs=[('go_upregulated','A','GO: higher in High_CP',(.56,.605,.29,.34)),
           ('go_downregulated','B','GO: lower in High_CP',(.56,.235,.29,.28))]
    for name,letter,title,rect in specs:
        d=table(name).sort_values(['p.adjust','Count'],ascending=[True,False],kind='stable').head(15)
        d.to_csv(DATA/f'Figure3_{letter}_displayed.csv',index=False)
        ax=fig.add_axes(rect);y=np.arange(len(d));ratio=d.GeneRatio.map(lambda s:eval(s,{'__builtins__':{}}))
        vmin,vmax=d['p.adjust'].min(),d['p.adjust'].max()
        sc=ax.scatter(ratio,y,s=10+3*d.Count,c=d['p.adjust'],cmap='RdYlBu',vmin=vmin,vmax=vmax,edgecolor='#333333',linewidth=.2,zorder=3)
        ax.set_yticks(y,d.Description);ax.invert_yaxis();ax.set_xlabel('Gene ratio');ax.grid(color='#e5e5e5',linewidth=.4);ax.set_axisbelow(True);ax.set_ylim(len(d)-.5,-.5)
        ax.set_xlim(ratio.min()-.018,ratio.max()+.02);clean(ax)
        label(fig,letter,.01,rect[1]+rect[3]+.048);fig.text(.065,rect[1]+rect[3]+.025,title,fontsize=9,weight='bold')
        cax=fig.add_axes((.89,rect[1]+.17,.025,.13));cb=fig.colorbar(sc,cax=cax);cb.ax.tick_params(labelsize=7);cb.ax.set_title('Adj. P',fontsize=8,pad=6)
        counts=[int(d.Count.min()),int(d.Count.max())]
        handles=[ax.scatter([],[],s=10+3*c,facecolor='#666666',label=str(c)) for c in counts]
        ax.legend(handles=handles,title='Count',loc='lower left',bbox_to_anchor=(1.07,-.02),frameon=False,borderpad=0,labelspacing=.6,handletextpad=.4,title_fontsize=8)
    d=table('kegg_downregulated');row=d.iloc[0]
    label(fig,'C',.01,.125);fig.text(.065,.102,'KEGG: lower in High_CP',fontsize=9,weight='bold')
    ax=fig.add_axes((.56,.055,.29,.04));x=int(row.GeneRatio.split('/')[0])/int(row.GeneRatio.split('/')[1]);ax.scatter([x],[0],s=50,color='#7B3294');ax.set_xlim(.18,.29);ax.set_yticks([0],[row.Description]);ax.set_ylim(-1,1);ax.set_xlabel('Gene ratio');clean(ax)
    fig.text(.88,.075,f'Count = {row.Count}\nAdj. P = {row["p.adjust"]:.4f}',fontsize=8,va='center')
    save(fig,3)

def fig4():
    nodes=table('ppi_node_centrality');edges=table('ppi_edges')
    fig=plt.figure(figsize=(170*MM,170*MM))
    for i,(direction,color) in enumerate([('Higher in High_CP','#CA5267'),('Lower in High_CP','#3B7EA8')]):
        nd=nodes[nodes.direction==direction];ed=edges[edges.direction==direction]
        G=nx.Graph();G.add_nodes_from(nd.string_id);G.add_edges_from(ed[['from','to']].itertuples(index=False,name=None))
        # Saved degree values are retained; only the network layout changes.
        ax=fig.add_axes((.06,.57 if i==0 else .065,.90,.35 if i==0 else .43));ax.axis('off')
        if i==0: pos=nx.circular_layout(G)
        else:
            order=nd.sort_values(['degree','gene']).string_id.to_list()
            outer=order[:max(0,len(order)-15)];inner=order[max(0,len(order)-15):]
            pos=nx.shell_layout(G,nlist=[inner,outer],rotate=.12)
        degrees=nd.set_index('string_id').degree;names=nd.set_index('string_id').gene.to_dict()
        nx.draw_networkx_edges(G,pos,ax=ax,width=.35,alpha=.18,edge_color='#555555')
        nx.draw_networkx_nodes(G,pos,ax=ax,node_size=[22+degrees[n]*2.6 for n in G],node_color=color,alpha=.65,linewidths=.45,edgecolors='white')
        texts=[ax.text(pos[n][0],pos[n][1]+.035,names[n],fontsize=8,ha='center',va='bottom') for n in G]
        ax.set_xlim(-1.16,1.16);ax.set_ylim(-1.16,1.16)
        adjust_text(texts,ax=ax,expand=(1.10,1.30),force_text=(.3,.5),iter_lim=500,arrowprops={'arrowstyle':'-','color':'#888888','lw':.4})
        label(fig,'AB'[i],.01,.98 if i==0 else .55);fig.text(.06,.96 if i==0 else .53,direction,fontsize=9,weight='bold')
    fig.text(.5,.018,'STRING v12.0  |  Combined score ≥ 400  |  Node size: degree',ha='center',fontsize=8)
    save(fig,4)

def fig5():
    long=table('immune_cell_proportions_long');raw=table('cibersort_results').set_index('Unnamed: 0');cells=raw.columns[:22].to_list()
    grp=table('sample_groups').set_index('sample')['group'];raw=raw[cells]
    fig=plt.figure(figsize=(170*MM,195*MM))
    ax=fig.add_axes((.08,.79,.90,.16));x=np.arange(22)
    for j,g in enumerate(['Low_CP','High_CP']):
        vals=[long[(long.group==g)&(long.cell_type==c)].proportion.to_numpy() for c in cells]
        bp=ax.boxplot(vals,positions=x+(j-.5)*.35,widths=.29,patch_artist=True,showfliers=True,
                      flierprops={'markersize':1.4,'markeredgewidth':.3},medianprops={'color':'black','linewidth':.6},
                      boxprops={'linewidth':.5},whiskerprops={'linewidth':.5},capprops={'linewidth':.5})
        for p in bp['boxes']:p.set_facecolor(COL[g])
    ax.set_xticks(x,[short(c) for c in cells],rotation=60,ha='right',fontsize=7);ax.set_ylabel('Estimated fraction');ax.set_xlim(-.7,21.7);clean(ax)
    ax.legend([plt.Rectangle((0,0),1,1,color=COL[g]) for g in COL],list(COL),ncol=2,frameon=False,loc='upper right')
    label(fig,'A',.01,.985)
    order=grp.sort_values(key=lambda s:s.map({'Low_CP':0,'High_CP':1}),kind='stable').index
    ax=fig.add_axes((.08,.455,.57,.19));bottom=np.zeros(len(order));palette=list(plt.get_cmap('tab20').colors)+['#4C4C4C','#E5C849']
    for j,c in enumerate(cells):
        vals=raw.loc[order,c].values;ax.bar(range(len(order)),vals,bottom=bottom,color=palette[j],width=.88,linewidth=0);bottom+=vals
    ids=[f'S{i+1:02d}' for i in range(len(order))]
    pd.DataFrame({'label':ids,'sample':order,'group':grp.loc[order].values}).to_csv(DATA/'Figure5_sample_key.csv',index=False)
    ax.set_xticks(range(len(order)),ids,rotation=90,fontsize=7);ax.set_xlim(-.7,len(order)-.3);ax.set_ylim(0,1);ax.set_ylabel('Estimated fraction');clean(ax)
    ax.axvline(4.5,color='black',lw=.6);ax.text(2,1.03,'Low_CP',ha='center',fontsize=8);ax.text(13,1.03,'High_CP',ha='center',fontsize=8)
    lax=fig.add_axes((.67,.435,.32,.235));lax.axis('off');lax.legend([plt.Rectangle((0,0),1,1,color=c) for c in palette],[short(c) for c in cells],ncol=1,frameon=False,loc='center left',fontsize=6.7,labelspacing=.04,handlelength=.8,handletextpad=.4,borderpad=0)
    label(fig,'B',.01,.68)
    z=raw.T.sub(raw.T.mean(axis=1),axis=0).div(raw.T.std(axis=1,ddof=1),axis=0).fillna(0)
    zr=linkage(z,method='complete');zc=linkage(z.T,method='complete');ir=leaves_list(zr);ic=leaves_list(zc)
    ax=fig.add_axes((.32,.04,.60,.32));im=ax.imshow(z.iloc[ir,ic],aspect='auto',cmap='RdBu_r',vmin=z.min().min(),vmax=z.max().max())
    ax.set_yticks(range(22),[short(s) for s in z.index[ir]],fontsize=7.5);ax.set_xticks([]);ax.tick_params(length=0);ax.spines[:].set_visible(False)
    da=fig.add_axes((.32,.365,.60,.037));dendrogram(zc,ax=da,no_labels=True,color_threshold=0,above_threshold_color='#555555');da.axis('off')
    ga=fig.add_axes((.32,.355,.60,.008));ga.imshow(np.array([matplotlib.colors.to_rgb(COL[grp[s]]) for s in z.columns[ic]])[None,:,:],aspect='auto');ga.axis('off')
    ca=fig.add_axes((.945,.08,.02,.20));cb=fig.colorbar(im,cax=ca,ticks=[-1,0,1,2,3,4]);cb.ax.tick_params(labelsize=7);cb.ax.set_title('Z',fontsize=8,pad=5)
    label(fig,'C',.01,.405);save(fig,5)

def fig6():
    ml=BASE/'python_ml/output';perf=pd.read_csv(source(ml/'ML_Model_Performance.csv'));corr=pd.read_csv(source(ml/'Gene_Immune_Spearman_Correlation.csv'),index_col=0);q=pd.read_csv(source(ml/'Gene_Immune_Adjusted_P_Values.csv'),index_col=0)
    fig=plt.figure(figsize=(170*MM,155*MM));ax=fig.add_axes((.24,.72,.69,.22));d=perf.iloc[::-1]
    ax.barh(d.Model,d.LOO_Balanced_Accuracy,color='#689AC2',edgecolor='black',linewidth=.4);ax.axvline(.5,ls='--',color='#888888');ax.set_xlim(0,1);ax.set_xlabel('Leave-one-out balanced accuracy');clean(ax)
    for i,v in enumerate(d.LOO_Balanced_Accuracy):ax.text(v+.01,i,f'{v:.3f}',va='center',fontsize=8)
    label(fig,'A',.01,.99)
    ax=fig.add_axes((.31,.15,.58,.45));im=ax.imshow(corr,aspect='auto',cmap='RdBu_r',vmin=-.8,vmax=.8)
    ax.set_xticks(range(corr.shape[1]),corr.columns,rotation=90,fontsize=8);ax.set_yticks(range(corr.shape[0]),[short(s) for s in corr.index],fontsize=8);ax.tick_params(length=0,pad=3)
    ax.spines[:].set_visible(False)
    for i in range(q.shape[0]):
        for j in range(q.shape[1]):
            t='**' if q.iloc[i,j]<.01 else '*' if q.iloc[i,j]<.05 else ''
            if t:ax.text(j,i,t,ha='center',va='center',fontsize=8,color='white' if abs(corr.iloc[i,j])>.55 else 'black')
    ca=fig.add_axes((.92,.19,.018,.36));cb=fig.colorbar(im,cax=ca);cb.ax.tick_params(labelsize=7,pad=2);cb.ax.set_title('ρ',fontsize=9)
    label(fig,'B',.01,.65);fig.text(.5,.02,'Spearman correlation; * adjusted P < 0.05, ** adjusted P < 0.01',ha='center',fontsize=8)
    save(fig,6)

def km_values(t,e):
    times=np.unique(t);s=1.;green=0.;out=[[0,1,1,1]]
    for ti in times:
        n=np.sum(t>=ti);d=np.sum((t==ti)&(e==1))
        if d:
            s*=1-d/n
            if n>d:green+=d/(n*(n-d))
        # Greenwood log confidence interval, matching survival::survfit default.
        se=np.sqrt(green);lo=s*np.exp(-1.95996398454*se);hi=min(1,s*np.exp(1.95996398454*se))
        out.append([ti,s,lo,hi])
    return np.array(out)
def fig7():
    genes=['IL2RG','CD247','THBD','ITGAX','CD3E','C5AR1','NLRP3','IL1A','CSF3R','LCK','MS4A1','HLA-DPB1','CXCR2','TREM1','S100A8','S100A9','IL1B']
    dat=pd.read_csv(DATA/'Figure7_survival_input.csv')
    fig,axs=plt.subplots(6,3,figsize=(170*MM,180*MM));fig.subplots_adjust(left=.10,right=.985,bottom=.065,top=.94,hspace=.60,wspace=.25)
    records=[]
    for i,g in enumerate(genes):
        ax=axs.flat[i];high=dat[g]>dat[g].median()
        for mask,c,group in [(high,'#377EB8','High'),(~high,'#E15A63','Low')]:
            a=dat[mask];k=km_values(a.time.values,a.event_code.values);ax.step(k[:,0],k[:,1],where='post',color=c,lw=.9);ax.fill_between(k[:,0],k[:,2],k[:,3],step='post',color=c,alpha=.10,linewidth=0)
            records.extend([{'panel':chr(65+i),'gene':g,'group':group,'time':r[0],'survival':r[1],'lower95':r[2],'upper95':r[3]} for r in k])
        ax.set_title(f'{chr(65+i)}  {g}',loc='left',fontsize=8.5,pad=3,weight='bold');ax.set_ylim(0,1.03);ax.set_xlim(0,55);ax.set_xticks([0,25,50]);ax.set_yticks([0,.5,1]);clean(ax)
        if i%3:ax.set_yticklabels([])
    axs.flat[-1].axis('off');fig.text(.025,.51,'Progression-free survival probability',rotation=90,va='center',fontsize=9);fig.text(.54,.015,'Time to event',ha='center',fontsize=9)
    fig.legend([Line2D([0],[0],color=c) for c in ['#377EB8','#E15A63']],['High expression','Low expression'],ncol=2,loc='upper center',frameon=False,bbox_to_anchor=(.54,.995),fontsize=8)
    pd.DataFrame(records).to_csv(DATA/'Figure7_KM_coordinates.csv',index=False);save(fig,7)

def fig9():
    # Preserve complete, unmodified raster panels recovered from the earlier source PDF.
    original0=source(ROOT/'work/fig9_original_0.png');original1=source(ROOT/'work/fig9_original_1.tif')
    prot=pd.read_csv(DATA/'Figure9_proteomics.csv')
    core=['FN1','S100P','EPHA2','ITGA3','THBS1','S100A10','STAT3','STMN1','STAT6','STAT1','PFN2'];p=prot[prot.Gene.str.upper().isin(core)].drop_duplicates('Gene').sort_values('log2FC')
    cyt=pd.read_csv(source(EXP/'EBV_recoding/cytokine_replicate.csv'),index_col=0);names=list(dict.fromkeys(c.rsplit('_R',1)[0] for c in cyt.columns));means=pd.DataFrame({c:cyt[[v for v in cyt if v.rsplit('_R',1)[0]==c]].mean(axis=1) for c in names}).T
    z=means.sub(means.mean(axis=1),axis=0).div(means.std(axis=1,ddof=1),axis=0);ir=leaves_list(linkage(z,method='complete'));ic=leaves_list(linkage(z.T,method='complete'))
    fig=plt.figure(figsize=(170*MM,185*MM))
    ax=fig.add_axes((.02,.755,.54,.235));ax.imshow(plt.imread(original0));ax.axis('off');label(fig,'A',.01,.99)
    ax=fig.add_axes((.63,.755,.25,.235));ax.imshow(plt.imread(original1));ax.axis('off');label(fig,'B',.56,.99)
    ax=fig.add_axes((.09,.44,.38,.275));ax.scatter(prot.log2FC,prot.negLog10Q,color='#d8d8d8',s=3,alpha=.5,rasterized=True)
    for _,r in p.iterrows():
        c='#B2182B' if r.log2FC>0 else '#2166AC';ax.scatter(r.log2FC,r.negLog10Q,s=20,c=c,edgecolors='black',linewidths=.4)
    # Draw labels in separated columns, with leader lines to exact data coordinates.
    for positive,xx,ha in [(False,-3.7,'left'),(True,3.3,'left')]:
        dd=p[p.log2FC.gt(0)==positive].sort_values('negLog10Q')
        ys=np.linspace(16,48,len(dd))
        for (_,r),yy in zip(dd.iterrows(),ys):ax.annotate(r.Gene,(r.log2FC,r.negLog10Q),xytext=(xx,yy),fontsize=7.5,ha=ha,va='center',arrowprops={'arrowstyle':'-','lw':.4,'color':'#666666'})
    ax.axhline(-np.log10(.05),ls='--',c='gray',lw=.5)
    for x in [-np.log2(1.5),np.log2(1.5)]:ax.axvline(x,ls='--',c='gray',lw=.5)
    ax.set_xlim(-5.3,7.4);ax.set_ylim(0,prot.negLog10Q.max()+2);ax.set_xlabel('log2 fold change');ax.set_ylabel('−log10 Q value');clean(ax);label(fig,'C',.01,.75)
    ax=fig.add_axes((.62,.44,.34,.275));ax.barh(p.Gene,p.log2FC,color=['#2166AC' if v<0 else '#B2182B' for v in p.log2FC],height=.6)
    for i,v in enumerate(p.log2FC):ax.text(v+(.08 if v>=0 else -.08),i,f'{v:+.2f}',ha='left' if v>=0 else 'right',va='center',fontsize=7.5)
    ax.axvline(0,c='black',lw=.6);ax.set_xlim(-2.3,3.6);ax.set_xlabel('log2 fold change');clean(ax);label(fig,'D',.53,.75)
    ax=fig.add_axes((.08,.08,.19,.285));im=ax.imshow(z.iloc[ir,ic],cmap='RdBu_r',vmin=-1.2,vmax=1.2,aspect='auto');ax.set_xticks(range(3),z.columns[ic],rotation=45,ha='right');ax.set_yticks(range(18),z.index[ir],fontsize=7.5);ax.yaxis.tick_right();ax.tick_params(length=0);ax.spines[:].set_visible(False)
    ca=fig.add_axes((.37,.13,.017,.18));cb=fig.colorbar(im,cax=ca,ticks=[-1,0,1]);cb.ax.set_title('Z',fontsize=8);label(fig,'E',.01,.40)
    top=(means.TPA/means.UT).sort_values(ascending=False).head(8).index
    for j,c in enumerate(top):
        rr,cc=divmod(j,4);ax=fig.add_axes((.51+cc*.12,.245-rr*.16,.09,.115));cols=[v for v in cyt if v.rsplit('_R',1)[0]==c]
        v=cyt.loc[['DMSO','TPA','UT'],cols];ax.bar(range(3),v.mean(axis=1),yerr=v.std(axis=1,ddof=1),color=['#aaaaaa','#E79232','#6AAED6'],width=.7,error_kw={'elinewidth':.5,'capsize':1.5})
        ax.set_xticks(range(3),['DMSO','TPA','UT'] if rr else ['','',''],rotation=65,ha='right',fontsize=7);ax.set_title(c,fontsize=8,pad=4);ax.tick_params(labelsize=7,pad=1);ax.yaxis.set_major_locator(plt.MaxNLocator(3));clean(ax)
    label(fig,'F',.46,.40);fig.text(.49,.38,'Cytokines (pg/mL); mean ± SD',fontsize=8);fig.text(.52,.015,'Proteomic contrast: NPC43-Re / NPC43-Ctl',ha='center',fontsize=8)
    save(fig,9)

if __name__=='__main__':
    import sys
    for n in sys.argv[1:] or ['3','4','5','6','7','9']:
        print('Drawing Figure',n,flush=True);globals()['fig'+n]()
    (DATA/'replot_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
