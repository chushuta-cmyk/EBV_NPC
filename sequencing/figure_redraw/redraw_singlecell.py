"""Replot Figure 8 from saved coordinates, metadata and summary tables."""
from redraw_figures import *
import re
allobs=pd.read_csv(DATA/'Figure8_all_cell_coordinates_metadata.csv')
ut=pd.read_csv(DATA/'Figure8_UT_coordinates_scores.csv')
orders=json.loads((DATA/'Figure8_category_order.json').read_text())
for dataset,key in [(allobs,'all'),(ut,'ut')]:
    for column,categories in orders[key].items():
        dataset[column]=pd.Categorical(dataset[column].astype(str),categories=categories)
palette=json.loads((DATA/'Figure8_palettes.json').read_text())
colors=palette['colors'];state_colors=palette['state_colors']
states=list(allobs.cellular_states.cat.categories)
tab=SC/'06_downstream_liana_v3/tables'
mit=pd.read_csv(source(tab/'mitotic_hub_gene_expression_by_ut_state.csv'))
clinical=pd.read_csv(source(tab/'clinical_bulk_and_proteomic_gene_detection_by_ut_state.csv'))
utstates=list(ut.cellular_states.cat.categories)

def compact_state(s):
    return s.replace(' cells','').replace('Cycling-','').replace('ECM-related (1)','ECM1').replace('ECM-related (2)','ECM2')

def umap(ax,data,key,continuous=False,title=None,size=.8):
    if continuous:
        order=data[key].sort_values().index
        im=ax.scatter(data.loc[order,'UMAP1'],data.loc[order,'UMAP2'],c=data.loc[order,key],s=size,cmap='Spectral_r',linewidths=0,rasterized=True)
    else:
        cats=list(data[key].cat.categories)
        palette=state_colors if key=='cellular_states' else dict(zip(cats,colors.get(key,plt.get_cmap('tab10').colors)))
        for cat in cats:
            dd=data[data[key]==cat];ax.scatter(dd.UMAP1,dd.UMAP2,s=size,color=palette[cat],linewidths=0,rasterized=True)
        im=None
    ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title or key,fontsize=7.5,pad=3)
    ax.set_aspect('equal',adjustable='box')
    ax.spines[['top','right']].set_visible(False);ax.spines[['bottom','left']].set_linewidth(.4)
    return im

def dot(ax,d,genes,labels=True,fs=7.5,scale=19):
    d=d[d.present_in_matrix.eq(True)].drop_duplicates(['gene','group'])
    genes=[g for g in genes if g in set(d.gene)]
    mean=d.pivot(index='group',columns='gene',values='mean_log_normalized_expression').reindex(index=utstates,columns=genes)
    frac=d.pivot(index='group',columns='gene',values='pct_detected').reindex(index=utstates,columns=genes)
    scaled=(mean-mean.min())/(mean.max()-mean.min());scaled=scaled.fillna(0)
    xx,yy=np.meshgrid(range(len(genes)),range(6))
    ax.scatter(xx.ravel(),yy.ravel(),s=(frac.to_numpy().ravel()/100)**1.5*scale,c=scaled.to_numpy().ravel(),vmin=0,vmax=1,cmap='viridis',linewidths=.25,edgecolors='#444444')
    ax.set_xlim(-.6,len(genes)-.4);ax.set_ylim(5.6,-.6);ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=fs);ax.set_yticks(range(6),[compact_state(s).replace('UT ','') for s in utstates] if labels else [],fontsize=fs)
    ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)
    return genes

mitgenes=list(dict.fromkeys(mit.gene))
cg=list(dict.fromkeys(clinical[clinical.present_in_matrix.eq(True)].gene))
fig=plt.figure(figsize=(170*MM,190*MM))
label(fig,'A',.005,.995)
for j,(k,t) in enumerate([('leiden','Leiden cluster'),('Treatment','Treatment'),('Batch','Batch'),('cellular_states','Cellular state')]):
    ax=fig.add_axes((.065+j*.235,.875,.185,.09));umap(ax,allobs,k,title=t,size=.45)
    if k=='leiden':
        for cat in allobs[k].cat.categories:
            dd=allobs[allobs[k]==cat];ax.text(dd.UMAP1.median(),dd.UMAP2.median(),cat,fontsize=6.8,ha='center',va='center',bbox={'facecolor':'white','edgecolor':'none','pad':.1,'alpha':.75})
    elif k in ['Treatment','Batch']:
        cats=list(allobs[k].cat.categories);cs=colors.get(k,plt.get_cmap('tab10').colors)
        handles=[Line2D([0],[0],marker='o',linestyle='',color=c,markersize=2) for c in cs[:len(cats)]]
        ax.legend(handles,[c.replace('npct','') for c in cats],loc='lower center',bbox_to_anchor=(.5,-.30),ncol=len(cats),frameon=False,fontsize=6.5,handlelength=.4,handletextpad=.1,columnspacing=.4,borderpad=0)
fig.text(.06,.842,'UMAP1 →; UMAP2 ↑. Batch labels omit the common “npct” prefix.',fontsize=7)
handles=[Line2D([0],[0],marker='o',linestyle='',color=state_colors[s],markersize=3) for s in states]
fig.legend(handles,[compact_state(s) for s in states],ncol=4,loc='upper left',bbox_to_anchor=(.035,.833),frameon=False,fontsize=7,labelspacing=.15,handlelength=.7,handletextpad=.3,columnspacing=1.1,borderpad=0)

label(fig,'B',.005,.755);ax=fig.add_axes((.075,.625,.18,.11))
comp=pd.crosstab(allobs.Treatment,allobs.cellular_states,normalize='index').mul(100).reindex(['UT','T24','T48']);bottom=np.zeros(3)
for s in states:
    ax.bar(range(3),comp[s],bottom=bottom,color=state_colors[s],width=.7);bottom+=comp[s]
ax.set_xticks(range(3),['UT','T24','T48']);ax.set_yticks([0,50,100]);ax.set_ylabel('Cells (%)');clean(ax)
label(fig,'C',.31,.755);ax=fig.add_axes((.50,.634,.43,.10));dot(ax,mit,mitgenes,scale=35)
ax.set_title('Mitotic markers',fontsize=8,pad=5)

label(fig,'D',.005,.58);ax=fig.add_axes((.225,.475,.735,.085));dot(ax,clinical,cg,scale=22)
fig.text(.225,.571,'Clinical bulk and proteomic candidates',fontsize=8)
# A shared key retains the source plot's within-gene 0–1 colour scaling.
ca=fig.add_axes((.76,.783,.15,.009));cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='viridis'),cax=ca,orientation='horizontal',ticks=[0,1]);cb.ax.tick_params(labelsize=6.5,length=2,pad=1)
fig.text(.74,.801,'Mean expression (scaled)',fontsize=6.8)
ka=fig.add_axes((.73,.813,.26,.016));ka.axis('off')
for x,pct in zip([.10,.35,.62],[20,60,100]):ka.scatter(x,.5,s=(pct/100)**1.5*25,c='gray');ka.text(x+.04,.5,str(pct),va='center',fontsize=6.5)
ka.set_xlim(0,1);ka.set_ylim(0,1);fig.text(.74,.834,'Cells expressing (%)',fontsize=6.8)

label(fig,'E',.005,.42);label(fig,'F',.285,.42)
keys=[('pseudotime','Pseudotime'),('cellular_states','UT state'),('S_score','S phase'),('G2M_score','G2/M'),('Inflammatory_Score','Inflammatory'),('EMT_score','EMT'),('Epithelial_score','Epithelial')]
for j,(k,t) in enumerate(keys):
    x=.045+j*.136;ax=fig.add_axes((x,.32,.108,.075));im=umap(ax,ut,k,k!='cellular_states',t,size=.5)
    if im:
        ca=fig.add_axes((x,.304,.108,.006));cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[float(ut[k].min()),float(ut[k].max())]);cb.ax.tick_params(labelsize=6.5,pad=1,length=1.5);cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v,p:f'{v:.2g}'))
        ticks=cb.ax.get_xticklabels();ticks[0].set_ha('left');ticks[-1].set_ha('right')

label(fig,'G',.005,.278)
enrichments=[]
for j,direction in enumerate(['increasing','decreasing']):
    path=SC/f'05_downstream_enrichment_v3/enrichment/pseudotime_{direction}/pseudotime_{direction}_all.csv'
    full=pd.read_csv(source(path))
    # Match the exact ten terms in draft13, including the original tie selection.
    ids=([ 'R-HSA-1640170','R-HSA-69278','R-HSA-69620','R-HSA-68886','R-HSA-68877',
            'R-HSA-2555396','R-HSA-2500257','R-HSA-68882','GO:0000070','R-HSA-2467813'] if j==0 else
         ['R-HSA-168249','R-HSA-6798695','R-HSA-156827','R-HSA-168256','R-HSA-72737',
          'R-HSA-72706','R-HSA-72689','R-HSA-6809371','R-HSA-72649','R-HSA-72662'])
    d=pd.DataFrame([full[full.Term.str.contains(term,regex=False)].iloc[0] for term in ids])
    labels=[]
    for s in d.Term:
        s=re.sub(r'\s*R-HSA-\d+','',s);s=re.sub(r'\s*\(GO:\d+\)','',s)
        replacements={'L13a-mediated Translational Silencing Of Ceruloplasmin Expression':'L13a-mediated\nceruloplasmin silencing',
        'mRNA Activation Upon Binding Of Cap-Binding Complex And eIFs, Subsequent Binding To 43S':'mRNA activation /\n43S recruitment',
        'Cap-dependent Translation Initiation':'Cap-dependent\ntranslation initiation',
        'GTP Hydrolysis And Joining Of 60S Ribosomal Subunit':'GTP hydrolysis / 60S joining',
        'Formation Of A Pool Of Free 40S Subunits':'Free 40S subunit pool',
        'Formation Of Cornified Envelope':'Cornified envelope formation',
        'Translation Initiation Complex Formation':'Translation initiation complex',
        'Resolution Of Sister Chromatid Cohesion':'Sister chromatid\ncohesion resolution',
        'Mitotic Sister Chromatid Segregation':'Mitotic sister\nchromatid segregation',
        'Separation Of Sister Chromatids':'Sister chromatid separation',
        'Mitotic Metaphase And Anaphase':'Mitotic metaphase / anaphase'}
        s=replacements.get(s,s)
        labels.append(s)
    d['display_label']=labels;d.to_csv(DATA/f'Figure8G_{direction}_displayed.csv',index=False);enrichments.append(d)
    ax=fig.add_axes((.345+j*.50,.045,.14,.205));ax.barh(range(10),-np.log10(d['Adjusted P-value']),color=plt.get_cmap('Reds' if j==0 else 'Blues')(np.linspace(.85,.25,10)),height=.7)
    ax.set_yticks(range(10),labels,fontsize=6.9);ax.invert_yaxis();ax.set_xlabel('−log10 adj. P',fontsize=7);ax.xaxis.set_major_locator(plt.MaxNLocator(3));ax.tick_params(axis='x',labelsize=7);clean(ax)
    fig.text(.06+j*.50,.263,'Pseudotime: '+direction,fontsize=8,weight='bold')
save(fig,8)

# Expanded panels preserve full labels and provide readable review/export masters.
for label_,data,genes in [('C',mit,mitgenes),('D',clinical,cg)]:
    f,ax=plt.subplots(figsize=(170*MM,65*MM));f.subplots_adjust(left=.21,right=.98,bottom=.32,top=.91)
    dot(ax,data,genes,fs=9,scale=45);ax.set_title('Untreated NPC43: '+('mitotic markers' if label_=='C' else 'clinical and proteomic candidates'),fontsize=10)
    f.savefig(OUT/f'panels/Figure8{label_}_expanded.pdf');plt.close(f)
for direction,d in zip(['increasing','decreasing'],enrichments):
    f,ax=plt.subplots(figsize=(170*MM,100*MM));f.subplots_adjust(left=.65,right=.98,bottom=.14,top=.9)
    ax.barh(range(10),-np.log10(d['Adjusted P-value']),color='#A74B4B' if direction=='increasing' else '#4278A5');ax.set_yticks(range(10),[textwrap.fill(s,55) for s in d.Term],fontsize=8);ax.invert_yaxis();ax.set_xlabel('−log10 adjusted P value');f.suptitle('Genes '+direction+' along inferred pseudotime',fontsize=10)
    f.savefig(OUT/f'panels/Figure8G_{direction}_expanded.pdf');plt.close(f)
(DATA/'singlecell_replot_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print('Figure 8 complete',flush=True)
