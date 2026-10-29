#!/usr/bin/env python3
"""Cobertura mensual de pozos por region (2020 en adelante), v1.1 vs v1.3.

Lee ../versions/{v1.1,v1.3}/*.csv, asigna cada pozo a una region por cruce
espacial con input/other_data/Chile_boundaries_simplified_from_SIIT.gpkg
(poligono mas cercano si cae fuera) y escribe el PNG en
../comparison_v1.1_vs_v1.3/. Requiere pandas, geopandas y matplotlib.
"""
import pandas as pd, geopandas as gpd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.dates as md
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
W=Path(__file__).resolve().parent.parent.as_posix()+'/'; R=Path(W).parent.as_posix()+'/'
def load(v):
    a=pd.read_csv(W+f'versions/{v}/cr2sub_v1.1_attributes.csv',dtype={'cr2sub_id':str})
    g=pd.read_csv(W+f'versions/{v}/cr2sub_v1.1_gwl_mon.csv',parse_dates=['date']).set_index('date')
    g.columns=g.columns.astype(str); return a,g[g.index>='2020-01-01'].notna()
a1,h1=load('v1.1'); a3,h3=load('v1.3')
att=pd.concat([a1,a3]).drop_duplicates('cr2sub_id')
reg=gpd.read_file(R+'input/other_data/Chile_boundaries_simplified_from_SIIT.gpkg')
pts=gpd.GeoDataFrame(att[['cr2sub_id','cr2sub_lat']],geometry=gpd.points_from_xy(att.cr2sub_lon,att.cr2sub_lat),crs=4326)
j=gpd.sjoin_nearest(pts.to_crs(32719),reg[['Region','geometry']].to_crs(32719),how='left').drop_duplicates('cr2sub_id')
j['name']=j.Region.str.replace('Región ','').str.replace('del ','').str.replace('de ','').str.replace('Libertador General Bernardo O\'Higgins',"O'Higgins").str.replace('Metropolitana Santiago','Metropolitana').str.replace('Libertador Bernardo O\'Higgins',"O'Higgins")
rid=j.set_index('cr2sub_id').name
lat=j.groupby('name').cr2sub_lat.mean().sort_values(ascending=False); order=list(lat.index)
def by_reg(h):
    D=pd.DataFrame({n:h[[c for c in h.columns if rid.get(c)==n]].sum(axis=1) for n in order}); return D
D1,D3=by_reg(h1),by_reg(h3)
# v1.1 solo llega a 2024-12: reindexar al rango completo
idx=D3.index; D1=D1.reindex(idx)
nw=j.groupby('name').size()
cm=LinearSegmentedColormap.from_list('s',['#EEF3F9','#C3D5EA','#8FB1D8','#5A8BC2','#2E64A5','#173F73'])
vmax=max(np.nanmax(D1.values),D3.values.max())
INK,MUT='#1F2937','#6B7280'; C1,C3='#E8833A','#2E64A5'
fig=plt.figure(figsize=(9.5,12))
gs=fig.add_gridspec(3,2,height_ratios=[1.1,5,5],width_ratios=[1,0.03],hspace=0.3,wspace=0.03)
at=fig.add_subplot(gs[0,0]); a=fig.add_subplot(gs[1,0]); b=fig.add_subplot(gs[2,0],sharex=a); cax=fig.add_subplot(gs[1:,1])
t1,t3=h1.sum(axis=1).reindex(idx),h3.sum(axis=1)
at.step(idx,t3,where='mid',color=C3,lw=4.5,label='v1.3',zorder=2); at.step(idx,t1,where='mid',color=C1,lw=1.6,label='v1.1 (superpuesta)',zorder=3)
at.axvline(pd.Timestamp('2025-01-01'),color=MUT,lw=0.8,ls=(0,(3,3)),zorder=1); at.text(pd.Timestamp('2025-01-20'),at.get_ylim()[1]*0.95,'v1.1 termina\n(dic 2024)',fontsize=8,color=MUT,va='top')
at.set_ylabel('Pozos\ncon dato',fontsize=9,color=MUT); at.legend(frameon=False,fontsize=9,loc='upper left',ncol=2,bbox_to_anchor=(0.0,1.35))
at.set_xlim(idx[0],idx[-1]); at.xaxis.set_major_locator(md.YearLocator(1)); at.xaxis.set_major_formatter(md.DateFormatter('%Y'))
for ax,D,t in ((a,D1,'v1.1  (input: 844 reportes xls)'),(b,D3,'v1.3  (input: export único)')):
    M=D.T.values.astype(float); M[M==0]=np.nan
    ext=[md.date2num(idx[0]),md.date2num(idx[-1])+31,len(order),0]
    im=ax.imshow(M,aspect='auto',cmap=cm,extent=ext,interpolation='nearest',vmin=0,vmax=vmax)
    ax.xaxis_date(); ax.xaxis.set_major_locator(md.YearLocator(1)); ax.xaxis.set_major_formatter(md.DateFormatter('%Y'))
    ax.set_title(t,loc='left',fontsize=10.5,fontweight='bold',color=INK,pad=6)
    ax.tick_params(colors=MUT,labelsize=8.5)
for ax in (a,b):
    ax.set_yticks(np.arange(len(order))+.5); ax.set_yticklabels([f'{n}  ({nw[n]})' for n in order],fontsize=9,color=INK)
for ax in (at,a,b):
    for s in ('top','right'): ax.spines[s].set_visible(False)
at.tick_params(colors=MUT,labelsize=8.5)
cb=fig.colorbar(im,cax=cax); cb.set_label('Pozos con dato en el mes',fontsize=9,color=MUT); cb.ax.tick_params(labelsize=8,colors=MUT)
fig.suptitle('CR2SUB — cobertura mensual de pozos por región, v1.1 vs v1.3 (2020–2026)',x=0.08,ha='left',fontsize=12,fontweight='bold',color=INK,y=0.955)
out=W+'comparison_v1.1_vs_v1.3/cobertura_mensual_por_region_v1.1_vs_v1.3_2020-2026.png'
fig.savefig(out,dpi=130,bbox_inches='tight',pad_inches=0.3,facecolor='white')
d=(D3-D1.fillna(0)); print('meses v1.1:',idx[D1.notna().any(axis=1)][[0,-1]].tolist())
print('pozos-mes extra por región:'); print(d.clip(lower=0).sum().astype(int).to_string()); print('perdidos',(-d.clip(upper=0)).sum().sum())
