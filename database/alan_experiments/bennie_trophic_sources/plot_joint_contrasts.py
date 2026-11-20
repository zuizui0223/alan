"""One experiment, six dependent contrasts; no pooled diamond."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE.parents[2]/'manuscript/biotic_synthesis/figures'
OUT.mkdir(exist_ok=True)
d=pd.read_csv(HERE/'joint_predator_contrasts.csv')
fig,ax=plt.subplots(figsize=(10.6,5.1))
fig.subplots_adjust(left=.34,right=.73,top=.79,bottom=.24)
labels={'aphid_count':'Aphids | final September 2013','slug_count':'Slugs | 15-month mean, 2013-2014','plant_biomass':'Plants | final biomass, 2016'}
colors={'Amber':'#bd7800','White':'#4773ad'}
for i,r in enumerate(d.itertuples()):
    y=5-i
    ax.errorbar(r.ratio_of_response_ratios,y,xerr=[[r.ratio_of_response_ratios-r.approx_t5_lower],[r.approx_t5_upper-r.ratio_of_response_ratios]],fmt='o',color=colors[r.light],capsize=3,markersize=6,lw=1.5)
    ax.text(1.035,y,f'{r.ratio_of_response_ratios:.2f} [{r.approx_t5_lower:.2f}, {r.approx_t5_upper:.2f}]',transform=ax.get_yaxis_transform(),va='center',fontsize=9)
ax.set_yticks(range(5,-1,-1),[f'{labels[r.endpoint]}\n{r.light} light' for r in d.itertuples()],fontsize=9)
ax.set_xscale('log');ax.set_xlim(.01,100);ax.set_ylim(-.7,5.7)
ax.axvline(1,color='.45',linestyle='--',lw=1,zorder=0)
ax.set_xticks([.01,.1,1,10,100],['0.01','0.1','1','10','100'])
ax.set_xlabel('Ratio of lighting response ratios\n(predators present / herbivores only)',fontsize=10)
for y in [1.5,3.5]:ax.axhline(y,color='.9',lw=.8,zorder=0)
for spine in ['top','right','left']:ax.spines[spine].set_visible(False)
ax.tick_params(axis='y',length=0,pad=9)
fig.suptitle('Does predator presence modify the response to lighting?',x=.05,ha='left',fontsize=14,y=.96)
fig.text(.05,.875,'Penryn experiment: 6 mesocosms per treatment cell; one shared experimental programme',fontsize=10)
fig.text(.05,.075,'Points: direct within-experiment contrasts. Bars: approximate delta-method 95% intervals using t(5).\nResponse periods differ; these are not a simultaneous food-chain or mediation analysis. No pooled estimate.',fontsize=9,color='.3')
for ext in ['png','svg','pdf']:fig.savefig(OUT/f'penryn_predator_contrasts.{ext}',dpi=200,facecolor='white')
svg=OUT/'penryn_predator_contrasts.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
plt.close(fig)
