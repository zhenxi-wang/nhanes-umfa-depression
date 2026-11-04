from analysis_runtime import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm,ticker
from matplotlib.lines import Line2D
from PIL import Image
FONT_FAMILY=os.environ.get('NHANES_FONT_FAMILY','Arial')
fontdir=os.environ.get('NHANES_FONT_DIR')
if fontdir is None and os.name=='nt':fontdir=str(Path(os.environ.get('WINDIR','C:/Windows'))/'Fonts')
if fontdir:
 for name in ['arial.ttf','arialbd.ttf','ariali.ttf','arialbi.ttf']:
  fp=Path(fontdir)/name
  if fp.exists():fm.fontManager.addfont(str(fp))
try:fm.findfont(FONT_FAMILY,fallback_to_default=False)
except ValueError as exc:raise RuntimeError('Requested font not found. Set NHANES_FONT_DIR or NHANES_FONT_FAMILY; see README.md.') from exc
plt.rcParams.update({'font.family':FONT_FAMILY,'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65,'xtick.major.size':3,'ytick.major.size':3,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
NAVY='#18354A';TEAL='#16857D';GREY='#687782';LIGHT='#E8EDF0';RED='#B45A48';INK='#18242B'
core=json.loads((OUT/'core_results.json').read_text());cm=json.loads((OUT/'curves_margins_validation.json').read_text());rep=json.loads((OUT/'replication_sensitivity_results.json').read_text());a=pd.read_pickle(OUT/'analytic_base.pkl')
figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
def save(fig,name):
 fig.canvas.draw()
 for ext in ['pdf','svg']:fig.savefig(figdir/(name+'.'+ext),facecolor='white')
 fig.savefig(figdir/(name+'.png'),dpi=1000,facecolor='white')
 with Image.open(figdir/(name+'.png')) as im:im.save(figdir/(name+'.tif'),compression='tiff_lzw',dpi=(1000,1000))
 fig.savefig(figdir/(name+'_preview.png'),dpi=180,facecolor='white')
 plt.close(fig);print('SAVED',name,flush=True)
def letter(fig,ax,t,x=-.09):ax.text(x,1.055,t,transform=ax.transAxes,fontsize=12,fontweight='bold',ha='left',va='bottom',color=INK)
# Figure 1: the three adjusted models are all recomputed using MICE.
fig=plt.figure(figsize=(7.2,4.15));gs=fig.add_gridspec(1,2,left=.17,right=.985,bottom=.23,top=.87,width_ratios=[1.4,1],wspace=.42)
ax=fig.add_subplot(gs[0]);labs=['5-MTHF','Total serum folate','Red-cell folate','MeFox','UMFA','UMFA fraction'];yy=np.arange(6)[::-1]
colors={1:'#A0A9AE',2:TEAL,3:NAVY};markers={1:'D',2:'s',3:'o'}
ax.axvline(1,color='#9BA5AA',lw=.8,ls=(0,(3,3)),zorder=0)
for mi,mo in enumerate([1,2,3]):
 for j,label in enumerate(labs):
  r=core[f'single_{label}_M{mo}'];y=yy[j]+(.22*(mi-1));ax.plot([r['lo'],r['hi']],[y,y],color=colors[mo],lw=1.2,solid_capstyle='butt');ax.plot(r['est'],y,markers[mo],color=colors[mo],ms=4,mec='white',mew=.55)
ax.set_yticks(yy,labs);ax.tick_params(axis='y',length=0,pad=7);ax.set_xscale('log');ax.set_xlim(.69,1.32);ax.set_xticks([.7,.8,.9,1,1.1,1.2,1.3],['0.7','0.8','0.9','1.0','1.1','1.2','1.3']);ax.xaxis.set_minor_locator(ticker.NullLocator());ax.set_ylim(-.6,5.6);ax.set_xlabel('Odds ratio per SD of log biomarker',labelpad=7);letter(fig,ax,'a',-.38)
fig.legend(handles=[Line2D([],[],marker=markers[m],color=colors[m],lw=1.2,ms=4,label=f'Model {m}') for m in [1,2,3]],loc='lower left',bbox_to_anchor=(.17,.065),frameon=False,ncol=3,columnspacing=1.3,handlelength=1.4,fontsize=8)
ax=fig.add_subplot(gs[1]);des=c.design(a);rawprev=[]
for lv,label,color,marker in [('l_mthf','5-MTHF',TEAL,'s'),('l_umfa','UMFA',RED,'o')]:
 vals=[]
 for q in [1,2,3,4]:
  pr=des.subset((a[lv+'_q']==q).to_numpy()).prop(a.dep.to_numpy());vals.append([pr[x]*100 for x in ['est','lo','hi']]);rawprev.append(dict(biomarker=label,quartile=q,**pr))
 vals=np.array(vals);ax.errorbar([1,2,3,4],vals[:,0],yerr=[vals[:,0]-vals[:,1],vals[:,2]-vals[:,0]],fmt=marker+'-',color=color,lw=1.25,ms=4.5,capsize=2.2,capthick=.8,label=label)
ax.set_xticks([1,2,3,4],['Q1','Q2','Q3','Q4']);ax.set_xlim(.8,4.2);ax.set_ylim(5.5,15);ax.set_yticks([6,8,10,12,14]);ax.set_ylabel('Unadjusted prevalence (%)',labelpad=5);ax.set_xlabel('Cycle-specific biomarker quartile',labelpad=7);ax.legend(frameon=False,loc='upper right',fontsize=8,handlelength=1.6);letter(fig,ax,'b',-.24)
fig.text(.17,.015,'Points and bars show estimates and 95% confidence intervals.',fontsize=7.3,color=GREY)
save(fig,'Figure_1_MICE');dump('figure1_unadjusted_prevalences.json',rawprev)
# Figure 2: both panels use the identical mutually adjusted spline model.
fig=plt.figure(figsize=(7.2,4.4));gs=fig.add_gridspec(2,2,left=.1,right=.985,bottom=.21,top=.9,width_ratios=[1.3,1],height_ratios=[3.2,1],wspace=.44,hspace=.13)
ax=fig.add_subplot(gs[0,0]);histax=fig.add_subplot(gs[1,0],sharex=ax);cv=pd.DataFrame(cm['curve'])
ax.fill_between(cv.umfa,cv.lo,cv.hi,color=NAVY,alpha=.14,lw=0);ax.plot(cv.umfa,cv.est,color=NAVY,lw=1.65);ax.axhline(1,color='#9BA5AA',lw=.75,ls=(0,(3,3)));ax.set_xscale('log');ax.set_xlim(.16,32);ax.set_ylim(.65,max(2.35,float(cv.hi.max())*1.03));ax.set_ylabel('Adjusted odds ratio (95% CI)',labelpad=6);ax.tick_params(labelbottom=False);ax.xaxis.set_minor_locator(ticker.NullLocator());letter(fig,ax,'a',-.18)
ax.plot(cm['reference'],1,'o',ms=3.5,mec=NAVY,mfc='white',mew=.8,zorder=5);ax.text(.04,.91,'Reference: 0.35 nmol/L',transform=ax.transAxes,fontsize=7.8,color=GREY);ax.text(.04,.82,f"P for non-linearity = {cm['nonlinearity']['p']:.3f}",transform=ax.transAxes,fontsize=7.8,color=GREY)
histax.hist(a.umfa,bins=np.logspace(np.log10(.16),np.log10(32),62),color='#AAB5BA',lw=0);histax.set_xscale('log');histax.set_xticks([.2,.5,1,2,4,8,20],['0.2','0.5','1','2','4','8','20']);histax.xaxis.set_minor_locator(ticker.NullLocator());histax.set_yticks([]);histax.spines['left'].set_visible(False);histax.set_xlabel('Serum UMFA (nmol/L; log scale)',labelpad=7)
histax.text(.98,.88,'>4 nmol/L: n = 1,166\n>8 nmol/L: n = 721',transform=histax.transAxes,ha='right',va='top',fontsize=7.0,color=GREY)
ax=fig.add_subplot(gs[:,1]);pp=pd.DataFrame(cm['prevalence']);ax.fill_between(pp.percentile,pp.lo,pp.hi,color=NAVY,alpha=.14,lw=0);ax.plot(pp.percentile,pp.prev,'o-',color=NAVY,lw=1.5,ms=4.1);ax.set_xlim(5,104);ax.set_xticks([10,25,50,75,90,99]);ax.set_ylim(6.8,max(17.6,float(pp.hi.max())+.3));ax.set_yticks([8,10,12,14,16]);ax.set_ylabel('Model-adjusted prevalence (%)',labelpad=5);ax.set_xlabel('Percentile of serum UMFA',labelpad=7);letter(fig,ax,'b',-.24)
for j in [0,-1]:
 r=pp.iloc[j];ax.annotate(f"{r.prev:.1f}%",(r.percentile,r.prev),xytext=(2,-13 if j==0 else 8),textcoords='offset points',ha='center',fontsize=8,fontweight='bold',color=NAVY)
rd=cm['risk_difference'];fig.text(.1,.06,f"10th–99th percentile difference: {rd['est']:.1f} percentage points (95% CI {rd['lo']:.1f}–{rd['hi']:.1f})",fontsize=8,color=GREY);fig.text(.1,.02,'Both panels: three-species Model 3, 30 imputations; n = 24,326.',fontsize=7.3,color=GREY)
save(fig,'Figure_2_MICE')
# Figure 3: separate imputations within discovery, replication, and each cycle.
fig=plt.figure(figsize=(7.2,4.35));ax=fig.add_axes([.45,.19,.25,.67]);rows=[('Discovery: Aug 2021–Aug 2023',rep['L'],True),('Earlier-cycle replication',rep['earlier'],True),('Pooled: 2011–2023',core['joint_M3']['UMFA'],True)]
rows += [('2011–2012',rep['G'],False),('2013–2014',rep['H'],False),('2015–2016',rep['I'],False),('2017–March 2020',rep['P'],False),('Aug 2021–Aug 2023',rep['L'],False)]
y=[8,7,6,4,3,2,1,0];ax.axvline(1,color='#9BA5AA',ls=(0,(3,3)),lw=.8);ax.axhspan(5.6,6.4,color='#F0F3F5',zorder=0)
for (lab,r,main),yy in zip(rows,y):
 col=NAVY if main else '#7E8C94';ax.plot([r['lo'],r['hi']],[yy,yy],color=col,lw=1.5 if main else 1.05);ax.plot(r['est'],yy,'D' if main else 'o',color=col,ms=5.1 if main else 4.1,mec='white',mew=.55)
 fig.text(.018,.19+.67*(yy+.55)/9.1,lab,fontsize=8.1,fontweight='bold' if main else 'normal',va='center',color=INK)
 fig.text(.752,.19+.67*(yy+.55)/9.1,f"{r['n']:,}",fontsize=8,ha='right',va='center',color=INK)
 fig.text(.787,.19+.67*(yy+.55)/9.1,f"{r['est']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})",fontsize=8,va='center',color=INK)
ax.set_ylim(-.55,8.55);ax.set_yticks([]);ax.spines['left'].set_visible(False);ax.set_xscale('log');ax.set_xlim(.8,1.75);ax.set_xticks([.8,1,1.25,1.5],['0.8','1.0','1.25','1.5']);ax.xaxis.set_minor_locator(ticker.NullLocator());ax.set_xlabel('Adjusted odds ratio',labelpad=6)
fig.text(.018,.91,'Analysis sample',fontsize=8.5,fontweight='bold');fig.text(.752,.91,'n',ha='right',fontsize=8.5,fontweight='bold');fig.text(.787,.91,'OR (95% CI)',fontsize=8.5,fontweight='bold');fig.text(.018,.19+.67*(5+.55)/9.1,'Individual survey cycles',color=GREY,fontsize=8,fontstyle='italic')
h=rep['heterogeneity'];fig.text(.018,.065,f"Heterogeneity: Q = {h['Q']:.2f} (df = 4), P = {h['p']:.3f}; I² = {h['I2']:.0f}%",fontsize=8,color=GREY);fig.text(.018,.025,'Three-species Model 3; OR per within-cycle SD of log UMFA; 30 imputations per analysis sample.',fontsize=7.2,color=GREY)
save(fig,'Figure_3_MICE')
