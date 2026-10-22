##### Packages
import sys,os

import xarray as xr
import pandas as pd
import numpy as np

import NSSEA as ns
import NSSEA.models as nsm
import scipy.stats as sc
from scipy.stats import gaussian_kde
import statsmodels.api as sm

from itertools import product
import string

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatch
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
from matplotlib.lines import Line2D
import matplotlib.backends.backend_pdf as mpdf

import seaborn as sns

#import matplotlib.colors as mcolors

#







### Added functions 
def weighted_median(values, weights):
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    return values[i[np.searchsorted(c, 0.5 * c[-1])]]


def header(label):
    return Line2D([], [], color="none", label=f"$\\bf{{{label}}}$")
plt.rcParams.update({
    # Font
    "font.size":          8,
    "axes.titlesize":     9,
    "axes.labelsize":     8,
    "xtick.labelsize":    7,
    "ytick.labelsize":    7,
    "legend.fontsize":    7,
    "legend.title_fontsize": 7,

    # Figure
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",

})

def _hist_kde(ax, values, weights=None, color="#C0392B", n_pts=512,
              bins=80, lw=1.0, alpha_bar=0.30, label=None):
    """
    Draw a step histogram (no bar edges) with a KDE curve on top.
    The histogram uses 'count' as stat; the KDE is scaled to match.
    Works on whichever axis is passed (left or right twinx).
    """
    w = None
    if weights is not None:
        w = np.asarray(weights, dtype=float)

    # Histogram — step style so bars are outlined, not filled solid
    n, bin_edges, patches = ax.hist(
        values, bins=bins, weights=w,
        histtype="stepfilled",
        color=color, alpha=alpha_bar,
        linewidth=0, label=label,
    )
    ax.hist(
        values, bins=bin_edges, weights=w,
        histtype="step",
        color=color, linewidth=0.6, alpha=0.7,
    )

    # KDE scaled to match count histogram
    bw = "scott"
    if w is not None:
        w_norm = w / w.sum() * len(w)
        kde = gaussian_kde(values, bw_method=bw, weights=w_norm)
    else:
        kde = gaussian_kde(values, bw_method=bw)

    bin_width = bin_edges[1] - bin_edges[0]
    total     = w.sum() if w is not None else len(values)

    x_min, x_max = values.min(), values.max()
    pad = (x_max - x_min) * 0.05
    xs  = np.linspace(x_min - pad, x_max + pad, n_pts)
    ys  = kde(xs) * total * bin_width   # scale to count

    ax.plot(xs, ys, color=color, lw=lw, zorder=5)
    return xs, ys


def _vline_with_label(ax, x, color, side="left", y_transform=None):
    """Draw a vertical dashed median line and annotate below the x-axis."""
    ax.axvline(x, color=color, ls="--", lw=0.8, zorder=4)
    ha  = "right" if side == "left" else "left"
    off = -0.003 if side == "left" else +0.003
    ax.text(
        x + off, -0.04, f"{x:.2f}",
        color=color, ha=ha, va="top", fontsize=5.5,
        transform=ax.get_xaxis_transform(),
    )


def _gcm_dots(ax, gcm, y_frac=0.04):
    """
    Scatter GCM estimates just above the x-axis on *ax* (should be the
    topmost twinx so dots sit above both KDEs).  Uses ax's own xlim.
    """
    x_min, x_max = ax.get_xlim()
    inside = gcm[(gcm >= x_min) & (gcm <= x_max)]
    ax.scatter(
        inside, np.full_like(inside, y_frac),
        transform=ax.get_xaxis_transform(),
        marker="|", s=18, color=BLACK, lw=0.7, zorder=10,
        clip_on=False,
    )
    if np.any(gcm < x_min):
        ax.annotate("", xy=(0.01, y_frac), xycoords=ax.transAxes,
                    xytext=(0.045, y_frac),
                    arrowprops=dict(arrowstyle="<-", color=BLACK, lw=0.7))
    if np.any(gcm > x_max):
        ax.annotate("", xy=(0.99, y_frac), xycoords=ax.transAxes,
                    xytext=(0.955, y_frac),
                    arrowprops=dict(arrowstyle="<-", color=BLACK, lw=0.7))


def _style_axes(ax_l, ax_r, xlabel=""):
    """Remove spines, set tick counts, colour left/right axes."""
    for ax in (ax_l, ax_r):
        ax.spines["top"].set_visible(False)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        ax.set_xlabel(xlabel, fontsize=6)

    # Left (red)
    ax_l.spines["left"].set_color(RED)
    ax_l.spines["right"].set_visible(False)
    ax_l.tick_params(axis="y", colors=RED, length=2.5, width=0.6)
    ax_l.set_ylabel("Count \n (all samples)", color=RED, labelpad=3)
    ax_l.yaxis.label.set_color(RED)
    ax_l.xaxis.grid(False)
    ax_l.yaxis.grid(False)

    # Right (blue)
    ax_r.spines["right"].set_color(BLUE)
    ax_r.spines["left"].set_visible(False)
    ax_r.tick_params(axis="y", colors=BLUE, length=2.5, width=0.6)
    ax_r.set_ylabel("Count \n (weighted subset)", color=BLUE, labelpad=3)
    ax_r.yaxis.label.set_color(BLUE)
    ax_r.xaxis.grid(False)
    ax_r.yaxis.grid(False)
def make_latex_table(summary, scenarios, years, caption="", label=""):
    ncol_per_scen = 2  # BE + 5-95%
    col_spec = "l" + "cc" * len(scenarios)

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{" + col_spec + "}")
    lines.append(r"\toprule")

    # Ligne d'en-tête 1 : noms des scénarios en multicolonnes
    header1 = ["Year"]
    for scen in scenarios:
        header1.append(r"\multicolumn{2}{c}{\textbf{" + scenario_labels.get(scen, scen) + "}}")
    lines.append(" & ".join(header1) + r" \\")

    # cmidrule sous chaque scénario
    cmidrules = []
    start = 2
    for _ in scenarios:
        end = start + 1
        cmidrules.append(rf"\cmidrule(lr){{{start}-{end}}}")
        start = end + 1
    lines.append("".join(cmidrules))

    # Ligne d'en-tête 2 : BE / 5-95%
    header2 = [""]
    for _ in scenarios:
        header2.extend(["BE", "2.5-97.5\\%"])
    lines.append(" & ".join(header2) + r" \\")
    lines.append(r"\midrule")

    # Corps du tableau
    for year in years:
        row_vals = [str(year)]
        for scen in scenarios:
            be = summary.loc[year, (scen, 'BE')]
            interval = summary.loc[year, (scen, '5-95\\%')]
            row_vals.append(f"{be:.2f}")
            row_vals.append(interval)
        lines.append(" & ".join(row_vals) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    if caption:
        lines.append(r"\caption{" + caption + "}")
    if label:
        lines.append(r"\label{" + label + "}")
    lines.append(r"\end{table}")

    return "\n".join(lines)
    
    
def make_latex_table_styled(summary, scenarios, row_labels, scenario_labels,
                              caption="", label=""):
    col_spec = "l" + "c" * len(scenarios)

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{" + col_spec + "}")
    lines.append(r"\toprule")

    header = ["Period"] + [r"\textbf{" + scenario_labels.get(s, s) + "}" for s in scenarios]
    lines.append(" & ".join(header) + r" \\")
    lines.append(r"\midrule")

    for row_label in row_labels:
        row_vals = [row_label]
        for scen in scenarios:
            val = summary.loc[row_label, scen]
            row_vals.append(f"{val:.2f}")
        lines.append(" & ".join(row_vals) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    if caption:
        lines.append(r"\caption{" + caption + "}")
    if label:
        lines.append(r"\label{" + label + "}")
    lines.append(r"\end{table}")

    return "\n".join(lines)

        
#External Arguments

Y_ref=sys.argv[1]
X_ref=sys.argv[2]

work_path=sys.argv[3]

pathObs=sys.argv[4]
pathGCM_X_separated=sys.argv[5]
pathGCM_Y_separated=sys.argv[6]

dt_string=sys.argv[7]


##### Paths
print("Figures creation for: "+Y_ref)

basepath=work_path
pathFig=os.path.join(basepath,'Figures/'+Y_ref+"_"+X_ref,dt_string)
if not os.path.exists(pathFig):
    os.makedirs(pathFig)
print("Working path: "+pathFig)

pathOut=os.path.join(basepath,'Outputs/'+Y_ref+"_"+X_ref,dt_string)
if not os.path.exists(pathOut):
    os.makedirs(pathOut)
    

### Fixed arguments
scenarios_all=["ssp119","ssp126","ssp245","ssp370","ssp585"]
sample_pred=10
ci          = 0.05
ns_law      = nsm.GEV()

#Time period of interest
deb=1850
fin=2101
time_period    = np.arange( 1850 , 2101 , 1 , dtype = int )
time_year=list(range(1850, 2101))

seuil=61.8 #For extrem dist
scen_seuil="ssp585"

T_RL=[2,100,1000] #For RL
prob_interest=[0.5,0.25,0.1,0.05,0.01] #For pred


sub_scenarios_figures=["ssp126","ssp245","ssp585"]
sub_scenarios_tables=["ssp126","ssp245","ssp370","ssp585"]

#color_EQR={"ssp119":'#009999',"ssp126":'#003d99',"ssp245":'#997a00',"ssp370":'#990000',"ssp585":'#4d0000'}
#color_RL={"ssp119":"#00c3ff","ssp126":"#0066ff","ssp245":"#ffcc00","ssp370":"#ff5050","ssp585":"#800000"}


    
SCENARIO_META = {
    "ssp126": dict(label="SSP1-2.6", color="#2166ac", lw=1.4, ls="-"),
    "ssp245": dict(label="SSP2-4.5", color="#d6a024", lw=1.4, ls="-"),
    "ssp585": dict(label="SSP5-8.5", color="#b2182b", lw=1.4, ls="-"),
}
scenario_labels = {
    'ssp119': 'SSP1-1.9',
    'ssp126': 'SSP1-2.6',
    'ssp245': 'SSP2-4.5',
    'ssp370': 'SSP3-7.0',
    'ssp585': 'SSP5-8.5',
}



#Observations Meteo
dXo = xr.open_dataset(os.path.join( pathObs,"Xo_Ano.nc")) #Deja en anomalies
Xo  = pd.DataFrame( dXo.Xo.values.squeeze() , columns = ["Xo"] , index = dXo.index.values )

dYo =  xr.open_dataset(os.path.join( pathObs,"Yo_Ano.nc"))
Yo  =  pd.DataFrame( dYo.Yo.values.squeeze() , columns = ["Yo"] , index = dYo.index.values )

dbias_multi=xr.open_dataset(os.path.join( pathObs,"bias_Obs.nc"))
bias_multi=float(dbias_multi["Multi_Synthesis"].values.squeeze())



climCXCB=ns.Climatology.from_netcdf(os.path.join( pathOut , ("clim_multiscenarios_CXCB_"+dt_string+".nc")  )  , ns_law  )
clim=ns.Climatology.from_netcdf( os.path.join( pathOut , ("clim_multiscenarios_multisynthesis_allmodels_"+dt_string+".nc"))  , ns_law   )
clim_MM=ns.Climatology.from_netcdf( os.path.join( pathOut , ("clim_multiscenarios_multisynthesis_"+dt_string+".nc"))  , ns_law   )

samples_MCMC = climCXCB.data.sample_MCMC
n_x=len(climCXCB.X.sample)-1
n_ess=int((len(climCXCB.law_coef.sample_MCMC)-1)/n_x)

XF_noBE=np.tile(climCXCB.X.loc[:,:,"F","Multi_Synthesis"][:,1:], (1, n_ess)).T
XF_data=np.vstack((XF_noBE,climCXCB.X.loc[:,"BE","F","Multi_Synthesis"]))
XF=xr.DataArray( XF_data.T, coords = [climCXCB.X.time , samples_MCMC  ] , dims = ["time","sample_MCMC"] )

locF  = climCXCB.law_coef.loc["loc0",:,"Multi_Synthesis"] + XF * climCXCB.law_coef.loc["loc1",:,"Multi_Synthesis"]
scaleF = np.exp(climCXCB.law_coef.loc["scale0",:,"Multi_Synthesis"] + XF * climCXCB.law_coef.loc["scale1",:,"Multi_Synthesis"])
shape = climCXCB.law_coef.loc["shape0",:,"Multi_Synthesis"] + xr.zeros_like(locF)
if not os.path.exists(os.path.join( pathOut , ("Parametres_vstime_posterior"+"_"+dt_string+".nc") )):
    print("Calculate Para")
    para=xr.concat([locF, scaleF,shape], pd.Index(["loc","scale","shape"], name='parametre'))
    para.to_netcdf(os.path.join( pathOut , ("Parametres_vstime_posterior"+"_"+dt_string+".nc") ))
para=xr.open_dataarray(os.path.join( pathOut , ("Parametres_vstime_posterior"+"_"+dt_string+".nc")))



#Only prior
n_alltimes=clim.n_time
draw=np.array(np.random.default_rng().multivariate_normal( mean=clim_MM.data.mm_mean,
                                                          cov = clim_MM.data.mm_cov,size=10000))
samples= draw[:,(2*n_alltimes):].T

ind=list(set(Yo.index).intersection(Xo.index))
ns_law.fit(Yo.loc[ind].values,Xo.loc[ind].values)
theta_obs=ns_law.get_params()


custom_params = {"axes.spines.bottom": True,"axes.spines.right": True,"axes.spines.left": True, "axes.spines.top":True}
sns.set_theme(style="whitegrid",rc=custom_params)
    
ofile=pathFig+"/"+dt_string+"_Comp_Distributions_Para_prior_posterior.pdf" 
pdf = mpdf.PdfPages( ofile )

qcoefX=climCXCB.law_coef[:,1:,:].quantile( [ci/2,1-ci/2,0.5] , dim = 'sample_MCMC' ).assign_coords( quantile = ["ql","qu","BE"] )

ymin = float( (climCXCB.law_coef ).min())
ymax = float( (climCXCB.law_coef ).max())
delta = 0.1 * (ymax-ymin)
ylim = (ymin-delta,ymax+delta)

kwargs = {  "showmeans" : False , "showextrema" : False , "showmedians" : True }
models=climCXCB.data.model.values
m='Multi_Synthesis'
fig = plt.figure( figsize = ( 16 , 10 ) )

sub_def=climCXCB.law_coef.coef.values#[2:4]
para_names=climCXCB.ns_law.get_params_names(True)#[2:4]
label = ["Prior","Posterior",'Observation ML fit']
legend = [mplpatch.Patch(facecolor = c , edgecolor = c , label = l , alpha = 0.5 ) for c,l in zip(["blue","red","green"],label)]


for i in range(len(sub_def)):
    
    ax = fig.add_subplot(1,len(sub_def),i+1)
   
    vplot = ax.violinplot( ((climCXCB.law_coef) )[:,1:,:].loc[sub_def[i],:,m]  , **kwargs )
    vplotf = ax.violinplot( samples[i]  , **kwargs )

    
    for pc in vplot["bodies"]:
        pc.set_facecolor("red")
        pc.set_edgecolor("red")
        pc.set_alpha(0.5)
    vplot["cmedians"].set_color("red")  

    for pc in vplotf["bodies"]:
            pc.set_facecolor("blue")
            pc.set_edgecolor("blue")
            pc.set_alpha(0.3)
    vplotf["cmedians"].set_color("blue")
    ax.set_xticks([1])
    
    xticks = [ para_names[i] ]
    ax.set_xticklabels( xticks , fontsize = 13 )
    if ((clim.law_coef) )[:,1:,:].loc[sub_def[i],:,m].min()*((clim.law_coef) )[:,1:,:].loc[sub_def[i],:,m].max()<0:
        ax.hlines( 0 , 1-0.3 , 1+0.3, color = "black" )
    gcm_para=clim.data.law_coef.loc[sub_def[i],'BE' ,:].values[:-1]

    ax.scatter([1]*len(gcm_para),gcm_para,color="blue")
    #Add obs
    ax.hlines( theta_obs[i] , 1-0.3 , 1+0.3, color = "green",linewidth=4 )
    
    if i ==0:
        ax.legend( handles = legend , fontsize = 20,loc='upper left' )


fig.set_tight_layout(True)

#fig.suptitle( "With X sampling", fontsize = 25 )
pdf.savefig(fig)
plt.close(fig)
pdf.close()
































print("Produce RLs and their tex tables")
scenarios = list(sub_scenarios_tables)# list(data.keys())
for T in T_RL:
    
    if not os.path.exists(os.path.join( pathOut , ("RL_F_posterior"+"_"+str(T)+"y_"+dt_string+".nc") )):
        print(" Calculate Return levels: "+ str(T)  )    
        
        RL_F=xr.DataArray(sc.genextreme.ppf((1-1/T),  loc = locF , scale = scaleF , c = - shape).T,
                   coords = [climCXCB.X.time , samples_MCMC  ] ,
                  dims = ["time","sample_MCMC"] )    
        RL_F.to_netcdf(os.path.join( pathOut , ("RL_F_posterior"+"_"+str(T)+"y_"+dt_string+".nc") ))
    RL_F=xr.open_dataarray(os.path.join( pathOut , ("RL_F_posterior"+"_"+str(T)+"y_"+dt_string+".nc") ))
    data={}    
    for scen in scenarios:
        
        time_year_rl=list([1850,2025,2050,2100])
        prod=list(product([scen],time_year_rl))  
        scen_times =   np.array([str(prod[i][1])+"_"+str(prod[i][0]) for i in range(len(prod))])
        quantiles=RL_F.loc[scen_times,:].quantile( [ ci / 2 , 1 - ci / 2 , 0.5 ] , dim = "sample_MCMC" ).assign_coords( quantile = ["ql","qu","BE"] )    
        data[scen]=pd.DataFrame(quantiles +bias_multi,["ql","qu","BE"],time_year_rl)

  

    
    years = list(data[scenarios[0]].columns)
    rows = []
    for year in years:
        row = {}
        for scen in scenarios:
                    be = data[scen].loc['BE', year]
                    ql = data[scen].loc['ql', year]
                    qu = data[scen].loc['qu', year]
                    row[(scen, 'BE')] = round(be, 2)
                    row[(scen, '5-95\\%')] = f"[{ql:.2f}, {qu:.2f}]"
        rows.append(row)

    summary = pd.DataFrame(rows, index=years)
    summary.columns = pd.MultiIndex.from_tuples(summary.columns, names=['Scénario', ''])
    summary.index.name = 'Année'

    latex_code = make_latex_table(
                summary, scenarios, years,
                caption="Instantaneous return levels' best estimates with 2.5-97.5\\% confidence intervals in parentheses, by scenario and year.",
                label="tab:scenarios"
        )
    with open(pathFig+"/"+str(T)+"RL_table_scenarios.tex", "w") as f:
        f.write(latex_code)
    





print("Produce PERs and their tex tables")
T1_l=2025
T2_l=2100
dist=T2_l-T1_l
#Predictive sample
for scen in scenarios_all:
    
    if not os.path.exists(os.path.join( pathOut , (scen+"_pred_tirages_EQR_posterior_"+str(T1_l)+"_"+str(T2_l)+"_"+dt_string+".nc") )):
        para=xr.open_dataarray(os.path.join( pathOut , ("Parametres_vstime_posterior"+"_"+dt_string+".nc")))
        samples_MCMC = para.sample_MCMC
        prod=list(product(list(range(sample_pred)),samples_MCMC.values)) 
        samples_pred =   np.array([str(prod[i][1])+"_"+str(prod[i][0]) for i in range(len(prod))])
        time_year_scen=list(range(T1_l, T2_l+1))
        n_time=len(time_year_scen)

        prod=list(product([scen],time_year_scen))  
        scen_times =   np.array([str(prod[i][1])+"_"+str(prod[i][0]) for i in range(len(prod))])
        shape_s=np.tile(para.loc["shape",:,scen_times ], ( sample_pred,1))
        loc_s=np.tile(para.loc["loc",:,scen_times ], ( sample_pred,1))
        scale_s=np.tile(para.loc["scale",:,scen_times ], ( sample_pred,1))
        tirages_pred=xr.DataArray(sc.genextreme.rvs(c=-shape_s,loc=loc_s,
                            scale=scale_s ).T,
                   coords = [scen_times , samples_pred  ] ,
                  dims = ["time","sample_pred"] ).max('time')
        tirages_pred.to_netcdf(os.path.join( pathOut , (scen+"_pred_tirages_EQR_posterior_"+str(T1_l)+"_"+str(T2_l)+"_"+dt_string+".nc") ))


T2_u=2000
T1_u=T2_u-dist

#Predictive sample hisorical
for scen in scenarios_all:
    
    if not os.path.exists(os.path.join( pathOut , (scen+"_pred_tirages_EQR_posterior_"+str(T1_u)+"_"+str(T2_u)+"_"+dt_string+".nc") )):
        para=xr.open_dataarray(os.path.join( pathOut , ("Parametres_vstime_posterior"+"_"+dt_string+".nc")))
        samples_MCMC = para.sample_MCMC
        prod=list(product(list(range(sample_pred)),samples_MCMC.values)) 
        samples_pred =   np.array([str(prod[i][1])+"_"+str(prod[i][0]) for i in range(len(prod))])
        time_year_hist=list(range(T1_u, T2_u+1))
        n_time=len(time_year_hist)

        prod=list(product([scen],time_year_hist))  
        scen_times =   np.array([str(prod[i][1])+"_"+str(prod[i][0]) for i in range(len(prod))])
        shape_s=np.tile(para.loc["shape",:,scen_times ], ( sample_pred,1))
        loc_s=np.tile(para.loc["loc",:,scen_times ], ( sample_pred,1))
        scale_s=np.tile(para.loc["scale",:,scen_times ], ( sample_pred,1))
        tirages_pred=xr.DataArray(sc.genextreme.rvs(c=-shape_s,loc=loc_s,
                            scale=scale_s ).T,
                   coords = [scen_times , samples_pred  ] ,
                  dims = ["time","sample_pred"] ).max('time')
        tirages_pred.to_netcdf(os.path.join( pathOut , (scen+"_pred_tirages_EQR_posterior_"+str(T1_u)+"_"+str(T2_u)+"_"+dt_string+".nc") ))



period_labels = {'hist': f"{T1_u}-{T2_u}",'future': f"{T1_l}-{T2_l}"}
for per in prob_interest:
    Eq_Reliability_T=1-per
    panda_data_per={}

    for scen in scenarios_all:
        pred_hist=xr.open_dataarray(os.path.join( pathOut , (scen+"_pred_tirages_EQR_posterior_"+str(T1_u)+"_"+str(T2_u)+"_"+dt_string+".nc") ))
        predict_EQR_hist=float(pred_hist.quantile(Eq_Reliability_T))+bias_multi
    
        pred=xr.open_dataarray(os.path.join( pathOut , (scen+"_pred_tirages_EQR_posterior_"+str(T1_l)+"_"+str(T2_l)+"_"+dt_string+".nc") ))
        predict_EQR=float(pred.quantile(Eq_Reliability_T))+bias_multi            
        panda_data_per[scen]=pd.DataFrame([predict_EQR_hist,predict_EQR],["hist", "future"],["Pred"])
        
    data2=panda_data_per
    scenarios = list(data2.keys())
    periods = list(data2[scenarios[0]].index)  # ['hist', 'future']

    rows = []
    for period in periods:
        row = {}
        for scen in scenarios:
            pred = data2[scen].loc[period, 'Pred']
            row[scen] = round(pred, 2)
        rows.append(row)

    summary2 = pd.DataFrame(rows, index=[period_labels[p] for p in periods])
    summary2.index.name = 'Period'
    
    row_labels = [period_labels[p] for p in periods]

    latex_code2 = make_latex_table_styled(
    summary2, scenarios, row_labels, scenario_labels,
    caption="Projected values of the Predictive Equivalent Reliability "+str(int(per*100))+"\\% for historical ("+str(T1_u)+"-"+str(T2_u)+") and future ("+str(T1_l)+"-"+str(T2_l)+") periods.",
    label="tab:pred_scenarios_"+str(int(1/per))+"per"
    )
    with open(pathFig+"/"+str(T1_l)+"_"+str(1/per)+"_per_table_pred_scenarios.tex", "w") as f:
        f.write(latex_code2)
        
        

#Figures
print("Create projection figures")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "lines.linewidth": 1.1,
    "pdf.fonttype": 42,   # keep text editable/searchable in the PDF
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


color_RL = {
    "ssp126": "#003466",
    "ssp245": "#f69320",
    "ssp585": "#980002",
}
color_EQR = color_RL  
HIST_COLOR = "#2E7D32" 
OBS_COLOR = "#1a1a1a"
LIFETIME_COLOR = "#6f9fd8"


for T in T_RL:
    for per in prob_interest:
        p = 1 / T
        ofile = os.path.join(
        pathFig, "PER_" + str(1 / per) + "_RL_" + str(T) + "_ci" + str(ci)
        + "_" + str(T1_l) + "_" + str(T2_l) + "_Multiscenario_Comp_RL_vrs_EQR_ssp585_dfast_NCCstyle.pdf"
        )

        mm = 1 / 25.4
        fig, ax = plt.subplots(figsize=(183 * mm, 100 * mm))

        Eq_Reliability_T = 1 - per

        scen = "ssp245"
        pred = xr.open_dataarray(os.path.join(
            pathOut, scen + "_pred_tirages_EQR_posterior_" + str(T1_u) + "_"
            + str(T2_u) + "_" + dt_string + ".nc"))
        predict_EQR = float(pred.quantile(Eq_Reliability_T))

        ax.hlines(predict_EQR + bias_multi, T1_u, T2_u,
          colors=HIST_COLOR, ls=(0, (4, 2)), lw=1.1, zorder=6)  
        ax.text(T2_l + 2, predict_EQR + bias_multi,
            f"{predict_EQR + bias_multi:.1f}\u00b0C", #round(x,1)
            va="center", ha="left", color=HIST_COLOR, fontsize=6.8)

        for scen in sub_scenarios_figures:
            c = color_RL[scen]
            prod_ = list(product([scen], time_year))
            scen_times = np.array([str(prod_[i][1]) + "_" + str(prod_[i][0])
                                for i in range(len(prod_))])

            RL_F = xr.open_dataarray(os.path.join(
            pathOut, "RL_F_posterior_" + str(T) + "y_" + dt_string + ".nc"))
            Q_1000_quantiles = RL_F.loc[scen_times, :].quantile(
            [ci / 2, 1 - ci / 2, 0.5], dim="sample_MCMC"
            ).assign_coords(quantile=["ql", "qu", "BE"])
            ax.plot(time_period, Q_1000_quantiles.loc["BE", scen_times] + bias_multi,
                 color=c, lw=1.2, zorder=4)
            ax.fill_between(time_period,
                         Q_1000_quantiles.loc["ql", scen_times] + bias_multi,
                         Q_1000_quantiles.loc["qu", scen_times] + bias_multi,
                         color=c, alpha=0.14, lw=0, zorder=2)

            pred = xr.open_dataarray(os.path.join(
            pathOut, scen + "_pred_tirages_EQR_posterior_" + str(T1_l) + "_"
            + str(T2_l) + "_" + dt_string + ".nc"))
            predict_EQR = float(pred.quantile(Eq_Reliability_T))
            ax.hlines(predict_EQR + bias_multi, T1_l, T2_l,
                  colors=color_EQR[scen], ls=(0, (4, 2)), lw=1.0, zorder=6)
            ax.text(T2_l + 2, predict_EQR + bias_multi,
                f"{predict_EQR + bias_multi:.1f}\u00b0C",
                va="center", ha="left", color=color_EQR[scen], fontsize=6.8)

        ax.scatter(Yo.index, Yo + bias_multi, s=7, facecolor=OBS_COLOR,
               edgecolor="none", alpha=0.75, zorder=5)

        ax.axvspan(T1_l, T2_l, color=LIFETIME_COLOR, alpha=0.07, zorder=1, lw=0)
        ax.set_xlim(np.min([1940, T1_u]) - 5, T2_l + 20)
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%d\u00b0C"))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", which="major", color="0.85", lw=0.5, zorder=0)
        ax.tick_params(direction="out", length=3)

        ax.set_xlabel("Year")
        ax.set_ylabel("Global land maximum temperature")  # replace with your actual variable/units


        rl_title = f"Instaneous return levels - "+"{:.3g}%".format(p*100)
        per_title = f"PER - "+"{:.3g}%".format(per*100)

        rl_handles = [Line2D([0], [0], color=color_RL[s], lw=1.2,
                          label=SCENARIO_META[s]["label"]) for s in sub_scenarios_figures]
        per_handles = [Line2D([0], [0], color=HIST_COLOR, lw=1.0, ls=(0, (4, 2)),
                           label=f"Historical ({T1_u}\u2013{T2_u})")] + \
                  [Line2D([0], [0], color=color_EQR[s], lw=1.0, ls=(0, (4, 2)),
                           label=f"{SCENARIO_META[s]['label']} ({T1_l}\u2013{T2_l})")
                   for s in sub_scenarios_figures]
        other_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=OBS_COLOR,
               markersize=4, label="Observations"),
        mplpatch.Patch(facecolor=LIFETIME_COLOR, alpha=0.15,
                        label=f"Project lifetime ({T1_l}\u2013{T2_l})"),
    ]

        leg1 = ax.legend(handles=rl_handles, title=rl_title,
                      loc="upper left", bbox_to_anchor=(1.01, 0.80),#1.01
                      frameon=False, handlelength=1.6, borderaxespad=0,
                      title_fontsize=6.8, labelspacing=0.4, alignment="left")
        leg2 = ax.legend(handles=per_handles, title=per_title,
                      loc="upper left", bbox_to_anchor=(1.01, 0.62),
                      frameon=False, handlelength=1.6, borderaxespad=0,
                      title_fontsize=6.8, labelspacing=0.4, alignment="left")
        ax.legend(handles=other_handles, loc="upper left", bbox_to_anchor=(1.01, 0.37),
              frameon=False, handlelength=1.6, borderaxespad=0, labelspacing=0.4, alignment="left")
        ax.add_artist(leg1)
        ax.add_artist(leg2)

        fig.subplots_adjust(right=0.72)

        fig.savefig(ofile, bbox_inches="tight")
        plt.close(fig)

print("Create parameter distribution figure for "+str(seuil)+"°C")
        




pred=xr.open_dataarray(os.path.join( pathOut , (scen_seuil+"_pred_tirages_EQR_posterior_"+str(T1_l)+"_"+str(T2_l)+"_"+dt_string+".nc") ))
data_theta=pred.to_pandas().T+bias_multi
samples_select=np.unique(data_theta[data_theta.values>seuil].index)


samples_kept_para=[str(samples_select[i]).split("_")[0]+"_"+str(samples_select[i]).split("_")[1] for i in range(len(samples_select))]
samples_para=np.unique(samples_kept_para, return_counts=True)[0]
samples_kept_para_weight=np.unique(samples_kept_para, return_counts=True)[1]

samples_kept_para_X=[str(samples_select[i]).split("_")[0] for i in range(len(samples_select))]
samples_para_X=np.unique(samples_kept_para_X, return_counts=True)[0]
samples_kept_para_X_weight=np.unique(samples_kept_para_X, return_counts=True)[1]


X_all=climCXCB.X.loc["2100_ssp585",:,'F',:].values.flatten()
X_select=climCXCB.X.loc["2100_ssp585",samples_para_X,'F',:].values.flatten()
clim_CXCB_select=climCXCB.law_coef.loc[:,samples_para,:]



# ── Typography & RC ───────────────────────────────────────────────────────────
mpl.rcParams.update({
    "font.family":        "Arial",          # Nature house font
    "font.size":          7,
    "axes.titlesize":     7,
    "axes.labelsize":     7,
    "xtick.labelsize":    6,
    "ytick.labelsize":    6,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "xtick.major.size":   2.5,
    "ytick.major.size":   2.5,
    "xtick.major.width":  0.6,
    "ytick.major.width":  0.6,
    "axes.linewidth":     0.6,
    "lines.linewidth":    1.0,
    "legend.fontsize":    6,
    "legend.frameon":     False,
    "pdf.fonttype":       42,               # embed fonts as Type 1 (Nature req.)
    "ps.fonttype":        42,
})

# ── Palette ───────────────────────────────────────────────────────────────────
RED   = "#C0392B"    # all posterior
BLUE  = "#2471A3"    # weighted subset
BLACK = "#1A1A1A"


# ── Layout constants ──────────────────────────────────────────────────────────
N_PANELS   = climCXCB.n_coef + 1          # parameter panels + X panel
FIG_W_IN   = 3.54                         # 90 mm — Nature single column
FIG_H_IN   = N_PANELS * 1.15
LABEL_FONT = dict(fontsize=7, fontweight="bold", va="top", ha="left")

fig, axes = plt.subplots(
    N_PANELS, 1,
    figsize=(FIG_W_IN, FIG_H_IN),
    constrained_layout=True,
)

# ── Parameter panels ──────────────────────────────────────────────────────────
panel_labels = list(string.ascii_lowercase)
for i in range(climCXCB.n_coef):
    ax_l = axes[i]
    ax_r = ax_l.twinx()

    param_name = climCXCB.ns_law.get_params_names(True)[i]

    # Data
    vals_all = climCXCB.law_coef[:, :, :].loc[
        climCXCB.law_coef.coef[i], :, "Multi_Synthesis"].values.ravel().astype(float)
    vals_sel = clim_CXCB_select[:, :, :].loc[
        climCXCB.law_coef.coef[i], :, "Multi_Synthesis"].values.ravel().astype(float)
    wts      = (samples_kept_para_weight / 10).ravel().astype(float)

    # Histograms + KDE overlays
    _hist_kde(ax_l, vals_all, color=RED,  lw=1.1, bins=80,
              label="All posterior samples")
    _hist_kde(ax_r, vals_sel, weights=wts, color=BLUE, lw=1.1, bins=80,
              label=f"Samples over {seuil} °C")

    # Medians
    med_all = float(np.median(vals_all))
    med_sel = float(weighted_median(vals_sel, wts))
    _vline_with_label(ax_l, med_all, RED,  side="left")
    _vline_with_label(ax_l, med_sel, BLUE, side="right")

    # GCM ticks 
    gcm = np.asarray(
        clim.data.law_coef.loc[
            climCXCB.law_coef.coef[i], "BE", :].values[:-1], dtype=float)
    _gcm_dots(ax_r, gcm)

    _style_axes(ax_l, ax_r, xlabel=param_name)

    # Panel label
    ax_l.text(0.02, 1.05, panel_labels[i],
              transform=ax_l.transAxes, **LABEL_FONT)

    # Legend — first panel only
    if i == 0:
        h1, lb1 = ax_l.get_legend_handles_labels()
        h2, lb2 = ax_r.get_legend_handles_labels()
        # Add GCM proxy
        gcm_proxy = mpl.lines.Line2D(
            [], [], marker="|", color=BLACK, lw=0, ms=5, mew=0.8,
            label="Individual GCM fits")
        ax_l.legend(
            h1 + h2 + [gcm_proxy],
            lb1 + lb2 + ["Individual GCM fits"],
            loc="upper left", handlelength=1.2,
            borderpad=0.45, labelspacing=0.3,
        )

# ── X panel (Global Tmean) ────────────────────────────────────────────────────
ax_l = axes[-1]
ax_r = ax_l.twinx()

vals_X_all = np.asarray(X_all, dtype=float).ravel()
vals_X_sel = np.asarray(X_select, dtype=float).ravel()
wts_X      = (samples_kept_para_X_weight / (10 * 100)).ravel().astype(float)

_hist_kde(ax_l, vals_X_all, color=RED,  lw=1.1, bins=60,
          label="All posterior samples")
_hist_kde(ax_r, vals_X_sel, weights=wts_X, color=BLUE, lw=1.1, bins=50,
          label=f"Samples over {seuil} °C")

med_X_all = float(np.median(vals_X_all))
med_X_sel = float(weighted_median(vals_X_sel, wts_X))
_vline_with_label(ax_l, med_X_all, RED,  side="left")
_vline_with_label(ax_l, med_X_sel, BLUE, side="right")
gcm_X = np.asarray(
    clim.data.X.loc["2100_ssp585", "BE", "F", :].values[:-1], dtype=float)
_gcm_dots(ax_r, gcm_X)
_style_axes(ax_l, ax_r,
            xlabel="Global mean temperature over land, SSP5-8.5, 2100 (°C)")
ax_l.text(0.02, 0.97, panel_labels[climCXCB.n_coef],
          transform=ax_l.transAxes, **LABEL_FONT)


ofile = os.path.join(
    pathFig,
    f"{T1_l}_{T2_l}_{dt_string}_Weighted_UP_pred{seuil}"
    f"_{dt_string}_Comp_Distributions_Para_X_Publication.pdf"
)
fig.savefig(ofile, dpi=300, bbox_inches="tight")
plt.close(fig)



##QQ Plot
print("Create QQ plot")


qparamsX =para[:,1:,:].quantile( [ ci / 2 , 1 - ci / 2 , 0.5 ] , dim = "sample_MCMC" ).assign_coords( quantile = ["ql","qu","BE"] )

fig = plt.figure( figsize = (15,5) )
j=1
for scen in scenarios:
	#KS test
	print(scen)
	prod=list(product([scen],Yo.index))  
	scen_times =   np.array([str(prod[i][1])+"_"+str(prod[i][0]) for i in range(len(prod))])

	loc   = qparamsX.loc["BE",'loc',scen_times].values
	scale = qparamsX.loc["BE",'scale',scen_times].values
	shape = qparamsX.loc["BE",'shape',scen_times].values
	Z = ( Yo.values.squeeze() - loc ) / scale
	print(sc.kstest( Z , lambda x : sc.genextreme.cdf( x , loc = 0 , scale = 1 , c = - shape )) )#p-value > 0.05, no rejection.


	ax = fig.add_subplot( 1 , len(scenarios), j )
	residuals=[]
	for i in range(0,len(Yo)):
    		shape=qparamsX.loc["BE",'shape',str(Yo.index[i])+"_"+scen].values
    		loc=qparamsX.loc["BE",'loc',str(Yo.index[i])+"_"+scen].values
    		scale=qparamsX.loc["BE",'scale',str(Yo.index[i])+"_"+scen].values
        
    		residuals.append(((Yo.iloc[i].values[0]-loc)/scale))

	sm.qqplot(np.asarray(residuals),dist=sc.genextreme,distargs=(-shape,), line="45",ax=ax)
	ax.set_ylabel("Sample Quantile for Observations" )
	#ax.set_xlabel("Théorique Posterior (GEV)" )
	ax.set_title(scen)
	j=j+1

plt.savefig(os.path.join( pathFig ,X_ref+"_"+Y_ref+"_QQplots.pdf"))
plt.show()

