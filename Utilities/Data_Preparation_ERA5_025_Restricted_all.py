
#####Packages

import sys,os
#print(sys.version)
import shutil
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
from matplotlib.lines import Line2D
import warnings
warnings.simplefilter("ignore")

#Arguments

Y_ref=sys.argv[1]
X_ref=sys.argv[2]

work_path=sys.argv[3]

path_X_Modeles_Brut=sys.argv[4]
path_Y_Modeles_Brut=sys.argv[5]

path_X_Obs_Brut=sys.argv[6]
path_Y_Obs_Brut=sys.argv[7]

station_type=sys.argv[8]


#Paths 

print("Data Preparation for: "+Y_ref)
basepath=os.path.join(work_path,"data")
print("Working path: "+basepath)
if os.path.exists(basepath):
    shutil.rmtree(basepath)

pathFig=os.path.join(basepath,"Figures/"+Y_ref)
if not os.path.exists(pathFig):
    os.makedirs(pathFig)
    
pathTemp=os.path.join(basepath,"temp")
if not os.path.exists(pathTemp):
    os.makedirs(pathTemp)
    
pathOut_X=os.path.join(basepath,"CMIP6/03_Post_Treatment/TasMean_"+X_ref)
if not os.path.exists(pathOut_X):
    os.makedirs(pathOut_X)
X_files="{}_tas_YearMean_"+X_ref+".nc"

pathOut_Y=os.path.join(basepath,"CMIP6/03_Post_Treatment/TXX_"+Y_ref)
if not os.path.exists(pathOut_Y):
    os.makedirs(pathOut_Y)
if station_type=='True':
	Y_files="{}_tasmax_YearMax_"+Y_ref+".nc"
else:
	Y_files="{}_tasmax_YearMax_TXX_"+Y_ref+".nc"
pathOut_Obs=os.path.join(basepath,"Observations/03_Post_Treatment/"+Y_ref+"_"+X_ref)
if not os.path.exists(pathOut_Obs):
    os.makedirs(pathOut_Obs)


    
#Functions
def correct_miss( X , lo =  100 , up = 365 ):
	mod = str(X.columns[0])
	bad = np.logical_or( X < lo , X > up )
	bad = np.logical_or( bad , np.isnan(X) )
	bad = np.logical_or( bad , np.logical_not(np.isfinite(X)) )
	if np.any(bad):
		print("Bad raw values "+bad.columns[0] +": "+str(bad.sum().values[0]))
		idx,_ = np.where(bad)
		idx_co = np.copy(idx)
		for i in range(idx.size):
			j = 0
			while idx[i] + j in idx:
				j += 1
			idx_co[i] += j
			if (idx_co[i]>=len(X)):
				idx_co[i] =idx_co[i-1] 
		X.iloc[idx] = X.iloc[idx_co].values
	return X

## Corrections Outliers
def is_outlier(Y,log_outliers):
    Y_means=[]
    Y=Y[Y.columns[0]]
    Y.index=pd.to_datetime(Y.index.astype(str), yearfirst=True)

    if sum(Y.index!=Y.sort_index().index)>0:
        print("index issues")
        print(Y.index[(Y.index!=Y.sort_index().index)])
        log_outliers.write(str(Y.index[(Y.index!=Y.sort_index().index)])+"\n")
        Y=Y.sort_index()
        Anomalies=(abs(Y)>4)*0
    else:

        #freq_roll="3650D" #10years
        #freq_roll="7300D" #20 ans
        #freq_roll="1825D" #5 ans
        freq_roll="10950D" #30 ans
        mean_rolling_Y=Y.rolling(freq_roll,min_periods=1, center=True).mean()
        std_rolling_Y=Y.rolling(freq_roll,min_periods=1, center=True).std()
        #std_rolling_Y[std_rolling_Y<1.6]=1.66666
        #Supposant une loi normale car n grand (min 30, voir 25*30)
        #Probabilité 3eq= 99.7%
        #Probabilité 4eq= 99.99%
        #pas de prise en compte de la queue
        Anomalies=(abs(Y-mean_rolling_Y)/std_rolling_Y>4)*1
    Anomalies.index=Anomalies.index.year.astype(int)
    return Anomalies
    
    
#General args
time_period    = np.arange( 1850 , 2101 , 1 , dtype = int )
time_reference = np.arange( 1986 , 2016 , 1 , dtype = int )

exclusion=['NorESM2-LM_i1p1f1','UKESM1-0-LL_i1p1f1','UKESM1-0-LL_i1p1f2','UKESM1-0-LL_i1p1f3','HadGEM3-GC31-LL_i1p1f3','HadGEM3-GC31-MM_i1p1f3','KACE-1-0-G_i1p1f1']
########### Treatment Obs
print("Start treatment of Observation Data")
#For X0, covariate
dXo = xr.open_dataset(path_X_Obs_Brut,decode_times=False) #Deja en anomalies
Xo  = pd.DataFrame( dXo.tas.values.squeeze() , columns = ["Xo"] , index = dXo.time.values )
#Obs Anomalies
Xo_ano =Xo-Xo.loc[time_reference].mean()
#Save obs
Xo_ano.to_xarray().to_netcdf( os.path.join( pathOut_Obs,"Xo_Ano.nc" ) )
print("Example for Tasmean: "+str(Xo_ano.iloc[-1]))

#For Y0, interest var
if station_type=='True':	
    bashCommand = "cdo -selvar,TX -yearmax "+path_Y_Obs_Brut+"_*.nc"+" "+pathTemp+"/Yo_test.nc"
    os.system(bashCommand)
    dYo = xr.open_dataset(os.path.join(pathTemp,"Yo_test.nc"))
    Yo  = pd.DataFrame( dYo.TX.values.squeeze() , columns = ["Yo"] , index = 		dYo.time["time.year"].values )	
    Yo=Yo[np.isfinite(Yo.values)]
elif path_Y_Obs_Brut[-3:]=="csv":
    dYo =  pd.read_csv(path_Y_Obs_Brut)
    Yo  = pd.DataFrame( dYo['max'].values.squeeze() , columns = ["Yo"] , index = dYo['extract'].values )
    #LImitation à 70 ans 
    Yo=Yo[:70]
    Yo=Yo[np.isfinite(Yo.values)]
elif path_Y_Obs_Brut[-2:]=="nc":
    dYo = xr.open_dataset(path_Y_Obs_Brut,decode_times=False)
    Yo  = pd.DataFrame( dYo.tas.values.squeeze() , columns = ["Yo"] , index = dYo.time.values+1900 )

# kelvin to celsius
if Yo.mean().values>273.15:
	Yo=Yo-273.15
##Obs Anomalies
bias = { "Multi_Synthesis" : Yo.loc[time_reference].mean().values }
Yo_ano =Yo- bias["Multi_Synthesis"]

#Save obs and bias Observations
Yo_ano.to_xarray().to_netcdf( os.path.join( pathOut_Obs,"Yo_Ano.nc" ) )
pd.DataFrame(bias).to_xarray().to_netcdf( os.path.join( pathOut_Obs,"bias_Obs.nc" ) )
print("Example for TXxx: "+str(Yo_ano.iloc[-1]))
print("Observation bias: "+str(bias))





########### Treatment CMIP
print("Start treatment of CMIP Data")
#Conserve si meme Modele+type+scenario pour tas et tasmax
## List of models X
modelsX = [  "_".join(f.split("/")[-1][:-3].split("_")[:3]) for f in os.listdir(path_X_Modeles_Brut) ] #Modele+type+scenario
modelsX.sort()

## List of models Y
modelsY =  [  "_".join(f.split("/")[-1][:-3].split("_")[:3]) for f in os.listdir(path_Y_Modeles_Brut) ] #Modele+type+scenario
modelsY.sort()

models = list(set(modelsX) & set(modelsY)) #Conserve si meme Modele+type+scenario pour tas et tasmax
models.sort()


#Only X treatment absolute
lX = []
for m in models:
        model, typeRun, scen=m.split("_")

        ## Load X
        
        df   = xr.open_dataset( os.path.join( path_X_Modeles_Brut , X_files.format(m) )  ,decode_times=False )
        time = df.time.values.astype(int)
        X    = pd.DataFrame( df.tas.values.ravel() , columns = [m] , index = time )
        X_miss=correct_miss(X)
        X=pd.DataFrame( np.array([X_miss.values.ravel(),
          np.array([model for i in range(len(time))]),
          np.array([typeRun for i in range(len(time))]),
          np.array([scen for i in range(len(time))])]).T, 
          columns = ['tas',"Model", "typeRun", "scenario"] , index = time )
        X['tas']=X['tas'].astype(float)
        lX.append(  X)
lX_out = pd.concat(lX)
lX_out.index.name= 'Year'

scenarios=np.unique(lX_out['scenario'])[1:] #liste scenarios
models_reduced=np.unique([m.split("_")[0]+"_"+m.split("_")[1] for m in models]) #Liste models/type 
#Only Y treatment absolute
lY = []
for m in models:
        model, typeRun, scen=m.split("_")
        ## Load Y
        df   = xr.open_dataset( os.path.join( path_Y_Modeles_Brut , Y_files.format(m) ) ,decode_times=False  )
        time = df.time.values.astype(int)
        Y    = pd.DataFrame( df.tasmax.values.ravel() , columns = [m] , index = time )
        Y_miss=correct_miss(Y)
        Y=pd.DataFrame( np.array([Y_miss.values.ravel(),
          np.array([model for i in range(len(time))]),
          np.array([typeRun for i in range(len(time))]),
          np.array([scen for i in range(len(time))])]).T, 
          columns = ['tasmax',"Model", "typeRun", "scenario"] , index = time )
        Y['tasmax']=Y['tasmax'].astype(float)
        lY.append( Y)
lY_out = pd.concat(lY)
lY_out.index.name= 'Year'

scenarios=np.unique(lY_out['scenario'])[1:] #liste scenarios

print("CMIP Scenarios available: "+str(scenarios))
models_reduced=np.unique([m.split("_")[0]+"_"+m.split("_")[1] for m in models]) #Liste models/type 

#Only modele with scenario (If only historical, no year past 2014 -> exclude)
#Only modele with historical (If only scenario, no year before 2014 -> exclude)

only_2015_Y=lY_out.loc[lY_out.index==2015  ]
only_2003_Y=lY_out.loc[lY_out.index==2003  ]
only_kept=pd.merge(only_2015_Y, only_2003_Y, how='inner', on=['Model', 'typeRun'])
models_kept=np.unique(only_kept['Model'])
typeRun_kept=np.unique(only_kept['typeRun'])
lY_out=lY_out.loc[(lY_out['Model'].isin(models_kept))&
                 (lY_out['typeRun'].isin(typeRun_kept)) ]

lX_out=lX_out.loc[(lX_out['Model'].isin(models_kept))&
                 (lX_out['typeRun'].isin(typeRun_kept)) ]




models_reduced_sce=np.unique([only_kept['Model'].iloc[i]+"_"+only_kept['typeRun'].iloc[i] for i in range(len(only_kept))]) #Liste models/type 
models_reduced_full=models_reduced
models_reduced=models_reduced_sce

main_list = list(set(models_reduced) - set(exclusion))
main_list.sort()
models_reduced=main_list

print("Number of CMIP Models available: "+str(len(models_reduced)))
print(models_reduced)
print("Models manually excluded :")
print(exclusion)


#### 




#Transformer en Anomalies for X
print("Start Anomalie transformation")
for m1 in models_reduced:
        X=lX_out.loc[(lX_out['Model']==m1.split("_")[0])&
                 (lX_out['typeRun']==m1.split("_")[1]) ]
        b_X=X['tas'].loc[time_reference].mean()
        lX_out.loc[(lX_out['Model']==m1.split("_")[0])&
                 (lX_out['typeRun']==m1.split("_")[1]) ,'tas']-=b_X
                 
#Transformer en Anomalies for Y
bias_GCM = {}
for m1 in models_reduced:        
        Y=lY_out.loc[(lY_out['Model']==m1.split("_")[0])&
                 (lY_out['typeRun']==m1.split("_")[1]) ]
        b_Y=Y['tasmax'].loc[time_reference].mean()
        bias_GCM [m1.split("_")[0]] = b_Y- 273.15 #biais en celsius
        lY_out.loc[(lY_out['Model']==m1.split("_")[0])&
                 (lY_out['typeRun']==m1.split("_")[1]) ,'tasmax']-=b_Y




    

#Delete sup years for X
scen=np.unique(lX_out['scenario'])
#Ssp 245 starting in 2014 for EC-Earth3-Veg
lX_out=lX_out.loc[-((lX_out['scenario'].isin( scen[1:])) &(lX_out.index==2014))]
# Historical up to 2016 for FGOALS-g3
lX_out=lX_out.loc[-((lX_out['scenario'].isin( [scen[0]])) &(lX_out.index>2014)) ]

#Delete sup years for Y
scen=np.unique(lY_out['scenario'])
#Ssp 245 starting in 2014 for EC-Earth3-Veg
lY_out=lY_out.loc[-((lY_out['scenario'].isin( scen[1:])) &(lY_out.index==2014))]
# Historical up to 2016 for FGOALS-g3
lY_out=lY_out.loc[-((lY_out['scenario'].isin( [scen[0]])) &(lY_out.index>2014)) ]



#### Outlier treatment
print("Start Outlier treatment")
scen=np.unique(lX_out['scenario'])
log_outliers=open(pathFig+'/X_Outliers_log.txt','a')
lX = []
#        lX.append(  X)
for m1 in models_reduced:
    
    for s in scen:
        X=lX_out.loc[(lX_out['Model']==m1.split("_")[0])&
                 (lX_out['typeRun']==m1.split("_")[1])&(lX_out['scenario']==s) ] 
        log_outliers.write(str(m1.split("_")[0])+"\n")
        

        Outliers=is_outlier(X,log_outliers)
        if sum(Outliers)>0:
            print("For Tasmean: "+m1.split("_")[0]+" " +s+" :"+str(sum(Outliers))+" Outliers")
            log_outliers.write(str(X[Outliers!=0].loc[:,['tas','scenario']])+"\n")
            log_outliers.write("\n")
            X=X[Outliers==0]
        lX.append(  X)
log_outliers.close()

lX_out_ano_clean = pd.concat(lX)
lX_out_ano_clean.index.name= 'Year'

scen=np.unique(lY_out['scenario'])
log_outliers=open(pathFig+'/Y_Outliers_log.txt','a')
lY = []
#        lX.append(  X)
for m1 in models_reduced:
    for s in scen:
        Y=lY_out.loc[(lY_out['Model']==m1.split("_")[0])&
                 (lY_out['typeRun']==m1.split("_")[1])&(lY_out['scenario']==s) ] 
        log_outliers.write(str(m1.split("_")[0])+"\n")
        
        Outliers=is_outlier(Y,log_outliers)
        if sum(Outliers)>0:
        #log
        

            print("For Txxx: "+m1.split("_")[0]+" " +s+" :"+str(sum(Outliers))+" Outliers")
            log_outliers.write(str(Y[Outliers!=0].loc[:,['tasmax','scenario']])+"\n")
            log_outliers.write("\n")
            Y=Y[Outliers==0]
        lY.append(  Y)
log_outliers.close()

lY_out_ano_clean = pd.concat(lY)
lY_out_ano_clean.index.name= 'Year'

    
    
    
######Save data
directory=os.path.join( pathOut_X,"separated")
if not os.path.exists(directory):
    os.makedirs(directory)
print("Start saving Tasmean data")
for m1 in models_reduced:
    scen_mod=np.unique(lX_out_ano_clean.loc[(lX_out_ano_clean['Model']==m1.split("_")[0])&
                 (lX_out_ano_clean['typeRun']==m1.split("_")[1]),'scenario'] ) #Only scenarios relevant for this modele
    for s in scen_mod:
        nameX="{}_{}_tas_YearMean_Europe.nc".format(m1,s) 
        X=lX_out_ano_clean.loc[(lX_out_ano_clean['Model']==m1.split("_")[0])&
                 (lX_out_ano_clean['typeRun']==m1.split("_")[1])&(lX_out_ano_clean['scenario']==s),'tas' ] 
        X.index.name='time'
        X.to_xarray().to_netcdf( os.path.join( directory,nameX ) )

directory=os.path.join( pathOut_Y,"separated" )
if not os.path.exists(directory):
    os.makedirs(directory)
    
print("Start saving TXxx data")
for m1 in models_reduced:
    scen_mod=np.unique(lY_out_ano_clean.loc[(lY_out_ano_clean['Model']==m1.split("_")[0])&
                 (lY_out_ano_clean['typeRun']==m1.split("_")[1]),'scenario'] ) #Only scenarios relevant for this modele
    for s in scen_mod:

    
        nameY="{}_{}_tasmax_YearMax.nc".format(m1,s) 
        Y=lY_out_ano_clean.loc[(lY_out_ano_clean['Model']==m1.split("_")[0])&
                 (lY_out_ano_clean['typeRun']==m1.split("_")[1])&(lY_out_ano_clean['scenario']==s) ,'tasmax'] 
        Y.index.name='time'
        Y.to_xarray().to_netcdf( os.path.join( directory,nameY ) )
        
        
#### Figures
print("Start Figures")

print("Create figure with all the GCM data")

sub_scenarios=["ssp126","ssp245","ssp585"]


SCENARIO_LABEL = {
    'ssp119': 'SSP1-1.9',
    'ssp126': 'SSP1-2.6',
    'ssp245': 'SSP2-4.5',
    'ssp370': 'SSP3-7.0',
    'ssp585': 'SSP5-8.5',
}
MODEL_COLOR = "0.55" 

fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey="row")

for col, sc in enumerate(sub_scenarios):
    Y_means=lY_out_ano_clean.loc[(lY_out_ano_clean['scenario'].isin(['historical', sc]))].groupby('Year').aggregate({lY_out_ano_clean.columns[0]:"mean"})


    #TXxx
    ax = axes[0, col]
    for m1 in models_reduced:
        Y = lY_out_ano_clean.loc[
            (lY_out_ano_clean['scenario'].isin(['historical', sc])) &
            (lY_out_ano_clean['Model'] == m1.split("_")[0]) &
            (lY_out_ano_clean['typeRun'] == m1.split("_")[1])
        ]
        ax.plot(Y.index, Y['tasmax'].values, linestyle="", marker=".", alpha=0.3)
    ax.plot(Yo_ano.index, Yo_ano.values, color="black", linestyle="",
             marker="x", markersize=6, markeredgewidth=1.5, zorder=5,
             label="Observations (ERA5)")
    ax.plot(Y_means.index,Y_means.values, label="Multi-model mean", color = "black"  )
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f°C'))
    ax.set_title(SCENARIO_LABEL[sc])

    #tasmean
    ax = axes[1, col]
    for m1 in models_reduced:
        X = lX_out_ano_clean.loc[
            (lX_out_ano_clean['scenario'].isin(['historical', sc])) &
            (lX_out_ano_clean['Model'] == m1.split("_")[0]) &
            (lX_out_ano_clean['typeRun'] == m1.split("_")[1])
        ]
        ax.plot(X.index, X['tas'].values, linestyle="", marker=".", alpha=0.3)
    ax.plot(Xo_ano.index, Xo_ano.values, color="black", linestyle="",
             marker="x", markersize=6, markeredgewidth=1.5, zorder=5,
             label="Observations (ERA5)")

    X_means=lX_out_ano_clean.loc[(lX_out_ano_clean['scenario'].isin(['historical', sc]))].groupby('Year').aggregate({X.columns[0]:"mean"})
    ax.plot(X_means.index,X_means.values, label="Multi-model mean", color = "black"  )
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f°C'))


axes[0, 0].set_ylabel("TXxx anomaly")
axes[1, 0].set_ylabel("LSAT anomaly")

fig.suptitle("Historical and scenario anomalies for 1986-2015", y=1.02)


model_proxy = Line2D([0], [0], linestyle="", marker=".", color=MODEL_COLOR,
                      alpha=0.6, markersize=8, label="Model realizations (annual)")
obs_handle, obs_label = axes[0, 0].get_legend_handles_labels()
fig.legend(handles=[model_proxy] + obs_handle, loc="upper right",
           bbox_to_anchor=(1.0, 1.05))


fig.tight_layout()
fig.savefig(os.path.join(pathFig, "Paper_Anomaly_XY_All_AllValues_grid.png"),
            bbox_inches="tight")
plt.close(fig)


print("Create table with number of run per models")
def count_runs_table(df, value_col):
    """
    Compte le nombre de runs disponibles par (Model, scenario),
    en se basant sur le nombre de lignes par (Model, scenario, Year).
    On prend le maximum observé sur les années comme estimation
    du nombre de runs (le plus robuste si certaines années
    ont des runs manquants).
    """
    counts_per_year = (
        df.groupby(["Model", "scenario", "Year"])
          .size()
          .reset_index(name="n_runs")
    )

    # Nombre "représentatif" de runs = max observé sur les années
    n_runs = (
        counts_per_year.groupby(["Model", "scenario"])["n_runs"]
        .max()
        .reset_index()
    )

    # Pivot : Model en lignes, scenario en colonnes
    table = n_runs.pivot(index="Model", columns="scenario", values="n_runs")

    # Optionnel : vérifier si le nombre de runs est bien constant sur les années
    consistency = (
        counts_per_year.groupby(["Model", "scenario"])["n_runs"]
        .nunique()
        .reset_index(name="n_distinct_counts")
    )
    inconsistent = consistency[consistency["n_distinct_counts"] > 1]
    if not inconsistent.empty:
        print(f"[{value_col}] Attention : nombre de runs non constant sur les années pour :")
        print(inconsistent)

    return table


# Tables pour chaque variable
table_tasmax = count_runs_table(lY_out_ano_clean, "tasmax")
table_tas    = count_runs_table(lX_out_ano_clean, "tas")

# Fusion en une seule table avec un niveau de colonnes par variable
combined = pd.concat(
    {"tasmax": table_tasmax, "tas": table_tas},
    axis=1
)
combined = combined.fillna(0).astype(int)  # ou "Int64", mais int suffit sans NaN

n_index = 1  # une colonne d'index (Model)
n_tasmax = table_tasmax.shape[1]  # nb de scénarios pour tasmax
n_tas = table_tas.shape[1]        # nb de scénarios pour tas

# Format : l pour l'index, puis colonnes tasmax, une barre verticale, puis colonnes tas
column_format = "l" + "c" * n_tasmax + "|" + "c" * n_tas

latex_str = combined.to_latex(
    caption="Number of a available runs per model and scenarios, for tasmax and tas.",
    label="tab:n_runs_per_model",
    multicolumn=True,
    multicolumn_format="c",
    column_format=column_format,
)
with open(pathFig+"/n_runs_table.tex", "w") as f:
    f.write(latex_str)    
