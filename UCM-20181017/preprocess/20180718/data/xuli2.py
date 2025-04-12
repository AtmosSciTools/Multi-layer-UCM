####
####
import pandas as pd
import sys, glob, os



idir = "Kugahara-full-datasets"
odir = 'ucm_data/'
if not os.path.exists(odir): os.makedirs(odir)
ifiles = glob.glob(idir+'/ts*')
dp = pd.read_csv('Kugahara-full-datasets/tmp_profile_20050901.csv', index_col=0, parse_dates=True)
dp.columns = [float(c.split('m')[0][4:]) for c in dp.columns]
              
for f in ifiles[2:]:
    print f
    d = pd.read_csv(f,index_col=0,parse_dates=True)

    ds = d[[u'Sdown (W m-2)', u'Ldown (W m-2)', u'U (m s-1)', u'Tsonic (oC)', u'q (g kg-1)']]
    ds.columns = ['swtop', 'lwtop', 'wstop', 'ttop', 'qtop']
    if f[-12:-4] == "20050901": ds.loc[:,'ttop'] = dp[28].values 
    ds.index = range(24)
    do = pd.concat([ds,ds,ds,ds,ds[0:1]],ignore_index=True)
    ofile = odir + 'kugahara_topvars' + f[-13:]
    print ofile
    do.to_csv(ofile, index=None)
    new_colors = ['#1f77b4', '#d62728','#2ca02c',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
              
    h = d[[u'Sdown (W m-2)', u'Sup (W m-2)', u'Ldown (W m-2)', u'Lup (W m-2)',u'H (W m-2)',u'LE (W m-2)' ]]
    h.columns = ["Sdown", "Sup", "Ldown", "Lup", "H", "LE"]
    h["Rnet"] = h.Sdown-h.Sup+h.Ldown-h.Lup
    fig = plt.figure(figsize=(5,3.3))
    ax = fig.add_axes([0.15,0.15,0.8,0.75])
    h.index = range(24)
    mk = ["x", "o", "", ""]
    for i, c in enumerate(h[["Rnet", "H", "LE"]].columns):
        h[c].plot(ax=ax, c = new_colors[i],marker=mk[i],ms=4,lw=1.2)
    ax.legend(["Rnet", "H", "Le"])
    ax.set_ylim(-110,750)
    ax.set_xlim(-0.5,23.5)
    ax.set_ylabel('Heat flux ($W/m^{2}$)')
    ax.set_xlabel("LST (hour)")
    ax.axhline(linewidth=1, color='k',ls="--")
    plt.savefig("heat_budget_top.png", dpi=300)







            





