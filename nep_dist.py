import dfmux_calc as d
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import nep_calc as nep



import matplotlib
matplotlib.use('Agg')


import sys

path = sys.argv[1]
sys.path.append(path)
from config import *


#band to calculate 
opts = [0.29184686,0.24187821,0.26855768,0.30598544,0.27094372,0.2958127,0.3278966,0.31625242,0.37653296,0.32719336,
                         0.30607807,0.35566829,0.35588565,0.43860747,0.42063367,0.39077456,0.35717258,0.62886291,0.4708167,0.37702823,
                         0.2996769,0.22026075]

seed = np.random.rand(n_sq_draw)
# changing some SQUID parameters to ones with spread, rescaling to appropriate array size
saa.zt = norm.ppf(seed,loc=saa.zt*sratio,scale=zt_spread)
saa.rdyn = norm.ppf(seed,loc=saa.rdyn*sratio*saa.n_parallel,scale=zdyn_spread)
saa.inoise  = norm.ppf(np.abs(1-seed),loc=saa.inoise / np.sqrt(nsqratio),scale=nei_spread)

#more parameters with no spread
saa.lin = saa.lin * nsqratio



#bolometer parameters and spreads
#draw one version of a comb of 68 bolos that gets used for all SQUIDs, CSF only gets calculated once
n=len(bias_f)
bolo.r  = np.random.normal(dist_rbolo_mean, dist_rbolo_spread, n)
bolo.rstray  = np.abs(np.random.normal(dist_rstray_mean, dist_rstray_spread, n))




dfm = d.dfmux_noise(saa,bolo,wh,para,nuller_cold=nuller_cold)
litebirb = nep.experiment('litebird',dfm)
dfm.init_freq(bias_f,skip_spice=False)

plt.plot(bias_f/1e6,dfm.bolo.r , label='$R_{bolo}$')
plt.plot(bias_f/1e6,dfm.bolo.r + dfm.bolo.rstray , label='$R_{bolo} + R_{stray}$')
plt.xlabel('Bias Frequency [MHz]')
plt.ylabel('Bolometer resistance [$\Omega$]')
plt.legend()
plt.savefig(path+'/nsq_'+str(nsq)+'_bolo_r.png')
plt.close()


plt.plot(bias_f/1e6, dfm.csf[0], label='CS')
plt.plot(bias_f/1e6, 1/dfm.tf[0],label='1/TF',c='tab:orange',alpha=0.2)
for i in range(1,100):
    plt.plot(bias_f/1e6, 1/dfm.tf[i],c='tab:orange',alpha=0.2)
plt.xlabel('Bias Frequency [MHz]')
plt.legend()
#plt.ylabel('Current sharing factor')
plt.title('Current sharing factor and SQUID transfer function')
plt.savefig(path+'/nsq_'+str(nsq)+'_csf.png')
plt.close()


for i in range(100):
    plt.plot(bias_f,dfm.total[i]*1e12,c='tab:blue',alpha=0.1)
plt.xlabel('Bias Frequency [Hz]')
plt.ylabel('Readout NEI [pA/rtHz]')
plt.title('Noise of '+str(n_sq_draw)+' SQ vs bias freq')
plt.savefig(path+'/nsq_'+str(nsq)+'_nei_lines.png')
plt.close()


plt.hist2d(np.array([bias_f for i in range(100)]).flatten(),   dfm.total.flatten()*1e12,
        cmap='Blues')
plt.xlabel('Bias freq [Hz]')
plt.colorbar()
plt.ylabel('Readout noise [pA/rtHz]')
plt.title('Estimated noise for '+str(n_sq_draw)+' SAA w/'+str(len(bias_f))+' ch each')
plt.savefig(path+'/nsq_'+str(nsq)+'_nei_hist.png')
plt.close()

plt.hist(dfm.total.flatten()*1e12)
plt.xlabel('Readout noise [pA/rtHz]')
plt.title(title+'Readout noise histogram')

plt.savefig(path+'/nsq_'+str(nsq)+'_nei_hist.png')
plt.close()


for band in range(len(opts)):
    popt = opts[band]
    


    nep_dat = (litebirb.nei_to_nep(band)*dfm.total)*1e18
    plt.hist(nep_dat.flatten(), alpha=0.5, bins=np.linspace(np.min(nep_dat),np.max(nep_dat),20))
    plt.xlabel('Readout NEP [aW/rtHz]')
    plt.title(title+'Band '+str(litebirb.opt_freqs[band])+' GHz')
    plt.savefig(path+'/band_'+str(band).zfill(2) +  '_nsq_'+str(nsq)+'_nep_hist.png')
    plt.close()


    np.save(path+'/band_'+str(band).zfill(2)+ '_nsq_'+str(nsq)+'_nep_dat.npy', nep_dat)