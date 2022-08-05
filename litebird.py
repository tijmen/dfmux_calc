import numpy as np
import matplotlib.pyplot as plt
import dfmux_calc as d
import nep_calc as n

#input parameters
lb = n.experiment('litebird',0)
#band = 1   #hard reqs
band = 11  #mid reqs
#band = 18  #easy reqs
frac = 0.05
popt = lb.popt[band]
psat = 2.5 * popt
loopgain = 10
nep = lb.ntot[band]  * np.sqrt((1+frac)**2 - 1)



#def r and stray to sweep through
rbolo, rstray = np.meshgrid(np.linspace(0.1, 1.0 , 100), np.linspace(0, 0.2, 100))
    
vbias = np.sqrt( (rbolo + rstray) * (psat - popt) )
    
resp = np.sqrt(2) / vbias * loopgain / (1 + loopgain * (rbolo - rstray)/(rbolo + rstray ) ) 
    
loop_atten = (rbolo - rstray)/(rbolo + rstray )
    
nei_req = nep * resp

req_power = nei_req.copy()
req_nsq = nei_req.copy()
fail = nei_req.copy()

#baseline SAA and other parasitics
#SA13 properties representative of the median SAA 
sa13 = d.squid(1750,                                     #Transimpedence [Ohms]
               350,                                     #Dynamic impedence [Ohms]
               1e-12,                                   #NEI [A/rtHz] (just taking the number JM used)
               70e-9,                                   #Input inductance [H]
               n_series=3*64,n_parallel=2,power=200e-9,  #array size and power dissipation 
               snubber=5,                           #if there is a snubber on the input
               t=0.3)                                     #what temperature the SAA is at [K]
#sa13.change_mutual_ind(3)
#Wiring harness properties
wh = d.wire(30,                            #resistance [Ohms]
            40e-12,                        #capacitance [F] (this is the important one)
            0.75e-6,                       #inductance  [H]
            rshunt=False, cshunt=False)    #if theres any resistive or capcitive shunts across the output of the SAA

#bolometer properties
bolo = d.bolo(1.0,            #operating impedence [Ohms]
              loopgain,              #loopgain 
              0.2,            #stray impedence [Ohms]
              5e-12,          #psat - this isn't used at all
              0.171,0.1)        #Tc and Tb [K]

#Other parasitics
para = d.parasitics(60e-9,0.7e-9,0) #stripline inductance[H], parasitic capacitance to ground[F] and R48 [Ohms]

dfm = d.dfmux_noise(sa13,bolo,wh,para,nuller_cold=True)

bias_f = [4.5e6]


#calc needed NEI at each combo
#step through each
for i in range(100):
    #print(i)
    for j in range(100):
        
        dfm.bolo.r = rbolo[i][j]
        dfm.bolo.rstray = rstray[i][j]
        target = nei_req[i][j]*1e12
        #print(target)
        n_sq = 1
        #increase number of SAA until NEI met
        while True:
            #print(i,j,n_sq)
            dfm.squid.scale_SAA(n_sq, 1)
            dfm.init_freq(bias_f)
            if max(dfm.total) <= target:
                req_power[i][j] = dfm.squid.power
                req_nsq[i][j] = n_sq
                fail[i][j] = 0
                #print(n_sq)
                break
            else:
                n_sq +=1
                if n_sq <=600:
                    continue
                else:
                    #print('fail!:',max(dfm.total))
                    req_power[i][j] = dfm.squid.power
                    req_nsq[i][j] = n_sq
                    fail[i][j] = 1
                    break
    
    #calculate power per mux
    #print(j)
    
    
#plot    
fig, ax = plt.subplots()

c = ax.pcolormesh(rbolo, rstray, req_nsq, cmap='jet',vmin=1,vmax = 100)
d = ax.pcolormesh(rbolo, rstray, np.where(fail == 1, 1, np.nan), cmap='binary', vmin=0, vmax=1,zorder=10)
CS = ax.contour(rbolo, rstray, loop_atten, 6, colors='w') 
ax.clabel(CS, fontsize=9, inline=True)
ax.set_title('# SQ requirement for {} GHz band: $P_{{sat}}$={}pW, $P_{{opt}}=${}pW, \n $\mathcal{{L}}=${}, NEP$_{{read}}$={}aW$/\sqrt{{\mathrm{{Hz}}}}$, {}% NEP increase'.format(
                    lb.opt_freqs[band] ,round(psat*1e12,2), round(popt*1e12,2), loopgain, round(nep*1e18,1), frac*100))
    # set the limits of the plot to the limits of the data
    #ax.axis([x.min(), x.max(), y.min(), y.max()])
cbar = fig.colorbar(c, ax=ax)
    
plt.xlabel('$R_{bolo}$ [$\Omega$]')
plt.ylabel('$R_{stray}$ [$\Omega$]')
cbar.set_label('Required SQUIDs in array')




fig, ax = plt.subplots()

c = ax.pcolormesh(rbolo, rstray, nei_req*1e12, cmap='jet',vmin=5,vmax = 15)
d = ax.pcolormesh(rbolo, rstray, np.where(fail == 1, 1, np.nan), cmap='binary', vmin=0, vmax=1,zorder=10)
CS = ax.contour(rbolo, rstray, loop_atten, 6, colors='w') 
ax.clabel(CS, fontsize=9, inline=True)
ax.set_title('NEI requirement for {} GHz band: $P_{{sat}}$={}pW, $P_{{opt}}=${}pW, \n $\mathcal{{L}}=${}, NEP$_{{read}}$={}aW$/\sqrt{{\mathrm{{Hz}}}}$, {}% NEP increase'.format(
                    lb.opt_freqs[band] ,round(psat*1e12,2), round(popt*1e12,2), loopgain, round(nep*1e18,1), frac*100))
    # set the limits of the plot to the limits of the data
    #ax.axis([x.min(), x.max(), y.min(), y.max()])
cbar = fig.colorbar(c, ax=ax)
    
plt.xlabel('$R_{bolo}$ [$\Omega$]')
plt.ylabel('$R_{stray}$ [$\Omega$]')
cbar.set_label('Required NEI [pA/$\sqrt{\mathrm{Hz}}$]')