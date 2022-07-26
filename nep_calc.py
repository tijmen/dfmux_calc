import dfmux_calc
import numpy as np
import matplotlib.pyplot as plt

colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']


class experiment:
    def __init__(self,name,dfmux_noise):
        if name == 'litebird':
            #optical power expected in each band as defined in IMo V1.2
            self.popt = [0.29184686,0.24187821,0.26855768,0.30598544,0.27094372,0.2958127,0.3278966,0.31625242,0.37653296,0.32719336,
                         0.30607807,0.35566829,0.35588565,0.43860747,0.42063367,0.39077456,0.35717258,0.62886291,0.4708167,0.37702823,
                         0.2996769,0.22026075]
            
            #assumes that the saturation power of a liteBIRD bolometer is tuned to be 2.5 times the optical power in the band
            self.psat = [ 2.5 * popt for popt in self.popt]
            
            #requirements on readout noise in aW/rtHz in each band
            self.rreq = [2.67399478,2.43434099,2.5650854,2.73799985,2.57645516,2.69210163,2.83433708,
                         2.78355604,3.03727798,2.83129603,2.73841425,2.95192701,2.95282887,3.27809138,3.21022195,3.09418433,
                         2.95816298,3.92519316,3.39632288,3.03927487,2.70962803,2.32301277]
            self.ntot = [8.31475732,7.83616637,8.64772011,8.63963629,8.50986423,9.36124323,9.5290737,9.72258017,11.27231753,9.69276546,9.77778099,
                         11.44233862,10.64069174,12.25155928,12.51592665,12.67883396,12.77011888,17.19246695,15.76686961,15.03623678,
                         14.41789678,13.34004439]
            self.nphon = [5.44441647,5.26100438,5.97783229,5.71954248,5.81005536,6.58377923,6.57637918,6.85404996,8.18115868,6.75923232,
                        6.97036698,8.45472357,7.61933154,8.9189256,9.26864051,9.55735219,9.78367735,13.22647082,12.31481045,11.91261558,
                          11.58180607,10.84704975]
            self.nphot = [3.94755931,3.59376372,3.7867788,4.04204858,3.80356373,3.97429005,4.1842691,4.10930217,4.48386626,
                          4.17977967,4.04266035,4.35786451,4.35919591,4.83937375,4.73917962,4.56787585,4.36707054,
                          5.7946758,5.0139163,4.48681422,4.00016388,3.42941234]
            #center frequency of each optical band in GHz
            self.opt_freqs = [40,60,78,50,68,89,68,89,119,78,100,140,100,119,140,166,195,195,235,280,337,402]
            
        #define any other experiments here...
        
        
        #stores the dfmux_noise object associated with the experiment
        self.readout = dfmux_noise
        
        
    def get_nep(self):
        self.nep = np.ones(len(self.opt_freqs))
        self.nep_range = np.ones((2,len(self.opt_freqs)))
        for i in range(len(self.opt_freqs)):
            self.nep[i] = self.nei_to_net(i)
            
        
        
    #function to take dfmux_noise object and experiment object
    def nei_to_nep(self,i,err=0):
        #err should be 0 if you want the target psat, and +1 or -1 if you want the psat to vary by a factor of 10%
        #estimating the voltage bias applied to the bolometer
        vbias = np.sqrt((self.readout.bolo.r + self.readout.bolo.rstray) * ( self.psat[i]*(1+err/10) -  self.popt[i] ) )
            
        #estimating the responsivity of the bolometer
        responsivity = np.sqrt(2)/vbias * self.readout.bolo.loopgain / ( 1 + 
              self.readout.bolo.loopgain *(self.readout.bolo.r - self.readout.bolo.rstray)/(self.readout.bolo.r + self.readout.bolo.rstray) )
        return 1/responsivity
    
    '''
    def plt_nei_req(self):
        resp = [self.nei_to_nep(i) for i in range(len(self.popt))]
        
        plt.figure()
        
        plt.scatter(self.opt_freqs, self.rreq, color=colors[3])
        for nei in np.range(5.5e-12, 10e-12, 5):
            plt.plot(self.opt_freqs, nei * resp, color=s
    '''
    
    
    
def plt_r_reqs(psat, popt, loopgain,nep):
    rbolo, rstray = np.meshgrid(np.linspace(0.1, 1.0 , 100), np.linspace(0, 0.2, 100))
    
    vbias = np.sqrt( (rbolo + rstray) * (psat - popt) )
    
    resp = np.sqrt(2) / vbias * loopgain / (1 + loopgain * (rbolo - rstray)/(rbolo + rstray ) ) 
    
    loop_atten = (rbolo - rstray)/(rbolo + rstray )
    
    nei_req = nep * resp
    
    fig, ax = plt.subplots()

    c = ax.pcolormesh(rbolo, rstray, nei_req, cmap='jet',vmin=5e-12,vmax = 15e-12)
    CS = ax.contour(rbolo, rstray, loop_atten, 6, colors='k') 
    ax.clabel(CS, fontsize=9, inline=True)
    ax.set_title('NEI requirement for $P_{{sat}}$={}, $P_{{opt}}=${}, $\mathcal{{L}}=${}, NEP$_{{read}}$={}'.format(
                    round(psat*1e12,2), popt*1e12, loopgain, nep))
    # set the limits of the plot to the limits of the data
    #ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)
    
    plt.xlabel('$R_{bolo}$ [$\Omega$]')
    plt.ylabel('$R_{stray}$ [$\Omega$]')
