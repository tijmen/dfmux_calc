#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import constants as c


#helper class to store SAA information
class squid:
    def __init__(self, zt, rdyn, inoise, lin, n_series=False, n_parallel = False, power = False):
        self.zt = zt            #Transimpedance of SAA
        self.rdyn = rdyn        #Dynamic impedance of SAA
        self.inoise = inoise    #SAA noise refered to input coil
        self.lin = lin          #input inductance of the SAA
        self.n_series = n_series #number of individual SQUIDs in series to form the SAA
        self.n_parallel = n_parallel #number of banks of SQUIDs in parallel to form the SAA
        self.power = power      #power dissipated by SAA when in operation

    #method to rescale the existing SAA to a new array
    def scale_SAA(self,new_series, new_parallel):
        if not self.n_series and not self.n_parallel and not self.power:
            print('squid has no specified array size and/or power!')
            return
        self.zt *= new_series / self.n_series
        self.rdyn *= new_series / new_parallel * self.n_parallel / self.n_series
        self.lin *= new_series * new_parallel / self.n_series / self.n_parallel
        self.inoise *= np.sqrt(self.n_series * self.n_parallel ) / np.sqrt(new_series * new_parallel)
        self.power *= new_series * new_parallel / self.n_series / self.n_parallel
        self.n_series = new_series
        self.n_parallel = new_parallel

#helper class to store bolometer information
class bolo:
    def __init__(self,r,loopgain,rstray,psat,tc,tb):
        self.r = r                    #Operating resistance of the bolometer in ohms
        self.loopgain = loopgain      #operating loopgain of the bolometer
        self.rstray = rstray          #stray resistance in series with the bolometer in ohms
        self.psat = psat              #saturation power of the bolometer in watts
        self.tc = tc                  #critical temperature of the bolometer in kelvin
        self.tb = tb                  #bath temperature the bolometer is operated at in kelvin
        
#helper class to store important parasitics
class parasitics:
    def __init__(self,stripline, c_gnd, r48=0):
        self.stripline = stripline    #self inductance in henries of the stripline conecting the SAA to the rest of the 
                                      #DfMUX circuit - or in its absence and stray inductance that looks like it
        self.c_gnd = c_gnd            #the parasitic capacitance in farads to ground from the focal plane, bolometers or LCs
        self.r48 = r48                #The value in ohms of R48 on the SQUID Controller Boards - this connected the negative
                                      #side of the nuller to ground and in 2022 its default value is 0 Ohms

#helper class to store wiring harness information and do some calculations    
class wire:
    def __init__(self, r, c, l, rshunt=False, cshunt=False):
        self.r = r              #wiring harness series resistance
        self.c = c              #capacitance
        self.l = l              #series self inductance
        self.rshunt = rshunt    #if a shunt resistor is places across the SAA output on the SQCB
        self.cshunt = cshunt    #if a shunt capacitor is placed across the SAA output on the SQCB
        
        
    #Calculate the series impedance of the wiring harness at a single given frequency
    #deceause of how transmission line models typically count wire resistance + inductance this is equivalent to
    #1 single wire of double length
    def series_imp(self,f):
        return 2j * np.pi * f * self.l + self.r
    
    
    #this calculates the complex ABCD or cascade matrix for the wiring harness using the transmision line model
    #at a given single frequency, taking into account resistive or capacitive shunts which are assumed to be
    #at the 300K input to a SQUID Controller board
    def get_abcd(self,f):
        wire_series = 2j * np.pi * f * self.l + self.r
        wire_shunt = 2j * np.pi * f * self.c
        
        gamma = np.sqrt(wire_series * wire_shunt)
        z0 = np.sqrt(wire_series / wire_shunt)
        
        b = z0 * np.sinh(gamma)
        if self.rshunt != False:
            a = np.cosh(gamma) + b/self.rshunt
        else:
            a = np.cosh(gamma)
        if self.cshunt != False:
            a = a + b * 2j*np.pi*f*self.cshunt
            
        
        d = np.cosh(gamma)
        if self.rshunt != False:
            c = 1/z0 * np.sinh(gamma) + d/self.rshunt
        else:
            c = 1/z0 * np.sinh(gamma)
        if self.cshunt != False:
            c = c + d * 2j*np.pi*f*self.cshunt
        return a, b, c, d
    
    
    #Calculates the effective resistance seen by the input of the 1st stage amplifier on the SQUID Controller board
    #this is used to calculate how the input current noise of that amplifer translates to voltage noise
    #at a given frequency and for a given squid with a dynamic impedance
    def reff(self,squid,f):
        [a,b,c,d] = self.get_abcd(f)
        return np.abs( (b + d * squid.rdyn)/( a + c * squid.rdyn))
    
    #Calculates the real effective resistance seen by the input of the 1st stage amplifier on the SQUID Controller board
    #without the SAA dynamic impedance, this is used to calculate the johnson noise of the wiring harness and any shunt elements
    def real_reff(f):
        [a,b,c,d] = self.get_abcd(f)
        return np.real( (b )/( a ))
        
            
    #calculates the voltage transfer function from the mK output of the SAA to the input of the 300K SQCB
    #for a list of frequencies and SAA with a dynamic impedance and stores the transfer function as a property of the wire
    def transfer_function(self, squid, frequencies):
        self.tf =[]        
        for f in frequencies:
            [a,b,c,d] = self.get_abcd(f) 
            self.tf.append( np.abs(1/(a + c * (squid.rdyn  ) - ( b + d * (squid.rdyn  ))/1e6 )) )
        return self.tf
    
    
    
    #Calculates the nuller(/carrier transfer function with inductive bias elements) from the wiring harness 
    #returns a ratio of current cold to current warm
    def nuller_tf(self,inductor,frequencies):
        itf = []
        for f in frequencies:
            
            [a,b,c,d] = self.get_abcd(f) 
            itf.append( np.abs(1/(c*2j*np.pi*f*inductor-d)))
        return itf



class dfmux_noise:
    
    def __init__(self, squid, bolo, wire, para, nuller_wire=None):
        self.squid = squid         #a squid object 
        self.bolo = bolo           #a bolometer object
        self.wire = wire           #a wiring harness object
        self.para = para           #a parasitics object
        self.nuller_wire = nuller_wire #a wiring harness object for the nuller lines if this is different than the ones used for the
                                   #SQUID output lines. By default this is None and it is assumed they are the same.
        
        #Warm electrnics noise
        carrier = 1.6e-12                                        #A/rtHz JM PhD Table 7.5
        nuller = 4.9e-12                                         #A/rtHz JM PhD Table 7.6
        self.warm_noise_nc = np.sqrt(carrier**2 + nuller**2)     #total noise from the carrier/nuller refered to SAA input
        

        #bolometer noise
        self.jnoise = np.sqrt(2) * 1/(1+self.bolo.loopgain)*np.sqrt(4*1.38e-23*self.bolo.tc / (self.bolo.r))  #JM masters section 5.6
        
        
    def init_freq(self, frequencies, #frequencies to calculate noise at
                  dan=True,          #if this is a dan on noise calculation or not
                  csf=None):         #if you want to calculate noise with a given csf (this must be the same size as frequencies)
        
        
        self.csf = []                    #current sharing factor
        self.tf = []                     #wiring harness transfer function
        self.saa_in_impedance = []       #impedance of the SAA imput coil
        self.on_res_comb_impedance = []  #impedance of the comb assuming on resonance
        self.c_r48 = []                  #impedance of the parasitic path through ground and R48
        self.demod = []                  #total demodulator noise refered to the input of the SAA A/rtHz
        self.total = []                  #total noise refered to the input of the SAA pA/rtHz
        #self.wire_j = []
        if not dan:
            self.jnoise_f = []
            self.warm_noise_nc_f = []
            self.saa_scale_f = []
        
        
        i=0
        for f in frequencies:
            #demod chain noise
            #see JM PhD table 7.7 1st stage amplifier
            if self.wire.rshunt==False:
                req=(1/(1/10 + 1/100 + 1/150) + 1/(1/4.22e3 + 1/ (self.wire.reff(self.squid,f))) ) #RSQCB as defined in table 7.7
                first_amp = np.sqrt(2) * np.sqrt((1.1e-9)**2  +
                                (2.2e-12* req)**2)
            
            
            #this is a modification to the previous line that adds the johnson noise of the shunt resistor
            #currently assuming that it is well represented by a current source in paraelle with the 1st stage amplifier
            #which adds in quadrature with that current noise and is then scaled by the impedance as defined in Joshua's thesis
            #the change in Reff due to the shunt will be automaticaly handled by the reff function
            else:
                req=(1/(1/10 + 1/100 + 1/150) + 1/(1/4.22e3 + 1/ (self.wire.reff(self.squid,f))) ) #RSQCB as defined in table 7.7
                first_amp = np.sqrt(2) * np.sqrt(   (1.1e-9)**2  +
                                ( np.sqrt(  (2.2e-12)**2 +  4* 1.38e-23 *300/self.wire.rshunt  ) * req)**2
                                    #+(4* 1.38e-23 *300 / (self.wire.r+self.wire.rshunt)*(self.wire.reff(self.squid,f))**2))
                                    #+(4* 1.38e-23 *300 * (self.wire.r+self.wire.rshunt))
                                   #(2.2e-12 * req)**2 +  
                                   #        1.38e-23 *300/self.wire.rshunt * (self.wire.reff(self.squid,f))**2 
                                                )

            #table 7.7- the rest of the terms
            demod_dc = np.sqrt(first_amp**2 + 
                               2*(0.23e-9**2 +     #ADC noise
                                  0.14e-9**2 +     #2nd stage amplifier
                               #0.36e-9**2 +        #Signal path johnson - think this is just wiring harness johnson
                                  (8.36e-9*self.wire.reff(self.squid,f)/(self.wire.reff(self.squid,f)+4.22e3))**2   #SQUID I bias johnson
                                  + 4* 1.38e-23 *300 * self.wire.r ))  #johnson noise of the wiring harness resistance
                                                                      #pessimisitcally assuming that the entire harness is at 300K
            
            
            
            if dan:
                #in this case estimate the current sharing factor based off of JM SPT-3G noise paper
                if csf == None:
                    ##Calculate the impedances necesary for current sharing
                    
                    #the input impedance of the SAA
                    self.saa_in_impedance.append(  2 * np.pi * f * self.squid.lin  )
                    
                    #the comb impedance- assuming this is for an on resonance frequency
                    self.on_res_comb_impedance.append(  2 * np.pi * f * self.para.stripline + self.bolo.r + self.bolo.rstray  )
                    
                    #the impedance of the path through ground and R48
                    self.c_r48.append(  1/(2 * np.pi * f * self.para.c_gnd) + self.para.r48  )
                    
                    #combining those terms to estimate the current sharing factor- from joshua's spt3g paper
                    if self.nuller_wire == None:
                        self.csf.append( 1/(self.c_r48[-1] / ( (self.on_res_comb_impedance[-1] + np.abs(self.wire.series_imp(f))) * 
                                         self.saa_in_impedance[-1] / 
                           (self.on_res_comb_impedance[-1] + self.saa_in_impedance[-1] + np.abs(self.wire.series_imp(f))) + 
                           np.abs(self.wire.series_imp(f)) + self.c_r48[-1] ) * 
                                self.on_res_comb_impedance[-1] / (self.on_res_comb_impedance[-1] + self.saa_in_impedance[-1])))
                    
                    #if nuller_wire is not none this has a different wiring than connected to the SAA output
                    else:
                        self.csf.append( 1/(self.c_r48[-1] / ( (self.on_res_comb_impedance[-1] + np.abs(self.nuller_wire.series_imp(f))) * 
                                         self.saa_in_impedance[-1] / 
                           (self.on_res_comb_impedance[-1] + self.saa_in_impedance[-1] + np.abs(self.nuller_wire.series_imp(f))) + 
                           np.abs(self.nuller_wire.series_imp(f)) + self.c_r48[-1] ) * 
                                self.on_res_comb_impedance[-1] / (self.on_res_comb_impedance[-1] + self.saa_in_impedance[-1])))
                
                #in this case take the measured CSF given as input
                else:
                    self.csf.append(csf[i])
            
            #this is an option to calculate SQCB input refered noise with DAN off forcing CSF to not impact the noise
            else:
                self.csf.append(1)
                
                
                
            #calculating the transfer function caused by the SAA Z_dyn and the wiring harness capacitance/any shunts across SAA
            self.tf.append(self.wire.transfer_function(self.squid,[f])[0])
            
            
            #scaling the noise from the warm electronics by the current sharing factor, the transfer function
            # and the transimpedance of the saa to refer it to the SAA input coil - units now A/rtHz
            if dan:
                self.demod.append( demod_dc * self.csf[-1] / self.tf[-1] / self.squid.zt )
            else:
                self.demod.append( demod_dc  )
                
            #self.wire_j.append(np.sqrt(4* 1.38e-23 *300 * self.wire.real_reff(self.squid,f)*np.sqrt(2)*self.csf[-1]/self.squid.zt))

            
            #scaling the noise of the SAA by the current sharing factor and the demodulation factor
            if dan:
                saa_scale = self.squid.inoise * self.csf[-1] * np.sqrt(2)
            #or if dan off noise by the transimpedance and tf to refer to SQCB
            else:
                self.saa_scale_f.append(self.squid.inoise * self.squid.zt * self.tf[-1] * np.sqrt(2))
            
            #total noise from all sources in order the noise from the carrier/nuller chain, the johnson noise of the bolo
            #the scaled demodulator chain noise, and the scaled SAA noise in PICOAMPS/rtHz
            if dan:
                self.total.append(np.sqrt(self.warm_noise_nc**2 + self.jnoise**2 + self.demod[-1]**2 + saa_scale**2
                                   )*1e12)
            
            #if dan is off refering remaining terms to SQCB input - outputs noise in nV/rtHz
            else:
                self.warm_noise_nc_f.append(self.warm_noise_nc* self.squid.zt* self.tf[-1])
                self.jnoise_f.append(self.jnoise* self.squid.zt* self.tf[-1])
                self.total.append(np.sqrt((self.warm_noise_nc_f[-1])**2  + (self.jnoise_f[-1])**2 + 
                                          self.demod[-1]**2 + self.saa_scale_f[-1]**2)*1e9)
                
            i+=1




#helper function which takes a SQUID object as input and returns another SQUID object as output
#it takes warm transimpecence dyn impedance etc measurements dont with a resistive shunt in place
#and calculates the expected parameters for the SQUID on its own, assuming the measurements were done
#with a 1.6 uA current bias step as is pydfmux's default behavior
def refer_squid(warm_squid, shunt):
    delta_v = warm_squid.rdyn*1.6e-6
    cold_zdyn = delta_v / (1.6e-6 - delta_v /shunt)
    
    cold_zt = warm_squid.zt * (cold_zdyn + shunt ) / shunt
    
    return squid(cold_zt, cold_zdyn, warm_squid.inoise, warm_squid.lin)
  





  
        
