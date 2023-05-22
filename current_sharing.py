import numpy as np
from scipy import signal


from PySpice.Spice.Netlist import Circuit, SubCircuit, SubCircuitFactory
from PySpice.Unit import *


#a subcircuit that is just a bolometer with a LC resonator
class rlc(SubCircuit):
    NODES = ('rlc_in', 'out','inter')
    def __init__(self, name, R=1@u_Ω, L=60@u_uH, C = 100@u_pF, para = 1 @u_pF):
        SubCircuit.__init__(self, name, *self.NODES)
        self.R(1, 'rlc_in', 'rc', R)
        self.L(1, 'cl', 'll', L/2)
        self.L(2, 'll', 'out', L/2)
        self.C(1, 'rc', 'cl', C)
        
        
        #self.C(3, 'rc', 'inter', 0.3e-12 @u_F)
        #self.C(4, 'cl', 'inter', 0.3e-12 @u_F)
        #self.C(5, 'll', 'inter', 1.28 @u_pF)
        
        #self.R(5,'inter', self.gnd, 15 @u_Ohm)
        
#a subcircuit that is a segment of twisted pair wire
class wire(SubCircuit):
    NODES = ('wire+_in', 'wire-_in', 'wire+_out', 'wire-_out')
    def __init__(self, name, R=1@u_Ω, L=1@u_uH, C = 100@u_pF):
        SubCircuit.__init__(self, name, *self.NODES)
        self.R(1, 'wire+_in', 'rl+', R)
        self.L(1, 'rl+', 'wire+_out', L)
        self.C(1, 'rl+', 'rl-', C)
        
        self.R(2, 'wire-_in', 'rl-', R)
        self.L(2, 'rl-', 'wire-_out', L)
        
        



#something just for diagnositcs and sanity checks
#should produce the responce you'd expect from a network analysis
#refered to the SAA input
def get_network_analysis(dfmux_noise):
    carrier, nuller = make_circuits(dfmux_noise)
    dfmux_noise.spice = {'carrier': carrier, 'nuller': nuller}
    
    
    simulator = carrier.simulator(temperature=1,nominal_temperature=1)

    cna = simulator.ac(start_frequency = dfmux_noise.f[0]/1e6*0.9@u_MHz, stop_frequency = dfmux_noise.f[-1]/1e6*1.1@u_MHz, 
                       number_of_points=(1.1*dfmux_noise.f[-1]-0.9*dfmux_noise.f[0])/36, variation = 'lin')
    
    
    simulator2 = nuller.simulator(temperature=1,nominal_temperature=1)

    nna = simulator2.ac(start_frequency = dfmux_noise.f[0]/1e6*0.9@u_MHz, stop_frequency = dfmux_noise.f[-1]/1e6*1.1@u_MHz, 
                        number_of_points=(1.1*dfmux_noise.f[-1]-0.9*dfmux_noise.f[0])/36, variation = 'lin')
    
    
    return cna, nna



#function that produces the pyspice circuit
#one version with a voltage source across the bias element
#another with a current source across the SAA input
def make_circuits(dfmux_noise):
    cs = [1/(2 * np.pi * f)**2 / 60e-6 /1e-12 for f in dfmux_noise.f]
    carrier = make_circuit(dfmux_noise,cs,'car')
    nuller = carrier.clone(title='nul')   
    
    
    carrier.SinusoidalVoltageSource('bias','bias_pos','bias_neg',amplitude=1@u_V)
    
    
    
    #nuller = make_circuit(dfmux_noise,cs,'nul')
       
    nuller.SinusoidalCurrentSource('nuller','wire+'+str(9), 'wire-'+str(9) ,amplitude=1@u_A)
    
    return carrier, nuller

def make_circuit(dfmux_noise,cs,cname):

    #if only a single bolo r or rstray is specified make in an array of the number repeated
    if len(dfmux_noise.bolo.r)== 1:
        dfmux_noise.bolo.r  = dfmux_noise.bolo.r * np.ones(len(cs))
    if len(dfmux_noise.bolo.rstray)== 1:
        dfmux_noise.bolo.rstray  = dfmux_noise.bolo.rstray * np.ones(len(cs))


    nuller = Circuit(cname)
    i=0
    for c in cs:
        name= 'rlc'+str(i)
        nuller.subcircuit(rlc(name,R=(dfmux_noise.bolo.r[i] + dfmux_noise.bolo.rstray[i] )@u_Ohm, 
                              C=c@u_pF, para = dfmux_noise.para.c_gnd/len(cs)/4 @u_F))
        nuller.X(str(i+30),name,'bias_pos','sqin_pos','rlc_inter')
        i+=1
    #nuller.R('LC_gdn','rlc_inter',nuller.gnd,150 @u_Ohm)
        
        
    nuller.R('bias','bias_pos','bias_neg',30@u_mOhm)
    nuller.L('sqin','sqin_pos','bias_neg',dfmux_noise.squid.lin[0] @u_H)
    if dfmux_noise.squid.snubber != False:
        if dfmux_noise.squid.snubber_c != False:
            nuller.C('snubc','sqin_pos','snub_mid',dfmux_noise.squid.snubber_c @u_F)
            nuller.R('snubr','snub_mid','bias_neg',dfmux_noise.squid.snubber @u_Ohm)
        else:
            nuller.R('snubr','sqin_pos','bias_neg',dfmux_noise.squid.snubber @u_Ohm)
        
    #nuller.C('wafer_para','bias_pos',nuller.gnd,dfmux_noise.para.c_gnd @u_F)
    #nuller.C('p_para','sqin_pos',nuller.gnd,dfmux_noise.para.c_gnd @u_F)
    
    for i in range(10):
        name = 'wire'+str(i)
        nuller.subcircuit(wire(name, R= dfmux_noise.wire.r/10@u_Ohm, C=dfmux_noise.wire.c/10@u_F, L=dfmux_noise.wire.l/10@u_H))
        if i == 0:
            nuller.X(str(i),name,'sqin_pos', 'bias_neg', 'wire+'+str(i), 'wire-'+str(i))
        else:
            nuller.X(str(i),name,'wire+'+str(i-1), 'wire-'+str(i-1), 'wire+'+str(i), 'wire-'+str(i))
            
    nuller.R('R48', 'wire-'+str(i),nuller.gnd,dfmux_noise.para.r48 @u_Ohm)
    
    return nuller


def get_csf(dfmux_noise):
    cna, nna = get_network_analysis(dfmux_noise)
    
    bias_fs = signal.argrelmax(np.abs(cna.branches['lsqin']/nna.branches['lsqin']))[0]
    
    csf = np.array([1/np.abs(nna.branches['lsqin'])[i] for i in bias_fs])

    csf = np.repeat(csf.reshape(1,len(bias_fs)),len(dfmux_noise.squid.zt),axis=0)
    return csf






