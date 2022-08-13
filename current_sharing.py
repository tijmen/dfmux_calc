import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
import numpy as np
from scipy import signal


from PySpice.Spice.Netlist import Circuit, SubCircuit, SubCircuitFactory
from PySpice.Unit import *


#a subcircuit that is just a bolometer with a LC resonator
class rlc(SubCircuit):
    NODES = ('in', 'out')
    def __init__(self, name, R=1@u_Ω, L=60@u_uH, C = 100@u_pF):
        SubCircuit.__init__(self, name, *self.NODES)
        self.R(1, 'in', 'rl', R)
        self.L(1, 'rl', 'lc', L)
        self.C(1, 'lc', 'out', C)
        
        



#something just for diagnositcs and sanity checks
#should produce the responce you'd expect from a network analysis
#refered to the SAA input
def get_network_analysis(dfmux_noise):
    carrier, nuller = make_circuits(dfmux_noise)
    
    
    simulator = carrier.simulator(temperature=1,nominal_temperature=1)

    cna = simulator.ac(start_frequency = dfmux_noise.f[0]/1e6*0.9@u_MHz, stop_frequency = dfmux_noise.f[-1]/1e6*1.1@u_MHz, number_of_points=100000, variation = 'lin')
    
    
    simulator2 = nuller.simulator(temperature=1,nominal_temperature=1)

    nna = simulator2.ac(start_frequency = dfmux_noise.f[0]/1e6*0.9@u_MHz, stop_frequency = dfmux_noise.f[-1]/1e6*1.1@u_MHz, number_of_points=100000, variation = 'lin')
    
    
    return cna, nna



#function that produces the pyspice circuit
#one version with a voltage source across the bias element
#another with a current source across the SAA input
def make_circuits(dfmux_noise):
    cs = [1/(2 * np.pi * f)**2 / 60e-6 /1e-12 for f in dfmux_noise.f]
    carrier = Circuit('car')
    i=0
    for c in cs:
        name= 'rlc'+str(i)
        carrier.subcircuit(rlc(name,R=dfmux_noise.bolo.r@u_Ohm,C=c@u_pF))
        carrier.X(str(i),name,'bias_pos','sqin_pos')
        i+=1
    carrier.R('bias','bias_pos','bias_neg',30@u_mOhm)
    carrier.L('sqin','sqin_pos','bias_neg',dfmux_noise.squid.lin @u_H)
    if dfmux_noise.squid.snubber != False:
        carrier.R('snubr','sqin_pos','bias_neg',dfmux_noise.squid.snubber @u_Ohm)
    if dfmux_noise.squid.snubber_c != False:
        carrier.C('snubc','sqin_pos','bias_neg',dfmux_noise.squid.snubber_c @u_F)
    carrier.R('short','bias_neg',carrier.gnd,0@u_mOhm)
    carrier.SinusoidalVoltageSource('bias','bias_pos','bias_neg',amplitude=1@u_V)
    
    
    
    nuller = Circuit('nul')
    i=0
    for c in cs:
        name= 'rlc'+str(i)
        nuller.subcircuit(rlc(name,R=dfmux_noise.bolo.r@u_Ohm, C=c@u_pF))
        nuller.X(str(i),name,'bias_pos','sqin_pos')
        i+=1
    nuller.R('bias','bias_pos','bias_neg',30@u_mOhm)
    nuller.L('sqin','sqin_pos','bias_neg',dfmux_noise.squid.lin @u_H)
    if dfmux_noise.squid.snubber != False:
        nuller.R('snubr','sqin_pos','bias_neg',dfmux_noise.squid.snubber @u_Ohm)
    if dfmux_noise.squid.snubber_c != False:
        nuller.C('snubc','sqin_pos','bias_neg',dfmux_noise.squid.snubber_c @u_F)
    nuller.R('short','bias_neg',nuller.gnd,0@u_mOhm)
    nuller.SinusoidalCurrentSource('nuller','sqin_pos','bias_neg',amplitude=1@u_A)
    
    return carrier, nuller


def get_csf(dfmux_noise):
    cna, nna = get_network_analysis(dfmux_noise)
    
    bias_fs = signal.argrelmax(np.abs(cna.branches['lsqin']/nna.branches['lsqin']))[0]
    
    csf = np.array([1/np.abs(nna.branches['lsqin'])[i] for i in bias_fs])
    
    return csf
    