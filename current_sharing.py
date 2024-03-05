import numpy as np
from scipy import signal


from PySpice.Spice.Netlist import Circuit, SubCircuit, SubCircuitFactory
from PySpice.Unit import *


# a subcircuit that is just a bolometer with a LC resonator
class rlc(SubCircuit):
    NODES = ("rlc_in", "out", "inter")

    def __init__(self, name, R=1 @ u_Ω, L=60 @ u_uH, C=100 @ u_pF, para=1 @ u_pF):
        SubCircuit.__init__(self, name, *self.NODES)
        self.R(1, "rlc_in", "rc", R)
        self.L(1, "cl", "ll", L / 2)
        self.L(2, "ll", "out", L / 2)
        self.C(1, "rc", "cl", C)

        # self.C(3, 'rc', 'inter', 0.3e-12 @u_F)
        # self.C(4, 'cl', 'inter', 0.3e-12 @u_F)
        # self.C(5, 'll', 'inter', 1.28 @u_pF)

        # self.R(5,'inter', self.gnd, 15 @u_Ohm)


# a subcircuit that is a segment of twisted pair wire
class wire(SubCircuit):
    NODES = ("wire+_in", "wire-_in", "wire+_out", "wire-_out")

    def __init__(self, name, R=1 @ u_Ω, L=1 @ u_uH, C=100 @ u_pF):
        SubCircuit.__init__(self, name, *self.NODES)
        self.R(1, "wire+_in", "rl+", R)
        self.L(1, "rl+", "wire+_out", L)
        self.C(1, "rl+", "rl-", C)

        self.R(2, "wire-_in", "rl-", R)
        self.L(2, "rl-", "wire-_out", L)


def get_network_analysis(dfmux_noise):
    carrier, nuller = make_circuits(dfmux_noise)
    dfmux_noise.spice = {"carrier": carrier, "nuller": nuller}

    simulator = carrier.simulator(temperature=1, nominal_temperature=1)

    f = np.atleast_1d(dfmux_noise.f)
    print(f"Input frequency(ies): {f}")

    if len(f) == 1:
        center_frequency = f[0] / 1e6 @ u_MHz
        start_frequency = center_frequency * 0.9
        stop_frequency = center_frequency * 1.1
        number_of_points = 101
        print(f"Single frequency input: {center_frequency}")
        print(
            f"Simulating from {start_frequency} to {stop_frequency} with {number_of_points} points"
        )
    else:
        start_frequency = f.min() / 1e6 * 0.9 @ u_MHz
        stop_frequency = f.max() / 1e6 * 1.1 @ u_MHz
        frequency_step = np.diff(f).min()
        number_of_points = int((1.1 * f.max() - 0.9 * f.min()) / frequency_step) + 1
        print(f"Multiple frequency input ranging from {f.min()} to {f.max()}")
        print(
            f"Simulating from {start_frequency} to {stop_frequency} with {number_of_points} points"
        )

    cna = simulator.ac(
        start_frequency=start_frequency,
        stop_frequency=stop_frequency,
        number_of_points=number_of_points,
        variation="lin",
    )

    simulator2 = nuller.simulator(temperature=1, nominal_temperature=1)

    nna = simulator2.ac(
        start_frequency=start_frequency,
        stop_frequency=stop_frequency,
        number_of_points=number_of_points,
        variation="lin",
    )

    # print("Verify outputs...")
    # print(f"Carrier branches: {cna.branches}")
    # print(f"Nuller branches: {nna.branches}")

    return cna, nna, simulator, simulator2


# function that produces the pyspice circuit
# one version with a voltage source across the bias element
# another with a current source across the SAA input
def make_circuits(dfmux_noise):
    cs = 1 / (2 * np.pi * dfmux_noise.f) ** 2 / 60e-6 / 1e-12
    carrier = make_circuit(dfmux_noise, cs, "car")
    nuller = carrier.clone(title="nul")

    carrier.SinusoidalVoltageSource("bias", "bias_pos", "bias_neg", amplitude=1 @ u_V)

    # nuller = make_circuit(dfmux_noise,cs,'nul')

    nuller.SinusoidalCurrentSource(
        "nuller", "wire+" + str(9), "wire-" + str(9), amplitude=1 @ u_A
    )

    return carrier, nuller


def make_circuit(dfmux_noise, cs, cname):
    nuller = Circuit(cname)

    # Convert cs to a numpy array to ensure consistent iteration and indexing
    cs_array = np.atleast_1d(cs)

    # Use numpy's broadcasting to handle scalar values of dfmux_noise.bolo.r and dfmux_noise.bolo.rstray
    r_array = np.atleast_1d(dfmux_noise.bolo.r)
    rstray_array = np.atleast_1d(dfmux_noise.bolo.rstray)

    # Ensure r_array and rstray_array are at least as long as cs_array, repeating their last element if necessary
    r_array = np.pad(r_array, (0, max(0, len(cs_array) - len(r_array))), "edge")
    rstray_array = np.pad(
        rstray_array, (0, max(0, len(cs_array) - len(rstray_array))), "edge"
    )

    for i, c in enumerate(cs_array):
        name = f"rlc{i}"
        nuller.subcircuit(
            rlc(
                name,
                R=(r_array[i] + rstray_array[i]) @ u_Ohm,
                C=c @ u_pF,
                para=dfmux_noise.para.c_gnd / len(cs_array) / 4 @ u_F,
            )
        )
        nuller.X(str(i + 30), name, "bias_pos", "sqin_pos", "rlc_inter")

    nuller.R("bias", "bias_pos", "bias_neg", 30 @ u_mOhm)
    nuller.L("sqin", "sqin_pos", "bias_neg", dfmux_noise.squid.lin @ u_H)
    if dfmux_noise.squid.snubber != False:
        if dfmux_noise.squid.snubber_c != False:
            nuller.C("snubc", "sqin_pos", "snub_mid", dfmux_noise.squid.snubber_c @ u_F)
            nuller.R("snubr", "snub_mid", "bias_neg", dfmux_noise.squid.snubber @ u_Ohm)
        else:
            nuller.R("snubr", "sqin_pos", "bias_neg", dfmux_noise.squid.snubber @ u_Ohm)

    # nuller.C('wafer_para','bias_pos',nuller.gnd,dfmux_noise.para.c_gnd @u_F)
    # nuller.C('p_para','sqin_pos',nuller.gnd,dfmux_noise.para.c_gnd @u_F)

    for i in range(10):
        name = "wire" + str(i)
        nuller.subcircuit(
            wire(
                name,
                R=dfmux_noise.wire.r / 10 @ u_Ohm,
                C=dfmux_noise.wire.c / 10 @ u_F,
                L=dfmux_noise.wire.l / 10 @ u_H,
            )
        )
        if i == 0:
            nuller.X(
                str(i), name, "sqin_pos", "bias_neg", "wire+" + str(i), "wire-" + str(i)
            )
        else:
            nuller.X(
                str(i),
                name,
                "wire+" + str(i - 1),
                "wire-" + str(i - 1),
                "wire+" + str(i),
                "wire-" + str(i),
            )

    nuller.R("R48", "wire-" + str(i), nuller.gnd, dfmux_noise.para.r48 @ u_Ohm)

    return nuller


def get_csf(dfmux_noise):
    print("Calling get_network_analysis...")
    cna, nna, sim, sim2 = get_network_analysis(dfmux_noise)

    carrier_na = cna.branches["lsqin"].as_ndarray()
    nuller_na = nna.branches["lsqin"].as_ndarray()
    print(f"Carrier NA: {carrier_na}")
    print(f"Nuller NA: {nuller_na}")

    print("Calculating ratio of lsqin branches...")
    ratio = np.abs(carrier_na / nuller_na)
    print(f"Ratio values: {ratio}")

    print("Finding bias_fs using signal.argrelmax...")
    bias_fs = signal.argrelmax(ratio)[0]
    print(f"bias_fs: {bias_fs}")

    print("Calculating csf...")
    csf = []
    for i in bias_fs:
        lsqin_value = np.abs(nna.branches["lsqin"].as_ndarray()[i])
        print(f"lsqin value at index {i}: {lsqin_value}")
        csf_value = 1 / lsqin_value
        csf.append(csf_value)
        print(f"CSF value at index {i}: {csf_value}")

    csf = np.array(csf)
    print(f"Final csf array: {csf}")

    print("Cleaning up simulation instances...")
    ngspice1 = sim.factory(cna).ngspice
    ngspice1.remove_circuit()
    ngspice1.destroy()
    ngspice2 = sim2.factory(nna).ngspice
    ngspice2.remove_circuit()
    ngspice2.destroy()

    return csf
