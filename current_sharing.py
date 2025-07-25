import numpy as np
from scipy import signal

try:
    from PySpice.Spice.Netlist import Circuit, SubCircuit, SubCircuitFactory
    from PySpice.Unit import u_Ω, u_uH, u_pF, u_F, u_Ohm, u_V, u_A, u_mOhm, u_H, u_MHz
except:
    print("\nWarning! Can't find PySpice! Using analytic calculation for current sharing. \
        \nInstall PySpice to prevent this from happening. \n")
    loaded_pyspice = False
else:
    loaded_pyspice = True

    
if loaded_pyspice:

    class RLC(SubCircuit):
        """Subcircuit representing a bolometer with an LC resonator."""

        NODES = ("rlc_in", "out", "inter")

        def __init__(self, name, R=1 @ u_Ω, L=60 @ u_uH, C=100 @ u_pF, para=1 @ u_pF):
            super().__init__(name, *self.NODES)
            self.R(1, "rlc_in", "rc", R)
            self.L(1, "cl", "ll", L / 2)
            self.L(2, "ll", "out", L / 2)
            self.C(1, "rc", "cl", C)


    class Wire(SubCircuit):
        """Subcircuit representing a segment of twisted pair wire."""

        NODES = ("wire+_in", "wire-_in", "wire+_out", "wire-_out")

        def __init__(self, name, R=1 @ u_Ω, L=1 @ u_uH, C=100 @ u_pF):
            super().__init__(name, *self.NODES)
            self.R(1, "wire+_in", "rl+", R)
            self.L(1, "rl+", "wire+_out", L)
            self.C(1, "rl+", "rl-", C)
            self.R(2, "wire-_in", "rl-", R)
            self.L(2, "rl-", "wire-_out", L)


    def get_network_analysis(dfmux_system):
        """Performs network analysis using PySpice."""

        carrier, nuller = make_circuits(dfmux_system)
        dfmux_system.spice = {"carrier": carrier, "nuller": nuller}

        simulator = carrier.simulator(temperature=1, nominal_temperature=1)
        f = np.atleast_1d(dfmux_system.f)

        if len(f) == 1:
            center_frequency = f[0] / 1e6 @ u_MHz
            start_frequency = center_frequency * 0.98
            stop_frequency = center_frequency * 1.02
            number_of_points = 1001
        else:
            start_frequency = f.min() / 1e6 * 0.9 @ u_MHz
            stop_frequency = f.max() / 1e6 * 1.1 @ u_MHz
            frequency_step = np.diff(f).min()/36.
            number_of_points = int((1.1 * f.max() - 0.9 * f.min()) / frequency_step) + 1

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

        return cna, nna, simulator, simulator2


    def make_circuits(dfmux_system):
        """Creates PySpice circuits for carrier and nuller analysis."""

        cs = 1 / (2 * np.pi * dfmux_system.f) ** 2 / 60e-6 / 1e-12
        carrier = make_circuit(dfmux_system, cs, "car")
        nuller = carrier.clone(title="nul")

        carrier.SinusoidalVoltageSource("bias", "bias_pos", "bias_neg", amplitude=1 @ u_V)
        nuller.SinusoidalCurrentSource(
            "nuller", "wire+" + str(9), "wire-" + str(9), amplitude=1 @ u_A
        )

        return carrier, nuller


    def make_circuit(dfmux_system, cs, cname):
        """Creates a PySpice circuit with the given parameters."""

        nuller = Circuit(cname)

        cs_array = np.atleast_1d(cs)
        r_array = np.atleast_1d(dfmux_system.operating_resistance)
        rstray_array = np.atleast_1d(dfmux_system.stray_resistance)

        # Ensure consistent array lengths for iteration
        r_array = np.pad(r_array, (0, max(0, len(cs_array) - len(r_array))), "edge")
        rstray_array = np.pad(
            rstray_array, (0, max(0, len(cs_array) - len(rstray_array))), "edge"
        )

        for i, c in enumerate(cs_array):
            name = f"rlc{i}"
            nuller.subcircuit(
                RLC(
                    name,
                    R=(r_array[i] + rstray_array[i]) @ u_Ohm,
                    C=c @ u_pF,
                    para=dfmux_system.parasitic_capacitance / len(cs_array) / 4 @ u_F,
                )
            )
            nuller.X(str(i + 30), name, "bias_pos", "sqin_pos", "rlc_inter")

        nuller.R("bias", "bias_pos", "bias_neg", 30 @ u_mOhm)
        nuller.L("sqin", "sqin_pos", "bias_neg", dfmux_system.squid_input_inductance @ u_H)

        if dfmux_system.snubber:
            if dfmux_system.snubber_capacitance:
                nuller.C(
                    "snubc", "sqin_pos", "snub_mid", dfmux_system.snubber_capacitance @ u_F
                )
                nuller.R("snubr", "snub_mid", "bias_neg", dfmux_system.snubber @ u_Ohm)
            else:
                nuller.R("snubr", "sqin_pos", "bias_neg", dfmux_system.snubber @ u_Ohm)

        for i in range(10):
            name = "wire" + str(i)
            nuller.subcircuit(
                Wire(
                    name,
                    R=dfmux_system.wire_harness_resistance / 10 @ u_Ohm,
                    C=dfmux_system.wire_harness_capacitance / 10 @ u_F,
                    L=dfmux_system.wire_harness_inductance / 10 @ u_H,
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

        nuller.R("R48", "wire-" + str(i), nuller.gnd, dfmux_system.r48 @ u_Ohm)

        return nuller


    def get_csf(dfmux_system):
        """Calculates the current sharing factor using PySpice."""

        cna, nna, sim, sim2 = get_network_analysis(dfmux_system)

        carrier_na = cna.branches["lsqin"].as_ndarray()
        nuller_na = nna.branches["lsqin"].as_ndarray()

        ratio = np.abs(carrier_na / nuller_na)

        bias_fs = signal.argrelmax(ratio)[0]

        csf = []
        for i in bias_fs:
            lsqin_value = np.abs(nna.branches["lsqin"].as_ndarray()[i])
            csf_value = 1 / lsqin_value
            csf.append(csf_value)

        csf = np.array(csf)

        # Clean up simulation instances
        ngspice1 = sim.factory(cna).ngspice
        ngspice1.remove_circuit()
        ngspice1.destroy()

        ngspice2 = sim2.factory(nna).ngspice
        ngspice2.remove_circuit()
        ngspice2.destroy()

        return csf


def get_csf_analytically(dfmux_system, debug=False):
    """ 
    Author: Nicole Farias
    Calculates the current sharing factor analytically.
    Does everythin in terms of absolute impedances otherwise the 
    results get affected for so far unknown reasons
    
    Note: dmfux_system.f must be initialized with some np.array of bias frequencies
    before calling this dunction
    
    Inputs needed in dfmux_system
    f : bias frequency [Hz]
    r48 : the reference to ground resistor [omhs]
    parasitic_capacitance : parasitic capacitance to ground [F]
    wire_harness_resistance : used to calculate series wire harness impedance [omhs]
    wire_harness_inductance : used to calculate series wire harness impedance [H]
    squid_input_inductance : [H]
    operating_resistance : the TES operating resistance [omhs]
    stray_resistance : the stray resistance in series with the TES [omhs]
    stripline_inductance : stripline inductance - presumably where most of the series inductance comes from [H]
    
    """
    
    omega_b = 2*np.pi*dfmux_system.f 
    
    # impedance of bias line. 
    Z_b = np.abs(dfmux_system.bias_resistance + 1j*omega_b*dfmux_system.bias_inductance) # np.abs(Rb + 1j*omega_b*Lb) 
    
    # impedance of TES leg (assuming on LC resonance)
    Z_tes = np.abs(dfmux_system.operating_resistance + dfmux_system.stray_resistance + 1j*omega_b*dfmux_system.stripline_inductance) 
    # impedance of squid input coil
    Z_sq = np.abs(1j*omega_b*dfmux_system.squid_input_inductance)
    
    # series impedance of wire harness
    Z_wire_harness = np.abs(1j*omega_b * dfmux_system.wire_harness_inductance + dfmux_system.wire_harness_resistance)
    
    # Following the definitions in my drawings in page 120 of reseaerch notebook 05/24/24
    Z1 = np.abs(Z_b + Z_tes )
    Z2 = Z1 * Z_sq / (Z1 + Z_sq)
    
    Z3 = np.abs( dfmux_system.r48 + 1/(1j*omega_b*dfmux_system.parasitic_capacitance) )
    Z4 = Z_wire_harness + Z2
    
    chi_cs = Z3/(Z3+Z4) * Z1/(Z1 + Z_sq)
    csf = 1/chi_cs
    abs_csf = np.abs(csf)

    if debug:    
        if dfmux_system.f[0] > 4.2e6:
            print("\nf: ", dfmux_system.f[0], "  R_ref: ", dfmux_system.r48, " Cpar_to_ground: ", dfmux_system.parasitic_capacitance)
            print(f"Z1: {np.abs(Z1):.2f}, Z2: {np.abs(Z2[0]):.2f} Z3: {np.abs(Z3):.2f} Z4: {np.abs(Z4[0]):.2f}")
            print("Zwh: ", Z_wire_harness)
        print("chi_cs: ", np.abs(chi_cs))
        print("csf: ", abs_csf, "\n")
    
    return abs_csf 
    
    
    
    
    
    
    