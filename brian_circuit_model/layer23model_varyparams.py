# layer23_model_varyparams_extra.py is a model of barrel cortex mouse layer 2/3 that includes
# one excitatory cell type, three inhibitory cell types and external inputs.
# Log-normal synaptic strengths.
# Condunctance-based synapses.
# Repeat input simtulus a few times to get average response probability.
# Average over multiple shuffles of input.
# Vary sparsity of input
# Vary all params by +20%

from brian2 import *
import fractions as frac
import scipy.io
dt = 0.01 * ms
duration = 50 * ms

import time
start_time = time.time()

# SIM PARAMETERS
nloop_fixinput = 10 # Number of loops for each input pattern
nloop_shuffleinput = 10 # Number of input pattern shuffles at fixed sparsity
pinput_vec = linspace(0.1,1,25) # Vector of input sparsities to try
#pinput_vec = linspace(0.1,1,46) # Vector of input sparsities to try
ninput_vec = len(pinput_vec) # Number of test input sparsities
fchange = 0.2 # Fractional change in parameters (can be +ve or -ve)

#####################
# DYNAMICS EQUATIONS
#####################

eqs_E = Equations('''
                    dV/dt = ( Rin_E*( ge*(Erev_e - V) + gi*(Erev_iE - V) ) - (V-vrest) )/taum_E : volt (unless refractory)
                    dge/dt = -ge/tau_Ee : siemens
                    dgi/dt = -gi/tau_Ei : siemens
                    vth : volt
                    vrest : volt
			''')

eqs_I5ht = Equations('''
                    dV/dt = ( Rin_I5ht*( ge*(Erev_e - V) + gi*(Erev_iI5ht - V) ) - (V-vrest))/taum_I5ht : volt (unless refractory)
                    dge/dt = -ge/tau_Ee : siemens
                    dgi/dt = -gi/tau_Ei : siemens
                    vth : volt
                    vrest : volt
			''')
eqs_Ipv = Equations('''
                    dV/dt = ( Rin_Isom*( ge*(Erev_e - V) + gi*(Erev_iIpv - V) ) -(V-vrest))/taum_Ipv : volt (unless refractory)
                    dge/dt = -ge/tau_Ee : siemens
                    dgi/dt = -gi/tau_Ei : siemens
                    vth : volt
                    vrest : volt
			''')
eqs_Isom = Equations('''
                    dV/dt = ( Rin_Ipv*( ge*(Erev_e - V) + gi*(Erev_iIsom - V) ) -(V-vrest))/taum_Isom : volt (unless refractory)
                    dge/dt = -ge/tau_Ee : siemens
                    dgi/dt = -gi/tau_Ei : siemens
                    vth : volt
                    vrest : volt
			''')
    

###############
# LIST PARAMETER MULTIPLIERS IN DICT
###############
pshift = {}
pshift['nochange'] = 0
pshift['N_E'] = 0
pshift['N_I5ht'] = 0
pshift['N_Ipv'] = 0
pshift['N_Isom'] = 0
pshift['N_EL4'] = 0

pshift['p_con_EL4E'] = 0
pshift['p_con_EL4Ipv'] = 0
pshift['p_con_EL4Isom'] = 0
pshift['prel_EL4E'] = 0
pshift['prel_EL4Ipv'] = 0
pshift['prel_EL4Isom'] = 0
pshift['w_EL4E'] = 0
pshift['w_EL4Ipv'] = 0
pshift['w_EL4Isom'] = 0

pshift['Vrest_E'] = 0
pshift['Vrest_I5ht'] = 0
pshift['Vrest_Ipv'] = 0
pshift['Vrest_Isom'] = 0
pshift['Vth_E'] = 0
pshift['Vth_I5ht'] = 0
pshift['Vth_Ipv'] = 0
pshift['Vth_Isom'] = 0
pshift['Rin_E'] = 0
pshift['Rin_I5ht'] = 0
pshift['Rin_Ipv'] = 0
pshift['Rin_Isom'] = 0
pshift['taum_E'] = 0
pshift['taum_I5ht'] = 0
pshift['taum_Ipv'] = 0
pshift['taum_Isom'] = 0
pshift['pcon_EE'] = 0
pshift['pcon_EI5ht'] = 0
pshift['pcon_EIpv'] = 0
pshift['pcon_EIsom'] = 0
pshift['pcon_I5htE'] = 0
pshift['pcon_I5htI5ht'] = 0
pshift['pcon_I5htIpv'] = 0
pshift['pcon_IpvE'] = 0
pshift['pcon_IpvI5ht'] = 0
pshift['pcon_IpvIpv'] = 0
pshift['pcon_IsomE'] = 0
pshift['prel_EE'] = 0
pshift['prel_EI5ht'] = 0
pshift['prel_EIpv'] = 0
pshift['prel_EIsom'] = 0
pshift['prel_I5htE'] = 0
pshift['prel_I5htI5ht'] = 0
pshift['prel_I5htIpv'] = 0
pshift['prel_IpvE'] = 0
pshift['prel_IpvI5ht'] = 0
pshift['prel_IpvIpv'] = 0
pshift['prel_IsomE'] = 0
pshift['w_EE'] = 0
pshift['w_EI5ht'] = 0
pshift['w_EIpv'] = 0
pshift['w_EIsom'] = 0
pshift['w_I5htE'] = 0
pshift['w_I5htI5ht'] = 0
pshift['w_I5htIpv'] = 0
pshift['w_IpvE'] = 0
pshift['w_IpvI5ht'] = 0
pshift['w_IpvIpv'] = 0
pshift['w_IsomE'] = 0

pshift['Erev_e'] = 0
pshift['Erev_iE'] = 0
pshift['Erev_iI5ht'] = 0
pshift['Erev_iIpv'] = 0
pshift['Erev_iIsom'] = 0
pshift['tau_Ee'] = 0
pshift['tau_Ei'] = 0
pshift['tau_I5hte'] = 0
pshift['tau_I5hti'] = 0
pshift['tau_Ipve'] = 0
pshift['tau_Ipvi'] = 0
pshift['tau_Isome'] = 0
pshift['tau_Isomi'] = 0


pshift2 = {}
pshift2['nochange'] = 0
pshift2['N_E'] = 0
pshift2['N_I5ht'] = 0
pshift2['N_Ipv'] = 0
pshift2['N_Isom'] = 0
pshift2['N_EL4'] = 0

pshift2['p_con_EL4E'] = 0
pshift2['p_con_EL4Ipv'] = 0
pshift2['p_con_EL4Isom'] = 0
pshift2['prel_EL4E'] = 0
pshift2['prel_EL4Ipv'] = 0
pshift2['prel_EL4Isom'] = 0
pshift2['w_EL4E'] = 0
pshift2['w_EL4Ipv'] = 0
pshift2['w_EL4Isom'] = 0

pshift2['Vrest_E'] = 0
pshift2['Vrest_I5ht'] = 0
pshift2['Vrest_Ipv'] = 0
pshift2['Vrest_Isom'] = 0
pshift2['Vth_E'] = 0
pshift2['Vth_I5ht'] = 0
pshift2['Vth_Ipv'] = 0
pshift2['Vth_Isom'] = 0
pshift2['Rin_E'] = 0
pshift2['Rin_I5ht'] = 0
pshift2['Rin_Ipv'] = 0
pshift2['Rin_Isom'] = 0
pshift2['taum_E'] = 0
pshift2['taum_I5ht'] = 0
pshift2['taum_Ipv'] = 0
pshift2['taum_Isom'] = 0
pshift2['pcon_EE'] = 0
pshift2['pcon_EI5ht'] = 0
pshift2['pcon_EIpv'] = 0
pshift2['pcon_EIsom'] = 0
pshift2['pcon_I5htE'] = 0
pshift2['pcon_I5htI5ht'] = 0
pshift2['pcon_I5htIpv'] = 0
pshift2['pcon_IpvE'] = 0
pshift2['pcon_IpvI5ht'] = 0
pshift2['pcon_IpvIpv'] = 0
pshift2['pcon_IsomE'] = 0
pshift2['prel_EE'] = 0
pshift2['prel_EI5ht'] = 0
pshift2['prel_EIpv'] = 0
pshift2['prel_EIsom'] = 0
pshift2['prel_I5htE'] = 0
pshift2['prel_I5htI5ht'] = 0
pshift2['prel_I5htIpv'] = 0
pshift2['prel_IpvE'] = 0
pshift2['prel_IpvI5ht'] = 0
pshift2['prel_IpvIpv'] = 0
pshift2['prel_IsomE'] = 0
pshift2['w_EE'] = 0
pshift2['w_EI5ht'] = 0
pshift2['w_EIpv'] = 0
pshift2['w_EIsom'] = 0
pshift2['w_I5htE'] = 0
pshift2['w_I5htI5ht'] = 0
pshift2['w_I5htIpv'] = 0
pshift2['w_IpvE'] = 0
pshift2['w_IpvI5ht'] = 0
pshift2['w_IpvIpv'] = 0
pshift2['w_IsomE'] = 0

pshift2['Erev_e'] = 0
pshift2['Erev_iE'] = 0
pshift2['Erev_iI5ht'] = 0
pshift2['Erev_iIpv'] = 0
pshift2['Erev_iIsom'] = 0
pshift2['tau_Ee'] = 0
pshift2['tau_Ei'] = 0
pshift2['tau_I5hte'] = 0
pshift2['tau_I5hti'] = 0
pshift2['tau_Ipve'] = 0
pshift2['tau_Ipvi'] = 0
pshift2['tau_Isome'] = 0
pshift2['tau_Isomi'] = 0



##############
# LOOP OVER ALL CHOSEN PARAMETERS
##############

allresults={} # Initialise results dictionary

for key in pshift2:
    pshift[key] = fchange
    
    # PARAMETERS
    N_E = 1700 + round(pshift['N_E']*1700)			# Number of excitatory neurons
    N_I5ht = 115 + round(pshift['N_I5ht']*115)			# Number of 5HT3A inhibitory neurons
    N_Ipv = 70 + round(pshift['N_Ipv']*70)			# Number of PV inhibitory neurons
    N_Isom = 45 + round(pshift['N_Isom']*45)				# Number of somatostatin inhibitory neurons
    
    
    ##########
    # L4 INPUT
    ##########
    
    N_EL4 = 1500 + pshift['N_EL4']*1500              # Number of Poisson input neurons from Layer 4
    p_con_EL4E = 0.15 + pshift['p_con_EL4E']*0.15
    p_con_EL4Ipv = 0.15 + pshift['p_con_EL4Ipv']*0.15
    p_con_EL4Isom = 0.15 + pshift['p_con_EL4Isom']*0.15          # Connection probability from L4 to L2/3
    w_EL4E_mean = 0.8 * mV + pshift['w_EL4E']*0.8 * mV
    w_EL4E_median = 0.48 * mV + pshift['w_EL4E']*0.48 * mV #Strength of L4E syn
    w_EL4Ipv_mean = 0.8 * mV + pshift['w_EL4Ipv']*0.8 * mV
    w_EL4Ipv_median = 0.48 * mV + pshift['w_EL4Ipv']*0.48 * mV
    w_EL4Isom_mean = 0.8 * mV + pshift['w_EL4Isom']*0.8 * mV
    w_EL4Isom_median = 0.48 * mV + pshift['w_EL4Isom']*0.48 * mV
    
    prel_EL4E = 0.25 + pshift['prel_EL4E']*0.25
    #prel_EL4I5ht = 0.25
    prel_EL4Ipv = 0.25 + pshift['prel_EL4Ipv']*0.25
    prel_EL4Isom = 0.25 + pshift['prel_EL4Isom']*0.25
    
    ###############
    # NEURONS
    ###############
    
    Vrest_E = -68 * mV + pshift['Vrest_E']*abs(-68+38) * mV 		# Resting potential of E neurons
    Vrest_I5ht = -62 * mV + pshift['Vrest_I5ht']*abs(-62+36) * mV	# Resting potential of 5HT3A I neurons
    Vrest_Ipv = -68	* mV + pshift['Vrest_Ipv']*abs(-68+37.4) * mV	# Resting potential of PV I neurons
    Vrest_Isom = -57 * mV + pshift['Vrest_Isom']*abs(-57+40) * mV	# Resting potential of SOM I neurons
    Vrest_sd = 0 * mV   # Resting potential s.d. all cell types
    
    Vref_i = - 55 * mV # potential to calculate IPSPs at
    
    Vth_E = -38	* mV + pshift['Vrest_E']*abs(-68+38) * mV 		# Spike threshold
    Vth_I5ht = -36 * mV + pshift['Vrest_I5ht']*abs(-62+36) * mV
    Vth_Ipv = -37.4 * mV + pshift['Vrest_Ipv']*abs(-68+37.4) * mV
    Vth_Isom = -40 * mV + pshift['Vrest_Isom']*abs(-57+40) * mV
    Vth_sd = 0 * mV   # Spike threshold s.d. all cell types
    
    Rin_E = 160	 * Mohm + pshift['Rin_E']*160 * Mohm		# Input resistance
    Rin_I5ht = 200 * Mohm + pshift['Rin_I5ht']*200 * Mohm
    Rin_Ipv = 100 * Mohm + pshift['Rin_Ipv']*100 * Mohm
    Rin_Isom = 250 * Mohm + pshift['Rin_Isom']*250 * Mohm
    
    taum_E = 28 * ms + pshift['taum_E']*28 * ms		# Membrane time constant
    taum_I5ht = 21 * ms + pshift['taum_I5ht']*21 * ms
    taum_Ipv = 10 * ms + pshift['taum_Ipv']*10 * ms
    taum_Isom = 30 * ms + pshift['taum_Isom']*30 * ms
    
    ref_E = frac.Fraction(1,18) * second	# Refractory period (set to inverse of max FR)
    ref_I5ht = frac.Fraction(1,47) * second
    ref_Ipv = frac.Fraction(1,185) * second
    ref_Isom = frac.Fraction(1,50) * second
    
    
    ##############
    # SYNPASES
    ##############
    
    tau_Ee = 2 * ms + pshift['tau_Ee']* 2 * ms			# Synaptic time constant
    tau_Ei = 40 * ms + pshift['tau_Ei']* 40 * ms
    tau_I5hte = 2 * ms + pshift['tau_I5hte']* 2 * ms
    tau_I5hti = 40 * ms + pshift['tau_I5hti']* 40 * ms
    tau_Ipve = 2 * ms + pshift['tau_Ipve']* 2 * ms
    tau_Ipvi = 16 * ms + pshift['tau_Ipvi']* 16 * ms
    tau_Isome = 2 * ms + pshift['tau_Isome']* 2 * ms
    tau_Isomi = 40 * ms + pshift['tau_Isomi']* 40 * ms
    
    Erev_e = 0 * mV  + pshift['Erev_e']*abs(-68 + 0) * mV  # Rev pot of synapses
    Erev_iE = Vrest_E + pshift['Erev_iE']*abs(-68+38) * mV
    Erev_iI5ht = Vrest_I5ht + pshift['Erev_iI5ht']*abs(-62+36) * mV
    Erev_iIpv = Vrest_Ipv + pshift['Erev_iIpv']*abs(-68+37.4) * mV
    Erev_iIsom = Vrest_Isom + pshift['Erev_iIsom']*abs(-57+40) * mV
    
    pcon_EE = 0.17 + pshift['pcon_EE']* 0.17 			# Connection probability
    pcon_EI5ht = 0.24 + pshift['pcon_EI5ht']* 0.24 
    pcon_EIpv = 0.575 + pshift['pcon_EIpv']* 0.575 
    pcon_EIsom = 0.5 + pshift['pcon_EIsom']* 0.5 
    pcon_I5htE = 0.465 + pshift['pcon_I5htE']* 0.465 
    pcon_I5htI5ht = 0.38 + pshift['pcon_I5htI5ht']* 0.38 
    pcon_I5htIpv = 0.38 + pshift['pcon_I5htIpv']* 0.38 
    pcon_I5htIsom = 0
    pcon_IpvE = 0.6 + pshift['pcon_IpvE']* 0.6 
    pcon_IpvI5ht = 0.24 + pshift['pcon_IpvI5ht']* 0.24 
    pcon_IpvIpv = 0.55 + pshift['pcon_IpvIpv']* 0.55 
    pcon_IpvIsom = 0
    pcon_IsomE = 0.5 + pshift['pcon_IsomE']* 0.5 
    pcon_IsomI5ht = 0
    pcon_IsomIpv = 0
    pcon_IsomIsom = 0
    
    prel_EE = 0.25 + pshift['prel_EE']* 0.25			# Connection probability
    prel_EI5ht = 0.25 + pshift['prel_EI5ht']* 0.25
    prel_EIpv = 0.25 + pshift['prel_EIpv']* 0.25
    prel_EIsom = 0.25 + pshift['prel_EIsom']* 0.25
    prel_I5htE = 0.25 + pshift['prel_I5htE']* 0.25
    prel_I5htI5ht = 0.25 + pshift['prel_I5htI5ht']* 0.25
    prel_I5htIpv = 0.25 + pshift['prel_I5htIpv']* 0.25
    #prel_I5htIsom = 0
    prel_IpvE = 0.25 + pshift['prel_IpvE']* 0.25
    prel_IpvI5ht = 0.25 + pshift['prel_IpvI5ht']* 0.25
    prel_IpvIpv = 0.25 + pshift['prel_IpvIpv']* 0.25
    #prel_IpvIsom = 0
    prel_IsomE = 0.25 + pshift['prel_IsomE']* 0.25
    #prel_IsomI5ht = 0
    #prel_IsomIpv = 0
    #prel_IsomIsom = 0
    
    w_EE_mean = 0.37 * mV + pshift['w_EE']*0.37*mV		# EPSP/IPSP amplitude
    w_EI5ht_mean = 0.39 * mV + pshift['w_EI5ht']*0.39*mV
    w_EIpv_mean = 0.82 * mV + pshift['w_EIpv']*0.92*mV
    w_EIsom_mean = 0.5 * mV + pshift['w_EIsom']*0.5*mV
    w_I5htE_mean = 0.49  * mV + pshift['w_I5htE']*0.49*mV
    w_I5htI5ht_mean = 0.49 * mV + pshift['w_I5htI5ht']*0.37*mV
    w_I5htIpv_mean = 0.37 * mV + pshift['w_I5htIpv']*0.39*mV
    #w_I5htIsom_mean = 0 * mV
    w_IpvE_mean = 0.52 * mV + pshift['w_IpvE']*0.52*mV
    w_IpvI5ht_mean = 0.83 * mV + pshift['w_IpvI5ht']*0.83*mV
    w_IpvIpv_mean = 0.56 * mV + pshift['w_IpvIpv']*0.56*mV
    #w_IpvIsom_mean = 0 * mV
    w_IsomE_mean = 0.5 * mV + pshift['w_IsomE']*0.5*mV
    #w_IsomI5ht_mean = 0 * mV
    #w_IsomIpv_mean = 0 * mV
    #w_IsomIsom_mean = 0 * mV
    
    w_EE_median = 0.2 * mV	 + pshift['w_EE']*0.2*mV	# EPSP/IPSP median amplitude
    w_EI5ht_median = 0.19 * mV + pshift['w_EI5ht']*0.19*mV
    w_EIpv_median = 0.68 * mV + pshift['w_EIpv']*0.68*mV
    w_EIsom_median = 0.4 * mV + pshift['w_EIsom']*0.4*mV
    w_I5htE_median = 0.3  * mV + pshift['w_I5htE']*0.3*mV
    w_I5htI5ht_median = 0.15 * mV + pshift['w_I5htI5ht']*0.15*mV
    w_I5htIpv_median = 0.23 * mV + pshift['w_I5htIpv']*0.23*mV
    #w_I5htIsom_median = 0 * mV
    w_IpvE_median = 0.29 * mV + pshift['w_IpvE']*0.29*mV
    w_IpvI5ht_median = 0.6 * mV + pshift['w_IpvI5ht']*0.6*mV
    w_IpvIpv_median = 0.44 * mV + pshift['w_IpvIpv']*0.44*mV
    #w_IpvIsom_median = 0 * mV
    w_IsomE_median = 0.4 * mV + pshift['w_IsomE']*0.4*mV
    #w_IsomI5ht_median = 0 * mV
    #w_IsomIpv_median = 0 * mV
    #w_IsomIsom_median = 0 * mV
    
    w_max = 8 # Maximum synaptic weight in mV
    
    
    # CALCULATE SYNAPTIC CONDUCTANCES AND TIME CONSTANTS FROM EPSP CONSTRAINTS
    def tpeakfunc(tau_m,tau_s): # Time of EPSP peak
        return log(tau_m/tau_s)/((1/tau_s)-(1/tau_m))
    
    def gestfunc(tpeak, Vpeak, tau_m, tau_s, R, Erev, Vrest):
        return Vpeak*(tau_m - tau_s) / ( tau_s*R*(exp(-tpeak/tau_m) - exp(-tpeak/tau_s) )*(Erev - Vrest) )
    
    gsyn_EE_mean = gestfunc(tpeakfunc(taum_E,tau_Ee), w_EE_mean, taum_E, tau_Ee, Rin_E, Erev_e, Vrest_E)         # Synaptic condunctance means
    gsyn_EI5ht_mean = gestfunc(tpeakfunc(taum_I5ht,tau_I5hte), w_EI5ht_mean, taum_I5ht, tau_I5hte, Rin_I5ht, Erev_e, Vrest_I5ht)
    gsyn_EIpv_mean = gestfunc(tpeakfunc(taum_Ipv,tau_Ipve), w_EIpv_mean, taum_Ipv, tau_Ipve, Rin_Ipv, Erev_e, Vrest_Ipv)
    gsyn_EIsom_mean = gestfunc(tpeakfunc(taum_Isom,tau_Isome), w_EIsom_mean, taum_Isom, tau_Isome, Rin_Isom, Erev_e, Vrest_Isom)
    gsyn_I5htE_mean = gestfunc(tpeakfunc(taum_E,tau_Ei), -w_I5htE_mean, taum_E, tau_Ei, Rin_E, Erev_iE, Vref_i)
    gsyn_I5htI5ht_mean = gestfunc(tpeakfunc(taum_I5ht,tau_I5hti), -w_I5htI5ht_mean, taum_I5ht, tau_I5hti, Rin_I5ht, Erev_iI5ht, Vref_i)
    gsyn_I5htIpv_mean = gestfunc(tpeakfunc(taum_Ipv,tau_Ipvi), -w_I5htIpv_mean, taum_Ipv, tau_Ipvi, Rin_Ipv, Erev_iIpv, Vref_i)
    #gsyn_I5htIsom_mean = gestfunc(tpeakfunc(taum_Isom,tau_Isomi), -w_I5htIsom_mean, taum_Isom, tau_Isomi, Rin_Isom, Erev_iIsom, Vref_i)
    gsyn_IpvE_mean = gestfunc(tpeakfunc(taum_E,tau_Ei), -w_IpvE_mean, taum_E, tau_Ei, Rin_E, Erev_iE, Vref_i)
    gsyn_IpvI5ht_mean = gestfunc(tpeakfunc(taum_I5ht,tau_I5hti), -w_IpvI5ht_mean, taum_I5ht, tau_I5hti, Rin_I5ht, Erev_iI5ht, Vref_i)
    gsyn_IpvIpv_mean = gestfunc(tpeakfunc(taum_Ipv,tau_Ipvi), -w_IpvIpv_mean, taum_Ipv, tau_Ipvi, Rin_Ipv, Erev_iIpv, Vref_i)
    #gsyn_IpvIsom_mean = gestfunc(tpeakfunc(taum_Isom,tau_Isomi), -w_IpvIsom_mean, taum_Isom, tau_Isomi, Rin_Isom, Erev_iIsom, Vref_i)
    gsyn_IsomE_mean = gestfunc(tpeakfunc(taum_E,tau_Ei), -w_IsomE_mean, taum_E, tau_Ei, Rin_E, Erev_iE, Vref_i)
    #gsyn_IsomI5ht_mean = gestfunc(tpeakfunc(taum_I5ht,tau_I5hti),-w_IsomI5ht_mean, taum_I5ht, tau_I5hti, Rin_I5ht, Erev_iI5ht, Vref_i)
    #gsyn_IsomIpv_mean = gestfunc(tpeakfunc(taum_Ipv,tau_Ipvi), -w_IsomIpv_mean, taum_Ipv, tau_Ipvi, Rin_Ipv, Erev_iIpv, Vref_i)
    #gsyn_IsomIsom_mean = gestfunc(tpeakfunc(taum_Isom,tau_Isomi), -w_IsomIsom_mean, taum_Isom, tau_Isomi, Rin_Isom, Erev_iIsom, Vref_i)
    
    gsyn_EE_median = gestfunc(tpeakfunc(taum_E,tau_Ee), w_EE_median, taum_E, tau_Ee, Rin_E, Erev_e, Vrest_E)         # Synaptic condunctance medians
    gsyn_EI5ht_median = gestfunc(tpeakfunc(taum_I5ht,tau_I5hte), w_EI5ht_median, taum_I5ht, tau_I5hte, Rin_I5ht, Erev_e, Vrest_I5ht)
    gsyn_EIpv_median = gestfunc(tpeakfunc(taum_Ipv,tau_Ipve), w_EIpv_median, taum_Ipv, tau_Ipve, Rin_Ipv, Erev_e, Vrest_Ipv)
    gsyn_EIsom_median = gestfunc(tpeakfunc(taum_Isom,tau_Isome), w_EIsom_median, taum_Isom, tau_Isome, Rin_Isom, Erev_e, Vrest_Isom)
    gsyn_I5htE_median = gestfunc(tpeakfunc(taum_E,tau_Ei), -w_I5htE_median, taum_E, tau_Ei, Rin_E, Erev_iE, Vref_i)
    gsyn_I5htI5ht_median = gestfunc(tpeakfunc(taum_I5ht,tau_I5hti), -w_I5htI5ht_median, taum_I5ht, tau_I5hti, Rin_I5ht, Erev_iI5ht, Vref_i)
    gsyn_I5htIpv_median = gestfunc(tpeakfunc(taum_Ipv,tau_Ipvi), -w_I5htIpv_median, taum_Ipv, tau_Ipvi, Rin_Ipv, Erev_iIpv, Vref_i)
    #gsyn_I5htIsom_median = gestfunc(tpeakfunc(taum_Isom,tau_Isomi), -w_I5htIsom_median, taum_Isom, tau_Isomi, Rin_Isom, Erev_iIsom, Vref_i)
    gsyn_IpvE_median = gestfunc(tpeakfunc(taum_E,tau_Ei), -w_IpvE_median, taum_E, tau_Ei, Rin_E, Erev_iE, Vref_i)
    gsyn_IpvI5ht_median = gestfunc(tpeakfunc(taum_I5ht,tau_I5hti), -w_IpvI5ht_median, taum_I5ht, tau_I5hti, Rin_I5ht, Erev_iI5ht, Vref_i)
    gsyn_IpvIpv_median = gestfunc(tpeakfunc(taum_Ipv,tau_Ipvi), -w_IpvIpv_median, taum_Ipv, tau_Ipvi, Rin_Ipv, Erev_iIpv, Vref_i)
    #gsyn_IpvIsom_median = gestfunc(tpeakfunc(taum_Isom,tau_Isomi), -w_IpvIsom_median, taum_Isom, tau_Isomi, Rin_Isom, Erev_iIsom, Vref_i)
    gsyn_IsomE_median = gestfunc(tpeakfunc(taum_E,tau_Ei), -w_IsomE_median, taum_E, tau_Ei, Rin_E, Erev_iE, Vref_i)
    #gsyn_IsomI5ht_median = gestfunc(tpeakfunc(taum_I5ht,tau_I5hti), -w_IsomI5ht_median, taum_I5ht, tau_I5hti, Rin_I5ht, Erev_iI5ht, Vref_i)
    #gsyn_IsomIpv_median = gestfunc(tpeakfunc(taum_Ipv,tau_Ipvi), -w_IsomIpv_median, taum_Ipv, tau_Ipvi, Rin_Ipv, Erev_iIpv, Vref_i)
    #gsyn_IsomIsom_median = gestfunc(tpeakfunc(taum_Isom,tau_Isomi), -w_IsomIsom_median, taum_Isom, tau_Isomi, Rin_Isom, Erev_iIsom, Vref_i)
    
    gmax = 20 # maximum synaptic conductance in nS
    
    ####################		
    # CREATE NEURONS
    ##################
    
    G_E = NeuronGroup(N_E, model=eqs_E, threshold='V > vth', reset='V = vrest', refractory=ref_E, method = 'euler')
    G_I5ht = NeuronGroup(N_I5ht, model=eqs_I5ht, threshold='V > vth', reset='V = vrest', refractory=ref_I5ht, method = 'euler' )
    G_Ipv = NeuronGroup(N_Ipv, model=eqs_Ipv, threshold = 'V > vth', reset='V = vrest', refractory=ref_Ipv, method = 'euler')
    G_Isom = NeuronGroup(N_Isom, model=eqs_Isom, threshold = 'V > vth', reset='V = vrest', refractory=ref_Isom, method = 'euler')
    
    # SYNAPSES (4 cell types x 4 = 16 types)
    S_EE = Synapses(G_E,G_E,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EE))')
    S_EI5ht = Synapses(G_E,G_I5ht,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EI5ht))')
    S_EIpv = Synapses(G_E,G_Ipv,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EIpv))')
    S_EIsom = Synapses(G_E,G_Isom,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EIsom))')
    
    S_I5htE = Synapses(G_I5ht,G_E,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_I5htE))')
    S_I5htI5ht = Synapses(G_I5ht,G_I5ht,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_I5htI5ht))')
    S_I5htIpv = Synapses(G_I5ht,G_Ipv,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_I5htIpv))')
    #S_I5htIsom = Synapses(G_I5ht,G_Isom,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_I5htIsom))')
    
    S_IpvE = Synapses(G_Ipv,G_E,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IpvE))')
    S_IpvI5ht = Synapses(G_Ipv,G_I5ht,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IpvI5ht))')
    S_IpvIpv = Synapses(G_Ipv,G_Ipv,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IpvIpv))')
    #S_IpvIsom = Synapses(G_Ipv,G_Isom,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IpvIsom))')
    
    S_IsomE = Synapses(G_Isom,G_E,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IsomE))')
    #S_IsomI5ht = Synapses(G_Isom,G_I5ht,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IsomI5ht))')
    #S_IsomIpv = Synapses(G_Isom,G_Ipv,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IsomIpv))')
    #S_IsomIsom = Synapses(G_Isom,G_Isom,'gi_syn : siemens',on_pre='gi_post += (gi_syn*(rand()<prel_IsomIsom))')
    
    # CONNECT THEM UP
    S_EE.connect('i!=j', p=pcon_EE)
    S_EI5ht.connect(True, p=pcon_EI5ht)
    S_EIpv.connect(True, p=pcon_EIpv)
    S_EIsom.connect(True, p=pcon_EIsom)
    
    S_I5htE.connect(True, p=pcon_I5htE)
    S_I5htI5ht.connect('i!=j', p=pcon_I5htI5ht)
    S_I5htIpv.connect(True, p=pcon_I5htIpv)
    #S_I5htIsom.connect(True, p=pcon_I5htIsom)
    
    S_IpvE.connect(True, p=pcon_IpvE)
    S_IpvI5ht.connect(True, p=pcon_IpvI5ht)
    S_IpvIpv.connect(True, p=pcon_IpvIpv)
    #S_IpvIsom.connect('i!=j', p=pcon_IsomIsom)
    
    S_IsomE.connect(True, p=pcon_IsomE)
    #S_IsomI5ht.connect(True, p=pcon_IsomI5ht)
    #S_IsomIpv.connect(True, p=pcon_IsomIpv)
    #S_IsomIsom.connect('i!=j', p=pcon_IsomIsom)
    
    
    # SET WEIGHTS FROM LOG NORMAL DISTRIBUTION
    S_EE.ge_syn = clip( exp(log(gsyn_EE_median/nS) + sqrt(2*log(gsyn_EE_mean/gsyn_EE_median))*randn(size(S_EE.ge_syn)) ) , 0, gmax ) * nS
    S_EI5ht.ge_syn = clip( exp(log(gsyn_EI5ht_median/nS) + sqrt(2*log(w_EI5ht_mean/w_EI5ht_median))*randn(size(S_EI5ht.ge_syn)) ) , 0, gmax ) * nS
    S_EIpv.ge_syn = clip( exp(log(gsyn_EIpv_median/nS) + sqrt(2*log(w_EIpv_mean/w_EIpv_median))*randn(size(S_EIpv.ge_syn)) ) , 0, gmax ) * nS
    S_EIsom.ge_syn = clip( exp(log(gsyn_EIsom_median/nS) + sqrt(2*log(w_EIsom_mean/w_EIsom_median))*randn(size(S_EIsom.ge_syn)) ) , 0, gmax ) * nS
    
    S_I5htE.gi_syn = clip( exp(log(gsyn_I5htE_median/nS) + sqrt(2*log(gsyn_I5htE_mean/gsyn_I5htE_median))*randn(size(S_I5htE.gi_syn)) ) , 0, gmax ) * nS
    S_I5htI5ht.gi_syn = clip( exp(log(gsyn_I5htI5ht_median/nS) + sqrt(2*log(gsyn_I5htI5ht_mean/gsyn_I5htI5ht_median))*randn(size(S_I5htI5ht.gi_syn)) ) , 0, gmax ) * nS
    S_I5htIpv.gi_syn = clip( exp(log(gsyn_I5htIpv_median/nS) + sqrt(2*log(gsyn_I5htIpv_mean/gsyn_I5htIpv_median))*randn(size(S_I5htIpv.gi_syn)) ) , 0, gmax ) * nS
    #S_I5htIsom.gi_syn = clip( exp(log(gsyn_I5htIsom_median/nS) + sqrt(2*log(gsyn_I5htIsom_mean/gsyn_I5htIsom_median))*randn(size(S_I5htIsom.gi_syn)) ) , 0, gmax ) * nS
    
    S_IpvE.gi_syn = clip( exp(log(gsyn_IpvE_median/nS) + sqrt(2*log(gsyn_IpvE_mean/gsyn_IpvE_median))*randn(size(S_IpvE.gi_syn)) ) , 0, gmax ) * nS
    S_IpvI5ht.gi_syn = clip( exp(log(gsyn_IpvI5ht_median/nS) + sqrt(2*log(gsyn_IpvI5ht_mean/gsyn_IpvI5ht_median))*randn(size(S_IpvI5ht.gi_syn)) ) , 0, gmax ) * nS
    S_IpvIpv.gi_syn = clip( exp(log(gsyn_IpvIpv_median/nS) + sqrt(2*log(gsyn_IpvIpv_mean/gsyn_IpvIpv_median))*randn(size(S_IpvIpv.gi_syn)) ) , 0, gmax ) * nS
    #S_IpvIsom.gi_syn = clip( exp(log(gsyn_IpvIsom_median/nS) + sqrt(2*log(gsyn_IpvIsom_mean/gsyn_IpvIsom_median))*randn(size(S_IpvIsom.gi_syn)) ) , 0, gmax ) * nS
    
    S_IsomE.gi_syn = clip( exp(log(gsyn_IsomE_median/nS) + sqrt(2*log(gsyn_IsomE_mean/gsyn_IsomE_median))*randn(size(S_IsomE.gi_syn)) ) , 0, gmax ) * nS
    #S_IsomI5ht.gi_syn = clip( exp(log(gsyn_IsomI5ht_median/nS) + sqrt(2*log(gsyn_IsomI5ht_mean/gsyn_IsomI5ht_median))*randn(size(S_IsomI5ht.gi_syn)) ) , 0, gmax ) * nS
    #S_IsomIpv.gi_syn = clip( exp(log(gsyn_IsomIpv_median/nS) + sqrt(2*log(gsyn_IsomIpv_mean/gsyn_IsomIpv_median))*randn(size(S_IsomIpv.gi_syn)) ) , 0, gmax ) * nS
    #S_IsomIsom.gi_syn = clip( exp(log(gsyn_IsomIsom_median/nS) + sqrt(2*log(gsyn_IsomIsom_mean/gsyn_IsomIsom_median))*randn(size(S_IsomIsom.gi_syn)) ) , 0, gmax ) * nS
    
    
    # SET THRESHOLDS
    G_E.vth = Vth_E + Vth_sd*randn(size(G_E.vth))
    G_I5ht.vth = Vth_I5ht + Vth_sd*randn(size(G_I5ht.vth))
    G_Ipv.vth = Vth_Ipv + Vth_sd*randn(size(G_Ipv.vth))
    G_Isom.vth = Vth_Isom + Vth_sd*randn(size(G_Isom.vth))
    
    # SET RESTINF POTENTIALS
    G_E.vrest = Vrest_E + Vrest_sd*randn(size(G_E.vth))
    G_I5ht.vrest = Vrest_I5ht + Vrest_sd*randn(size(G_I5ht.vth))
    G_Ipv.vrest = Vrest_Ipv + Vrest_sd*randn(size(G_Ipv.vth))
    G_Isom.vrest = Vrest_Isom + Vrest_sd*randn(size(G_Isom.vth))
    
    
    # CREATE INPUT
    cellvec = numpy.linspace(0,N_EL4-1,N_EL4)
    cellvec = cellvec.astype(int)
    idx = numpy.arange(cellvec.size)
    numpy.random.shuffle(idx)
    indices = sort(idx[:2])
    times = array(10*ones(size(indices)) + randn(2))*ms
    G_INPUT = SpikeGeneratorGroup(N_EL4, indices, times, name = 'Inputgroup')
    
    gsyn_EL4E_mean = gestfunc(tpeakfunc(taum_E,tau_Ee), w_EL4E_mean, taum_E, tau_Ee, Rin_E, Erev_e, Vrest_E)
    gsyn_EL4E_median = gestfunc(tpeakfunc(taum_E,tau_Ee), w_EL4E_median, taum_E, tau_Ee, Rin_E, Erev_e, Vrest_E)
    # gsyn_EL4I5ht_mean = gestfunc(tpeakfunc(taum_I5ht,tau_I5hte), w_EL4E_mean, taum_I5ht, tau_I5hte, Rin_I5ht, Erev_e, Vrest_I5ht)
    # gsyn_EL4I5ht_median = gestfunc(tpeakfunc(taum_I5ht,tau_I5hte), w_EL4E_median, taum_I5ht, tau_I5hte, Rin_I5ht, Erev_e, Vrest_I5ht)
    gsyn_EL4Ipv_mean = gestfunc(tpeakfunc(taum_Ipv,tau_Ipve), w_EL4Ipv_mean, taum_Ipv, tau_Ipve, Rin_Ipv, Erev_e, Vrest_Ipv)
    gsyn_EL4Ipv_median = gestfunc(tpeakfunc(taum_Ipv,tau_Ipve), w_EL4Ipv_median, taum_Ipv, tau_Ipve, Rin_Ipv, Erev_e, Vrest_Ipv)
    gsyn_EL4Isom_mean = gestfunc(tpeakfunc(taum_Isom,tau_Isome), w_EL4Isom_mean, taum_Isom, tau_Isome, Rin_Isom, Erev_e, Vrest_Isom)
    gsyn_EL4Isom_median = gestfunc(tpeakfunc(taum_Isom,tau_Isome), w_EL4Isom_median, taum_Isom, tau_Isome, Rin_Isom, Erev_e, Vrest_Isom)
    
    S_INPUTE = Synapses(G_INPUT,G_E,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EL4E))')
    S_INPUTE.connect(True, p=p_con_EL4E)
    S_INPUTE.ge_syn = exp(log(gsyn_EL4E_median/nS) + sqrt(2*log(gsyn_EL4E_mean/gsyn_EL4E_median))*randn(size(S_INPUTE.ge_syn)) ) * nS
    #S_INPUTI5ht = Synapses(G_INPUT,G_I5ht,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EL4I5ht))')
    #S_INPUTI5ht.connect(True, p=p_con_EL4I5ht)
    #S_INPUTI5ht.ge_syn = exp(log(gsyn_EL4I5ht_median/nS) + sqrt(2*log(gsyn_EL4I5ht_mean/gsyn_EL4I5ht_median))*randn(size(S_INPUTI5ht.ge_syn)) ) * nS
    S_INPUTIpv = Synapses(G_INPUT,G_Ipv,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EL4Ipv))')
    S_INPUTIpv.connect(True, p=p_con_EL4Ipv)
    S_INPUTIpv.ge_syn = exp(log(gsyn_EL4Ipv_median/nS) + sqrt(2*log(gsyn_EL4Ipv_mean/gsyn_EL4Ipv_median))*randn(size(S_INPUTIpv.ge_syn)) ) * nS
    S_INPUTIsom = Synapses(G_INPUT,G_Isom,'ge_syn : siemens',on_pre='ge_post += (ge_syn*(rand()<prel_EL4Isom))')
    S_INPUTIsom.connect(True, p=p_con_EL4Isom)
    S_INPUTIsom.ge_syn = exp(log(gsyn_EL4Isom_median/nS) + sqrt(2*log(gsyn_EL4Isom_mean/gsyn_EL4Isom_median))*randn(size(S_INPUTIsom.ge_syn)) ) * nS
    
    
    
    # RECORD
    M_E = SpikeMonitor(G_E)
    M_I5ht = SpikeMonitor(G_I5ht)
    M_Ipv = SpikeMonitor(G_Ipv)
    M_Isom = SpikeMonitor(G_Isom)
    
    # INITIALIZE STATE VARIABLES
    G_E.V = Vrest_E # + (Vth_E-Vrest_E) * rand(len(G_E))
    G_I5ht.V = Vrest_I5ht
    G_Ipv.V = Vrest_Ipv
    G_Isom.V = Vrest_Isom
    
    # COLLECT NETWORK
    net = Network(collect())
    net.store()
    
    
    # INITIALISE SOME STUFF
    pspike_E = zeros([ninput_vec,N_E])
    pspike_I5ht = zeros([ninput_vec,N_I5ht])
    pspike_Ipv = zeros([ninput_vec,N_Ipv])
    pspike_Isom = zeros([ninput_vec,N_Isom])
    
    # RUN
    for k in range (ninput_vec):
        NL4_ON = int(pinput_vec[k]*N_EL4)
        
        nspikeloop_E = zeros([1,N_E])
        nspikeloop_I5ht = zeros([1,N_I5ht])
        nspikeloop_Ipv = zeros([1,N_Ipv])
        nspikeloop_Isom = zeros([1,N_Isom])
        
        print "Fraction L4 active = %s" %(pinput_vec[k]),"\n"
    
        for j in range(nloop_shuffleinput): # Loop over input shuffles
            
            indices = sort(idx[:NL4_ON])
            times = array(10*ones(size(indices)) + randn(NL4_ON))*ms
            G_INPUT.set_spikes(indices,times)
            
            net.store()
                
            for i in range(nloop_fixinput): # loop over fixed input probability
                net.run(duration)
                
                nspikeloop_E = nspikeloop_E + clip(M_E.count,0,1)
                nspikeloop_I5ht = nspikeloop_I5ht + clip(M_I5ht.count,0,1)
                nspikeloop_Ipv = nspikeloop_Ipv + clip(M_Ipv.count,0,1)
                nspikeloop_Isom = nspikeloop_Isom + clip(M_Isom.count,0,1)
                
                net.restore()
                net.store()
        
        pspike_E[k,:] = nspikeloop_E/(nloop_fixinput*nloop_shuffleinput)
        pspike_I5ht[k,:] = nspikeloop_I5ht/(nloop_fixinput*nloop_shuffleinput)
        pspike_Ipv[k,:] = nspikeloop_Ipv/(nloop_fixinput*nloop_shuffleinput)
        pspike_Isom[k,:] = nspikeloop_Isom/(nloop_fixinput*nloop_shuffleinput)
    
        runresults = dict(pspike_E=pspike_E, pspike_I5ht=pspike_I5ht, pspike_Ipv=pspike_Ipv, pspike_Isom=pspike_Isom, pinput_vec=pinput_vec)
        
    allresults[key] = runresults
    scipy.io.savemat('layer23_pspike_varyparam_20pc_10pats_pos_extra2.mat', allresults)
    pshift = dict.fromkeys(pshift, 0) # Reset shift parameters


# POST ANALYSIS

##########
# PLOT
##########

#subplot(411)
#plot(pinput_vec,pspike_E[:,:10],'b')
#ylabel('Exc')
#ylim([0,1])
#subplot(412)
#plot(pinput_vec,pspike_I5ht[:,:10],'g')
#ylabel('5ht')
#ylim([0,1])
#subplot(413)
#plot(pinput_vec,pspike_Ipv[:,:10],'r')
#ylabel('PV')
#ylim([0,1])
#subplot(414)
#plot(pinput_vec,pspike_Isom[:,:10],'c')
#ylabel('SOM')
#ylim([0,1])
#xlabel('Fraction of L4 active')


print("--- %s minutes ---" % ((time.time() - start_time)/60))