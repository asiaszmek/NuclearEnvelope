c_unit = 6.0221409e5
ip3rtau = 2000 #wagner
ip3degTau = 1000
gip3r = 120400

gleak = 15.905e-3
ip3_init = 25 # 50 nM? signalling pathways model

calr_tot = 86 #  calreticulin
calr_bound = 7.11945
kf_calr = 0.1 # 1/mM/ms
kb_calr = 0.2 #  1/ms


ca_init_cyt = 100e-6
ca_init_ER = .180512
caDiff = 0.174
ip3Diff = 0.1
camDiff = 0.066
calbDiff = 0.066

AtoN_ratio = 2.1
gAMPA = 25e-3
gNMDA = gAMPA/AtoN_ratio


Kip3 = 0.15e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
Kact = 0.8e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
k_inh = 2e-3 #Wagner 2004 oocytes and later Ausra's astrocyte model
Km_Serca = 0.13e-3 # Michaelis constant for SERCA pump
kcat_Serca = 7.5e-3 #ms
gSerca = 1#signaling pathways model
ip3r_gate_state = 0.8
calbindin_tot = 0.150 #(150 uM)
calmodulin_tot = 0.03 #(30 uM)
fixed_buffer_tot = 2 #(2 mM)


kf_calbindin = 2.8e-2
kb_calbindin = 0.0196

kf_camn = 2*7.7e2
kb_camn = 1.6e2
kf_camc = 2*8.4e1
kb_camc = 2.6
kf_fixed_b = 400
kb_fixed_b = 40
fixed_buffer_ca = kf_fixed_b*ca_init*fixed_buffer_tot/kb_fixed_b

camn = kf_camn*ca_init*calmodulin_tot/kb_camn
camc = kf_camc*ca_init*calmodulin_tot/kb_camc
calbca = kf_calbindin*ca_init*calbindin_tot/kb_calbindin

g_leak_ECS = 1
g_leak_spine = 1

kf_pmca = 5000
kb_pmca = 0.7
kcat_pmca = 1
Km_pmca = kb_pmca/kf_pmca
gpmca = 100e-5*ca_factor # {"apical_dendrite[10]":0.1e-5*ca_factor}
gpmca_spine = 50e-5*ca_factor # {"apical_dendrite[10]": 0.1e-5*ca_factor}

ncx_pow = 1
kf_ncx = 168
kb_ncx = 0.112
kcat_ncx = 5
Km_ncx = kb_ncx/kf_ncx
#  this dynamics is more similar to quasi-steady state approx
gncx = 1e-6*ca_factor#0.6e-8*ca_factor#1.035e-5*ca_factor#{"apical_dendrite[10]": 1.035e-5*ca_factor}#4.6875e-3}
gncx_spine = gncx/7 #{}
# for key in gncx:
#     gncx_spine[key] = gncx[key]/7
# 1/7 of dends  https://doi.org/10.1073/pnas.0605412104 
Ca_Ext = 2
n_seg = 1


#Ca indicators

#For Sabatinis life cycle of a Ca2+ ion Ca in the spine tuning
tot_magnesium_green_BS = 0.100
magnesium_green_bound = 0.0014972
kf_magnesium_green = 9e1#1/mM/ms from doi: 10.1016/j.ymeth.2008.09.025
kb_magnesium_green = 0.6 #1/ms, Kd in vitro 6uM, doi: 10.1016/j.ymeth.2008.09.025 and Sabatini
mggreenDiff = 15e-3   #15 um^2/s doi:10.1016/S0006-3495(96)79633-9 

#For Marsden et al. 10.1073/pnas.1010346107
#Fluo-3, constants Pflügers Arch – Eur J Physiol (1997) 434:615–631
tot_fluo3 = 1e-3 #  1 uM
kf_fluo3 = 0.0236e3 # 0.236 1/uM/ms
kb_fluo3 = 0.0236   # 1/ms
fluo3Diff = 0.0 # um^2/ms take Fluo5's


# For Fricks Ca imaging with bAP (Normalization of Ca 2! Signals by Small Oblique Dendrites of
# CA1 Pyramidal Neurons)
#Constants https://doi.org/10.3390/biom11030343
tot_OGB1 = 0.200 #200 uM
kf_OGB1 =  365 # 1/mM/ms
kb_OGB1 = 73e-3 #1/ms, Kd roughly 200 nM
OGB1Diff = 0.27 # um^2/ms taken from https://doi.org/10.1016/j.ceca.2007.11.007
diffusions = {"CaM": camDiff,
              "Calb": calbDiff,
              "Mg Green": mggreenDiff,
              "Fluo3":fluo3Diff,
              "OGB1":OGB1Diff}
              
membrane_shell_width = .1


