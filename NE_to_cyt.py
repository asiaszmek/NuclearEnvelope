from neuron import h

ER = h.Section(name="endoplasmic_reticulum")
#cm the same as cell membrane cm
#rm ~ same as dendritic membrane
#Ra ~ 100 ohmcm

#the computational unit

ER.insert("pas")
ER.L = 100
ER.diam = 0.1
ER.Ra = 100
ER.e_pas = 90 #  https://doi.org/10.1073/pnas.1015686108 (in 37C)
# in 20 C probably ~ 80 mV
ER.nseg = 101
