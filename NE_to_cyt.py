import os
from collections import OrderedDict
import numpy as np
from neuron import rxd, h

my_loc = os.path.dirname(os.path.abspath(__file__))
params_file = os.path.join(my_loc, "ca_params.py")
c_unit = 6.0221409e5


class InsideER:
    sections = []
    def __init__(self, params=params_file):
        ER = h.Section(name="endoplasmic_reticulum")
        #cm the same as cell membrane cm
        #rm ~ same as dendritic membrane
        #Ra ~ 100 ohmcm
        
        #the computational unit
        
        ER_comp.insert("pas")
        ER_comp.L = 100
        ER_comp.diam = 0.1
        ER_comp.Ra = 100
        ER_comp.e_pas = 90 #  https://doi.org/10.1073/pnas.1015686108 (in 37C)
        # in 20 C probably ~ 80 mV
        ER_comp.nseg = 101
        self.sections.append(ER)
        
        self.reactions = []
        self.add_rxd_regions()
        self.add_species()
        self.add_diffusion()

    def outermost_shell(diam, shell_width):
        return 2*shell_width/diam


    def add_inside_shells(self, section, start_from, first_shell):
        name = section.name()
        which_dend = name.replace("[", "").replace("]", "")
        i = start_from
        new_factor = first_shell
        last_shell = False
        while True:
            if i == 0:
                self.membrane[name] = rxd.Region(section,
                                                 name='%s_membrane' %
                                                 which_dend,
                                                 geometry=rxd.membrane())
                self.factors[name] = []
                self.inside_shells[name] = []
                
            inner = 1-(sum(self.factors[name])+new_factor)
            if inner < 0:
                inner = 0
            outer = 1-sum(self.factors[name])
            if i == 0 and inner == 0:
                self.inside_shells[name].append(rxd.Region(section, nrn_region='i',
                                                    name="%s_InsideShell_%d" %
                                                    (which_dend, i)))
            else:
                if i:
                    self.inside_shells[name].append(rxd.Region(section,
                                                        geometry=rxd.Shell(inner,
                                                                           outer),
                                                        name="%s_InsideShell_%d" %
                                                        (which_dend, i)))
                else:
                    self.inside_shells[name].append(rxd.Region(section,
                                                        nrn_region='i',
                                                        geometry=rxd.Shell(inner,
                                                                           outer),
                                                        name="%s_InsideShell_%d" %
                                                        (which_dend, i)))
            if i == 1:
                self.inside_borders[name] = []
            if i > 0:
                self.inside_borders[name].append(rxd.Region(section,
                                                     name='%s_InsideBorder_%d'
                                                     % (which_dend, i-1),
                                                     geometry=rxd.ScalableBorder(diam_scale=outer)))

            self.factors[name].append(new_factor)
            if last_shell:
                break
            if not inner:
                break
            if sum(self.factors[name]) + 2*new_factor >= self.fc:
                last_shell = True
                new_factor = (self.fc - sum(self.factors[name]))
            else:
                new_factor = 2*new_factor 
            i = i+1

    def add_outside_shells(self, section, start_from, first_shell, how_many_shells):
        name = section.name()
        which_dend = name.replace("[", "").replace("]", "")
        for i in range(how_many_shells):
            if i == 0:
                self.outside_shells[name] = []
            inner = (i+1)*first_shell
            outer = (i+2)*first_shell
            if i == 0:
                self.outside_shells[name].append(rxd.Region(section,
                                                            nrn_region='o',
                                                            geometry=rxd.Shell(inner,
                                                                               outer),
                                                            name="%s_OutsideShell_%d" %
                                                            (which_dend, i)))
            else:
 
                self.outside_shells[name].append(rxd.Region(section,
                                                           geometry=rxd.Shell(inner,
                                                                              outer),
                                                           name="%s_OutsideShell_%d" %
                                                           (which_dend, i)))
            if i == 1:
                self.outside_borders[name] = []
            if i > 0:
                self.outside_borders[name].append(rxd.Region(section,
                                                     name='%s_OutsideBorder_%d'
                                                     % (which_dend, i-1),
                                                     geometry=rxd.ScalableBorder(diam_scale=outer)))

            
    def add_rxd_regions(self):
        self.membrane = OrderedDict()
        self.inside_shells = OrderedDict()
        self.outside_shells = OrderedDict()
        self.inside_borders = OrderedDict()
        self.outside_borders = OrderedDict()
        self.factors = OrderedDict()  #  membrane_shell_width/max_radius
        secs_spines = OrderedDict()
        sec_list = []
        for sec in self.sections_rxd:
            sec_list.append(sec)
        dendrites = sec_list[:]
        all_geom = OrderedDict()
        membrane_shell_width = self.params["membrane_shell_width"]
        for sec in dendrites:
            print("Add rxd to %s" % sec.name())
            factor = 2*membrane_shell_width/sec.diam

            self.add_inside_shells(sec, 0, factor)
            self.add_outside_shells(sec, 0, membrane_shell_width, 3)
  
    def add_species(self):
        regions = self.shell_list 
        caDiff = self.params["caDiff"]
        self.ca = rxd.Species(regions, d=caDiff, name='ca', charge=2,
                              initial=lambda nd:
                              self.get_ca_init(nd.region.name),
                              atolscale=1e-9)
        
        self.k = rxd.Species(regions, d=caDiff, name='k', charge=1,
                             initial=lambda nd:
                             self.get_k_init(nd.region.name),
                             atolscale=1e-9)
        ip3Diff = self.params["ip3Diff"]
        ip3_init = self.params["ip3_init"]
        calr_tot = self.params["calr_tot"]
        calr_bound = self.params["calr_bound"]
        self.ip3 = rxd.Species(self.shell_list, d=ip3Diff, initial=ip3_init,
                               atolscale=1e-9)
        self.calr = rxd.Species(self.ER_regions,
                                initial=calr_tot-calr_bound)
        self.calrca = rxd.Species(self.ER_regions, initial=calr_bound)

        self.serca = rxd.Species(self.membrane_list, initial=serca)
        self.sercaca = rxd.Species(self.membrane_list, initial=sercaca)

    def _make_object_lists(self):
        self.shell_list = []
        self.ER_regions = []
        self.membrane_list = []
        for key in self.inside_shells.keys():
            self.shell_list.extend(self.inside_shells[key])
            self.ER_regions.extend(self.inside_shells[key])
        for key in self.outside_shells.keys():
            self.shell_list.extend(self.outside_shells[key])
            self.membrane_list.append(self.outside_shells[key][0])

    def add_serca_and_leak(self):
        self.leak = []
        for key in self.membrane.keys():
            pump1 = rxd.Reaction(2*self.ca[self.outside_shells[key][0]]
                                 +self.serca[self.outside_shells[key][0]],
                                 self.sercaca[self.outside_shells[key][0]],
                                 kfserca, kbserca)
            pump2 = rxd.Reaction(self.sercaca[self.outside_shells[key][0]],
                                 2*self.ca[self.inside_shells[key][0]]
                                 +self.serca[self.outside_shells[key][0]],
                                 kcatserca, 0)
            self.reactions.extend([pump1, pump2])
            gleak = rxd.Parameter(ca_region, value=lambda nd:
                                  self.params["gleak"])
            self.leak.append(rxd.MultiCompartmentReaction(self.ca_ER[er_region],
                                                          self.ca[ca_region],
                                                          gleak,
                                                          gleak,
                                                          membrane=membrane))

    def add_ip3r(self):
        """
        hopefully corresponding to Neymotin's model however extremely 
        difficult to judge. Apparently (triple checked) Neymotin used
        Wagner et al 2004's model and this implementation corresponds
        to Wagner and hopefully is correct.


        This either needs ghk added, to actually account for the ER potential,
        which is roughly 95 mV, or a rewrite to mod
        """
        self.gip3r = rxd.Parameter(self.cyt_er_membrane_list,
                                   initial=self.params["gip3r"])
        self.m = []
        self.n = []
        self.h_inf = []
        self.h_gate = rxd.State(self.cyt_er_membrane_list, initial=0.8)
        self.ip3rg = []
        self.ip3r = []
        ip3_init = self.params["ip3_init"]
        for i, key in enumerate(self.ER.keys()):
            er_region =	self.outside_shells[key][0]
            ca_region =	self.inside_shells[key][0]
            membrane = self.membrane[key]
            self.m.append(self.ip3[ca_region]/(self.ip3[ca_region] +
                                               self.params["Kip3"]))
            self.n.append(self.ca[ca_region]/(self.ca[ca_region]  +
                                              self.params["Kact"]))
            self.h_inf.append(self.params["k_inh"]/(self.ca[ca_region] +
                                                    self.params["k_inh"]))
            minf = self.m[i] * self.n[i]
        

            self.ip3rg.append(rxd.Rate(self.h_gate[membrane],
                                       (self.h_inf[i] - self.h_gate[membrane])
                                       /self.params["ip3rtau"]))
            self.k = self.gip3r[membrane] * (minf * self.h_gate[membrane]) ** 3

            self.ip3r.append(rxd.MultiCompartmentReaction(self.ca[er_region],
                                                          self.ca[ca_region],
                                                          self.k, self.k,
                                                          membrane=membrane))
        
            # IP3 degradation
            for shell in self.shells[key]:
                self.reactions.append(rxd.Rate(self.ip3,
                                               (ip3_init-self.ip3[shell])
                                               /ip3degTau,
                                               regions=shell,
                                               membrane_flux=False))
                
    def add_diffusion(self):
        self.diffusions = []
        caDiff = self.params["caDiff"]
        diffusions = self.params["diffusions"]
        self.drs = []
        for sec_name in self.inside_shells.keys():
            for i, shell in enumerate(self.inside_shells[sec_name][:-1]):
                dname = sec_name.replace("[", "").replace("]", "")
                f = self.factors[sec_name][i]
                dr = rxd.Parameter(self.inside_borders[sec_name][i],
                                   name=dname, 
                                   value=lambda nd:
                                   nd.segment.diam/2/f)
                self.drs.append(dr)
                rxn = rxd.MultiCompartmentReaction(self.ca[shell],
                                                   self.ca[self.inside_shells[sec_name][i+1]], 
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.inside_borders[sec_name][i])
                self.diffusions.append(rxn)
                rxn = rxd.MultiCompartmentReaction(self.k[shell],
                                                   self.k[self.inside_shells[sec_name][i+1]], 
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.inside_borders[sec_name][i])
                self.diffusions.append(rxn)

                
                f = i+1
                dr = rxd.Parameter(self.outside_borders[sec_name][i],
                                   name=dname, 
                                   value=lambda nd:
                                   nd.segment.diam/2/f)
                self.drs.append(dr)
                rxn = rxd.MultiCompartmentReaction(self.ca[self.outside_shells[sec_name][i+1]],
                                                   self.ca[self.outside_shells[sec_name][i]],
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.outside_borders[sec_name][i])
                self.diffusions.append(rxn)
                rxn = rxd.MultiCompartmentReaction(self.k[self.outside_shells[sec_name][i+1]],
                                                   self.k[self.outside_shells[sec_name][i]],
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.outside_borders[sec_name][i])
                self.diffusions.append(rxn)
                
                
