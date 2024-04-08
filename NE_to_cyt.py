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
        self.cyt_buffers = []
        self.ER_buffers = []
        self.reactions = []
        self.add_rxd_regions()
        self.add_species()
        self.add_diffusion()
        self.add_buffer_reactions(self.buffers)

    def outermost_shell(diam, shell_width):
        return 2*shell_width/diam


    def add_ER_shells(self, section, start_from, first_shell):
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
                self.ER_shells[name] = []
                
            inner = 1-(sum(self.factors[name])+new_factor)
            if inner < 0:
                inner = 0
            outer = 1-sum(self.factors[name])
            if i == 0 and inner == 0:
                self.ER_shells[name].append(rxd.Region(section, nrn_region='i',
                                                    name="%s_InsideShell_%d" %
                                                    (which_dend, i)))
            else:
                if i:
                    self.ER_shells[name].append(rxd.Region(section,
                                                        geometry=rxd.Shell(inner,
                                                                           outer),
                                                        name="%s_InsideShell_%d" %
                                                        (which_dend, i)))
                else:
                    self.ER_shells[name].append(rxd.Region(section,
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

    def add_cyt_shells(self, section, start_from, first_shell, how_many_shells):
        name = section.name()
        which_dend = name.replace("[", "").replace("]", "")
        for i in range(how_many_shells):
            if i == 0:
                self.cyt_shells[name] = []
            inner = (i+1)*first_shell
            outer = (i+2)*first_shell
            if i == 0:
                self.cyt_shells[name].append(rxd.Region(section,
                                                            nrn_region='o',
                                                            geometry=rxd.Shell(inner,
                                                                               outer),
                                                            name="%s_OutsideShell_%d" %
                                                            (which_dend, i)))
            else:
 
                self.cyt_shells[name].append(rxd.Region(section,
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
        self.ER_shells = OrderedDict()
        self.cyt_shells = OrderedDict()
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

            self.add_ER_shells(sec, 0, factor)
            self.add_cyt_shells(sec, 0, membrane_shell_width, 3)
  
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
        self.ip3 = rxd.Species(self.outside_shell_list, d=ip3Diff,
                               initial=ip3_init,
                               atolscale=1e-9)
        self.calr = rxd.Species(self.ER_regions,
                                initial=calr_tot-calr_bound)
        self.calrca = rxd.Species(self.ER_regions, initial=calr_bound)

        self.serca = rxd.Species(self.membrane_list, initial=serca)
        self.sercaca = rxd.Species(self.membrane_list, initial=sercaca)
        self.add_cytosol_buffers()

    def add_cytosol_buffers(self):
        self.buffers = OrderedDict()
        self.indicator = None
        self.add_calmodulin()
        calbDiff = self.params["calbDiff"]
        calbindin_tot = self.params["calbindin_tot"]
        calbca = self.params["calbca"]
        self.calb = rxd.Species(self.outside_shell_list, d=calbDiff,
                                initial=calbindin_tot-calbca,
                                name='Calbindin',
                                charge=0, atolscale=1e-9)
        self.calbca = rxd.Species(self.outside_shell_list, d=calbDiff,
                                  initial=calbca,
                                  name='CalbindinCa',
                                  charge=0, atolscale=1e-9)
        self.buffers["Calb"] = [self.calb, self.calbca]
        fixed_buffer_tot = self.params["fixed_buffer_tot"]
        fixed_buffer_ca = self.params["fixed_buffer_ca"]
        self.fixed = rxd.Species(self.outside_shell_list,
                                 initial=fixed_buffer_tot
                                 -fixed_buffer_ca,
                                 name='FixedBuffer',
                                 charge=0, atolscale=1e-9)
        self.fixedca = rxd.Species(self.outside_shell_list,
                                   initial=fixed_buffer_ca,
                                   name='FixedBufferCa',
                                   charge=0, atolscale=1e-9)
        camDiff = self.params["camDiff"]
        calmodulin_tot = self.params["calmodulin_tot"]
        camn = self.params["camn"]
        camc = self.params["camc"]
        self.cam = rxd.Species(self.outside_shell_list, d=camDiff,
                               initial=calmodulin_tot-camn-camc,
                               name='CaM', charge=0, atolscale=1e-9)
        self.camn = rxd.Species(self.outside_shell_list, d=camDiff, initial=camn,
                                name='CaMN',
                                charge=0, atolscale=1e-9)
        self.camc = rxd.Species(self.outside_shell_list, d=camDiff, initial=camc,
                                name='CaMC',
                                charge=0, atolscale=1e-9)
        self.buffers["CaM"] = [self.cam, self.camn, self.camc]


    def add_buffer_reactions(self, buffer_list):
        kf_fixed_b = self.params["kf_fixed_b"]
        kb_fixed_b = self.params["kb_fixed_b"]
        kf_camn = self.params["kf_camn"]
        kb_camn = self.params["kb_camn"]
        kf_camc = self.params["kf_camc"]
        kb_camc = self.params["kb_camc"]
        kf_calbindin = self.params["kf_calbindin"]
        kb_calbindin = self.params["kb_calbindin"]
        kf_magnesium_green = self.params["kf_magnesium_green"]
        kb_magnesium_green = self.params["kb_magnesium_green"]
        
        fixed_rxn = rxd.Reaction(self.fixed + self.ca, self.fixedca,
                                 kf_fixed_b,
                                 kb_fixed_b)
        self.reactions.append(fixed_rxn)
        rn = rxd.Reaction(self.cam + self.ca, self.camn, kf_camn,
                          kb_camn)
        rc = rxd.Reaction(self.cam + self.ca, self.camc, kf_camc, kb_camc)
        self.reactions.extend([rn, rc])
        calb_rxn = rxd.Reaction(self.calb + self.ca, self.calbca,
                                kf_calbindin,
                                kb_calbindin)
        self.reactions.append(calb_rxn)

        
    def _make_object_lists(self):
        self.shell_list = []
        self.outside_shell_list []
        self.ER_regions = []
        self.membrane_list = []
        for key in self.ER_shells.keys():
            self.shell_list.extend(self.ER_shells[key])
            self.ER_regions.extend(self.ER_shells[key])
        for key in self.cyt_shells.keys():
            self.shell_list.extend(self.cyt_shells[key])
            self.outside_shell_list.extend(self.cyt_shells[key])
            self.membrane_list.append(self.cyt_shells[key][0])

    def add_serca_and_leak(self):
        self.leak = []
        for key in self.membrane.keys():
            pump1 = rxd.Reaction(2*self.ca[self.cyt_shells[key][0]]
                                 +self.serca[self.cyt_shells[key][0]],
                                 self.sercaca[self.cyt_shells[key][0]],
                                 kfserca, kbserca)
            pump2 = rxd.Reaction(self.sercaca[self.cyt_shells[key][0]],
                                 2*self.ca[self.ER_shells[key][0]]
                                 +self.serca[self.cyt_shells[key][0]],
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
        which is roughly 30 mV (ghk with K and Na, ER membrane is pearmeable 
        to both to similar extent), or a rewrite to mod (better option, JJS)
        pca/pk of IP3R is 6 
        Foskett JK, White C, Cheung KH, Mak DO 2007. 
        Inositol trisphosphate receptor Ca2+ release channels. 
        Physiol Rev 87: 593–658
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
            er_region =	self.cyt_shells[key][0]
            ca_region =	self.ER_shells[key][0]
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
            self.ip3_rate = self.gip3r[membrane] * (minf * self.h_gate[membrane]) ** 3

            self.ip3r.append(rxd.MultiCompartmentReaction(self.ca[er_region],
                                                          self.ca[ca_region],
                                                          0.85*self.ip3_rate,
                                                          0.85*self.ip3_rate,
                                                          membrane=membrane))
            self.ip3r.append(rxd.MultiCompartmentReaction(self.k[er_region],
                                                          self.k[ca_region],
                                                          0.85*self.ip3_rate,
                                                          0.85*self.ip3_rate,
                                                          membrane=membrane))
        
            # IP3 degradation
            for shell in self.cyt_shells[key]:
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
        for sec_name in self.ER_shells.keys():
            for i, shell in enumerate(self.ER_shells[sec_name][:-1]):
                dname = sec_name.replace("[", "").replace("]", "")
                f = self.factors[sec_name][i]
                dr = rxd.Parameter(self.inside_borders[sec_name][i],
                                   name=dname, 
                                   value=lambda nd:
                                   nd.segment.diam/2/f)
                self.drs.append(dr)
                rxn = rxd.MultiCompartmentReaction(self.ca[shell],
                                                   self.ca[self.ER_shells[sec_name][i+1]], 
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.inside_borders[sec_name][i])
                self.diffusions.append(rxn)
                rxn = rxd.MultiCompartmentReaction(self.k[shell],
                                                   self.k[self.ER_shells[sec_name][i+1]], 
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.inside_borders[sec_name][i])
                self.diffusions.append(rxn)

            for i, shell in enumerate(self.cyt_shells[sec_name][:-1]):
                f = i+1
                dr = rxd.Parameter(self.outside_borders[sec_name][i],
                                   name=dname, 
                                   value=lambda nd:
                                   nd.segment.diam/2/f)
                self.drs.append(dr)
                rxn = rxd.MultiCompartmentReaction(self.ca[self.cyt_shells[sec_name][i+1]],
                                                   self.ca[shell],
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.outside_borders[sec_name][i])
                self.diffusions.append(rxn)
                rxn = rxd.MultiCompartmentReaction(self.k[self.cyt_shells[sec_name][i+1]],
                                                   self.k[shell],
                                                   c_unit*caDiff/dr, 
                                                   c_unit*caDiff/dr,
                                                   border=self.outside_borders[sec_name][i])
                self.diffusions.append(rxn)
                
                
