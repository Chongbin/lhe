import sys
sys.path.insert(0,'../functions')
import math
import ROOT
import graphing
import buckets_all

data = graphing.split('../tt_hadronic.lhe')

mass_tw = []
mass_t_ = []
mass_t0 = []
mass_extra = []
pt_tw = []
pt_t_ = []
pt_t0 = []
pt_extra = []
eta_tw = []
eta_t_ = []
eta_t0 = []
eta_extra = []
m_jk = []
m_b = []
ratio = []

for i in range(len(data)):
        result = buckets_all.top_buckets(data[i])
        if len(result) > 0:
                labels = result[0]
                bucket = result[1]
                ratios = buckets_all.label_buckets(data[i])[1]
		if len(ratios) != 0:
			m_jk.append(ratios[0][0])
			m_b.append(ratios[0][1])
			ratio.append(ratios[0][2])
		if len(ratios) == 2:
			m_jk.append(ratios[1][0])
                        m_b.append(ratios[1][1])
                        ratio.append(ratios[1][2])
                for n in range(len(bucket)-1):
                        p = [0.0, 0.0, 0.0, 0.0]
                        for j in range(len(bucket[n])):
                                p = [sum(x) for x in zip(p, bucket[n][j][0])]
                        m = buckets_all.mass(p)
                        pt = math.sqrt(p[0] ** 2 + p[1] ** 2)
                        pmag = math.sqrt(pt ** 2 + p[2] ** 2)
                        pangle = math.acos(p[2] / pmag)
                        eta = - math.log(math.tan(pangle/2))
                        if labels[n] == "tw":
                                mass_tw.append(m)
                                pt_tw.append(pt)
                                eta_tw.append(eta)
                        if labels[n] == "t_":
                                mass_t_.append(m)
                                pt_t_.append(pt)
                                eta_t_.append(eta)
                        if labels[n] == "t0":
                                mass_t0.append(m)
                                pt_t0.append(pt)
                                eta_t0.append(eta)
                if len(bucket[-1]) != 0:
                        for j in range(len(bucket[-1])):
                                p = bucket[-1][j][0]
                                m = buckets_all.mass(p)
                                pt = math.sqrt(p[0] ** 2 + p[1] ** 2)
                                pmag = math.sqrt(pt ** 2 + p[2] ** 2)
                                pangle = math.acos(p[2] / pmag)
                                mass_extra.append(m)
                                pt_extra.append(pt)
                                if pt < 0.001:
                                        if p[2] > 0:
                                                eta = 9.7
                                        if p[2] < 0:
                                                eta = -9.7
                                else:
                                        eta = - math.log(math.tan(pangle/2))
                                eta_extra.append(eta)
                if i % 1000 == 0:
                        print "DONE " + str(i)

# Mass
f_m_tw = ROOT.TFile("reconstructed_top_mass_tw.root",'RECREATE')
h_m_tw = ROOT.TH1F("h",'Mass of tw Buckets',150,0,300)
graphing.graph(h_m_tw, f_m_tw, mass_tw, 'Mass (GeV)', 'Frequency', 'reconstructed_top_mass_tw')

f_m_t_ = ROOT.TFile("reconstructed_top_mass_t_.root",'RECREATE')
h_m_t_ = ROOT.TH1F("h",'Mass of t_ Buckets',150,0,300)
graphing.graph(h_m_t_, f_m_t_, mass_t_, 'Mass (GeV)', 'Frequency', 'reconstructed_top_mass_t_')

f_m_t0 = ROOT.TFile("reconstructed_top_mass_t0.root",'RECREATE')
h_m_t0 = ROOT.TH1F("h",'Mass of t0 Buckets',150,0,300)
graphing.graph(h_m_t0, f_m_t0, mass_t0, 'Mass (GeV)', 'Frequency', 'reconstructed_top_mass_t0')

f_m_x = ROOT.TFile("reconstructed_top_mass_x.root",'RECREATE')
h_m_x = ROOT.TH1F("h",'Mass of the Extra Buckets',110,-1,10)
graphing.graph(h_m_x, f_m_x, mass_extra, 'Mass (GeV)', 'Frequency', 'reconstructed_top_mass_x')

# Transverse Momentum
f_pt_tw = ROOT.TFile("reconstructed_top_pt_tw.root",'RECREATE')
h_pt_tw = ROOT.TH1F("h",'Transverse Momentum of tw Buckets',250,0,1200)
graphing.graph(h_pt_tw, f_pt_tw, pt_tw, 'Pt (GeV)', 'Frequency', 'reconstructed_top_pt_tw')

f_pt_t_ = ROOT.TFile("reconstructed_top_pt_t_.root",'RECREATE')
h_pt_t_ = ROOT.TH1F("h",'Transverse Momentum of t_ Buckets',250,0,1200)
graphing.graph(h_pt_t_, f_pt_t_, pt_t_, 'Pt (GeV)', 'Frequency', 'reconstructed_top_pt_t_')

f_pt_t0 = ROOT.TFile("reconstructed_top_pt_t0.root",'RECREATE')
h_pt_t0 = ROOT.TH1F("h",'Transverse Momentum of t0 Buckets',250,0,1200)
graphing.graph(h_pt_t0, f_pt_t0, pt_t0, 'Pt (GeV)', 'Frequency', 'reconstructed_top_pt_t0')

f_pt_x = ROOT.TFile("reconstructed_top_pt_x.root",'RECREATE')
h_pt_x = ROOT.TH1F("h",'Transverse Momentum of the Extra Buckets',100,0,500)
graphing.graph(h_pt_x, f_pt_x, pt_extra, 'Pt (GeV)', 'Frequency', 'reconstructed_top_pt_x')

# Pseudorapidity
f_eta_tw = ROOT.TFile("reconstructed_top_eta_tw.root",'RECREATE')
h_eta_tw = ROOT.TH1F("h",'Pseudorapidity of tw Buckets',100,-10,10)
graphing.graph(h_eta_tw, f_eta_tw, eta_tw, 'Pseudorapidity', 'Frequency', 'reconstructed_top_eta_tw')

f_eta_t_ = ROOT.TFile("reconstructed_top_eta_t_.root",'RECREATE')
h_eta_t_ = ROOT.TH1F("h",'Pseudorapidity of t_ Buckets',100,-10,10)
graphing.graph(h_eta_t_, f_eta_t_, eta_t_, 'Pseudorapidity', 'Frequency', 'reconstructed_top_eta_t_')

f_eta_t0 = ROOT.TFile("reconstructed_top_eta_t0.root",'RECREATE')
h_eta_t0 = ROOT.TH1F("h",'Pseudorapidity of t0 Buckets',100,-10,10)
graphing.graph(h_eta_t0, f_eta_t0, eta_t0, 'Pseudorapidity', 'Frequency', 'reconstructed_top_eta_t0')

f_eta_x = ROOT.TFile("reconstructed_top_eta_x.root",'RECREATE')
h_eta_x = ROOT.TH1F("h",'Pseudorapidity of the Extra Buckets',100,-10,10)
graphing.graph(h_eta_x, f_eta_x, eta_extra, 'Pseudorapidity', 'Frequency', 'reconstructed_top_eta_x')

# Masses and ratios in the W candidate test
f_m_jk = ROOT.TFile("reconstructed_top_m_jk.root",'RECREATE')
h_m_jk = ROOT.TH1F("h",'Mass of the (possible) W candidate',150,0,300)
graphing.graph(h_m_jk, f_m_jk, m_jk, 'Mass (GeV)', 'Frequency', 'reconstructed_top_m_jk')

f_m_b = ROOT.TFile("reconstructed_top_m_b.root",'RECREATE')
h_m_b = ROOT.TH1F("h",'Mass of the Entire Buckets before Recalculation',150,0,300)
graphing.graph(h_m_b, f_m_b, m_b, 'Mass (GeV)', 'Frequency', 'reconstructed_top_m_b')

f_ratio = ROOT.TFile("reconstructed_top_ratio.root",'RECREATE')
h_ratio = ROOT.TH1F("h",'Mass Ratio Difference',120,-0.1,1.1)
graphing.graph(h_ratio, f_ratio, ratio, 'Mass (GeV)', 'Frequency', 'reconstructed_top_ratio')
