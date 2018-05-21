import ROOT
import math

f = open("v1_events.lhe", "r")
data = f.readlines()

def graph(h, rt, dataset, xlabel, ylabel, filename):
        for i in range(len(dataset)):
                h.Fill(dataset[i])
        c = ROOT.TCanvas('c', 'c', 800, 600)
        h.Draw()                
        h.GetXaxis().SetTitle(xlabel)
        h.GetYaxis().SetTitle(ylabel)
        c.Update()
        rt.Write()
        rt.Close()
        c.Print(filename + ".eps")

pt1 = []
pt2 = []
pz1 = []
pz2 = []
pangle1 = []
pangle2 = []
energy1 = []
energy2 = []
invamass1 = []
invamass2 = []

for i in range(len(data)):
	lst = data[i].strip().split()
        if lst[0] == "6" or lst[0] == "-6":
		if len(lst) < 10:
			continue
		pt = math.sqrt(float(lst[6])**2 + float(lst[7])**2)
		pz = float(lst[8])
		p = math.sqrt(pt ** 2 + pz ** 2)
		pangle = pz / p
		energy = float(lst[9])
		if lst[2] == "1":
			pt1.append(pt)
			pz1.append(pz)
			pangle1.append(pangle)
			energy1.append(energy)
		if lst[2] == '3':
			pt2.append(pt)
			pz2.append(pz)
			pangle2.append(pangle)
			energy2.append(energy)

for i in range(len(data)):
	if data[i].strip()  == "<event>":
		num = int(data[i+1].split()[0])
		momentum1 = [0.0,0.0,0.0,0.0]
                momentum2 = [0.0,0.0,0.0,0.0]
		for j in range(i+2, i+2+num):
			lst = data[j].strip().split()
			if (lst[0] == "6" or lst[0] == "-6") and lst[2] == "1":
				for k in range(len(momentum1)):
					momentum1[k] += float(lst[6+k])
			if (lst[0] == "6" or lst[0] == "-6") and lst[2] == "3":
				for k in range(len(momentum2)):
                                        momentum2[k] += float(lst[6+k])
		mass1 = math.sqrt(momentum1[3] ** 2 - momentum1[2] ** 2 - momentum1[1] ** 2 - momentum1[0] ** 2)
		mass2 = math.sqrt(momentum2[3] ** 2 - momentum2[2] ** 2 - momentum2[1] ** 2 - momentum2[0] ** 2)
		invamass1.append(mass1)
		invamass2.append(mass2)		
	
# Transverse Momentum
f_pt1 = ROOT.TFile("pt_pair1.root",'RECREATE')
h_pt1 = ROOT.TH1F("h",'Transverse Momentum of the First Pair',250,0,1200)
graph(h_pt1, f_pt1, pt1, 'Momentum (GeV)', 'Frequency', 'pt_pair1')

f_pt2 = ROOT.TFile("pt_pair2.root",'RECREATE')
h_pt2 = ROOT.TH1F("h",'Transverse Momentum of the Second Pair',250,0,1200)
graph(h_pt2, f_pt2, pt2, 'Momentum (GeV)', 'Frequency', 'pt_pair2')

# Z-momentum
f_pz1 = ROOT.TFile("pz_pair1.root",'RECREATE')
h_pz1 = ROOT.TH1F("h",'Z-Momentum of the First Pair',300,0,1500)
graph(h_pz1, f_pz1, pz1, 'Momentum (GeV)', 'Frequency', 'pz_pair1')

f_pz2 = ROOT.TFile("pz_pair2.root",'RECREATE')
h_pz2 = ROOT.TH1F("h",'Z-Momentum of the Second Pair',300,0,1500)
graph(h_pz2, f_pz2, pz2, 'Momentum (GeV)', 'Frequency', 'pz_pair2')

# Polar Angle
f_pa1 = ROOT.TFile("pangle_pair1.root",'RECREATE')
h_pa1 = ROOT.TH1F("h",'Polar Angle of the First Pair',200,-1,1)
graph(h_pa1, f_pa1, pangle1, 'Cosine of Polar Angle', 'Frequency', 'pangle_pair1')

f_pa2 = ROOT.TFile("pangle_pair2.root",'RECREATE')
h_pa2 = ROOT.TH1F("h",'Polar Angle of the Second Pair',200,-1,1)
graph(h_pa2, f_pa2, pangle2, 'Cosine of Polar Angle', 'Frequency', 'pangle_pair2')

# Energy
f_e1 = ROOT.TFile("energy_pair1.root",'RECREATE')
h_e1 = ROOT.TH1F("h",'Energy of the First Pair',300,0,1500)
graph(h_e1, f_e1, energy1, 'Energy (GeV)', 'Frequency', 'energy_pair1')

f_e2 = ROOT.TFile("energy_pair2.root",'RECREATE')
h_e2 = ROOT.TH1F("h",'Energy of the Second Pair',300,0,1500)
graph(h_e2, f_e2, energy2, 'Energy (GeV)', 'Frequency', 'energy_pair2')

# Invariant Mass
f_m1 = ROOT.TFile("mass_pair1.root",'RECREATE')
h_m1 = ROOT.TH1F("h",'Invariant Mass of the First Pair',300,0,3000)
graph(h_m1, f_m1, invamass1, 'Mass (GeV)', 'Frequency', 'mass_pair1')

f_m2 = ROOT.TFile("mass_pair2.root",'RECREATE')
h_m2 = ROOT.TH1F("h",'Invariant Mass of the Second Pair',300,0,3000)
graph(h_m2, f_m2, invamass2, 'Energy (GeV)', 'Frequency', 'mass_pair2')

'''
#Transverse Momentum
c_pt1 = ROOT.TCanvas('c', 'c', 800, 600)
h_pt1.Draw()		
h_pt1.GetXaxis().SetTitle('Momentum (GeV)')
h_pt1.GetYaxis().SetTitle('Frequency')
c_pt1.Update()
f_pt1.Write()
f_pt1.Close()
c_pt1.Print("pt_pair1.eps")

c_pt2 = ROOT.TCanvas('c', 'c', 800, 600)
h_pt2.Draw()
h_pt2.GetXaxis().SetTitle('Momentum (GeV)')
h_pt2.GetYaxis().SetTitle('Frequency')
c_pt2.Update()
f_pt2.Write()
f_pt2.Close()
c_pt2.Print("pt_pair2.eps")
'''
