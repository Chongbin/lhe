import ROOT
f = open("v1_events.lhe", "r")
data = f.readlines()
masses = []
f = ROOT.TFile("mass_top_quark.root",'RECREATE')
h = ROOT.TH1F("h",'Mass of Top Quark',100,0,300)       
for i in range(len(data)):
	lst = data[i].strip().split()
        if lst[0] == "6" or lst[0] == "-6":
		if len(lst) < 10:
			continue
		mass = float(lst[10])
		masses.append(mass)
		h.Fill(mass)
print masses[0:100]
c = ROOT.TCanvas('c', 'c', 800, 600)
h.Draw()		
h.GetXaxis().SetTitle('Mass (GeV)')
h.GetYaxis().SetTitle('Frequency')
c.Update()
f.Write()
f.Close()
c.Print("mass_top_quark.eps")
