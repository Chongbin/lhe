import ROOT

# Return a three-dimensional list. [events][particles][properties]
def split(filename):
	f = open(filename, "r")
        data = f.readlines()
	lst = []
	for i in range(len(data)):
        	if data[i].strip()  == "<event>":
                	num = int(data[i+1].split()[0])
			event = []
                	for j in range(i+2, i+2+num):
                        	particle = data[j].strip().split()
                        	for k in range(len(particle)):
                                	particle[k] = float(particle[k])
                        	event.append(particle)
			lst.append(event)
	return lst

def graph(h, rt, dataset, xlabel, ylabel, filename):
        for i in range(len(dataset)):
                h.Fill(dataset[i])
        c = ROOT.TCanvas('c', 'c', 800, 600)
        h.Draw()                
        h.GetXaxis().SetTitle(xlabel)
        h.GetYaxis().SetTitle(ylabel)
        c.Update()
        rt.Write()
	c.Print(filename + ".eps")
        rt.Close()

def graph2d(h, rt, dataset, xlabel, ylabel, filename):
        for i in range(len(dataset)):
                h.Fill(dataset[i][0],dataset[i][1])
        c = ROOT.TCanvas('c', 'c', 800, 600)
        h.Draw('TEXT')
        h.GetXaxis().SetTitle(xlabel)
        h.GetYaxis().SetTitle(ylabel)
        c.Update()
        rt.Write()
        c.Print(filename + ".eps")
        rt.Close()

def graph2d2(h, rt, dataset1, dataset2, xlabel, ylabel, filename, form = 'TEXT',statbin = True):
	if len(dataset1) != len(dataset2):
		print "Lengths of the two data sets are not equal"
	else:
        	for i in range(len(dataset1)):
                	h.Fill(dataset1[i],dataset2[i])
	        c = ROOT.TCanvas('c', 'c', 800, 600)
	        h.Draw(form)
		h.SetStats(statbin)
	        h.GetXaxis().SetTitle(xlabel)
	        h.GetYaxis().SetTitle(ylabel)
	        c.Update()
	        rt.Write()
	        c.Print(filename + ".eps")
	        rt.Close()
