import ROOT

def readfile(filename):
	f = open(filename, "r")
	data = f.readlines()
	return data

# Return a three-dimensional list. [events][particles][properties]
def split(data):
	lst = []
	for i in range(len(data)):
        	if data[i].strip()  == "<event>":
                	num = int(data[i+1].split()[0])
                	for j in range(i+2, i+2+num):
				event = []
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

