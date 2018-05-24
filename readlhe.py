f = open("v1_events.lhe", "r")
data = f.readlines()
momenta = []
for i in range(len(data)):
	if data[i].strip()  == "<event>":
		num = int(data[i+1].split()[0])
		for j in range(i+2, i+2+num):
			momentum = data[j].strip().split()[6:11]
			for k in range(len(momentum)):
				momentum[k] = float(momentum[k])	 
			momenta.append(momentum)
print momenta[6]
result = []
count = 0
diffs = []
print "There are " + str(len(momenta)) + " particles in the dataset."
for i in range(len(momenta)):
	cal_mass = momenta[i][3] ** 2 - momenta[i][2] ** 2 - momenta[i][1] ** 2 - momenta[i][0] ** 2
	mass = momenta[i][4] ** 2
	diff = abs(cal_mass - mass)
	if diff < 1e-3:
		result.append(True)
		count += 1
	else:
		result.append(False)
		diffs.append(diff)
print "Number of particles in agreement with relavitity: " + str(count)
print result[0:100]
print diffs[0:100]
