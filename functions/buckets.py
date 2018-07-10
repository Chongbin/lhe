import math

def mass(momentum):
        mass = math.sqrt(momentum[3]**2-momentum[2]**2-momentum[1]**2-momentum[0]**2)
        return mass

# Filter out the quarks in the data
def quarks(event):
        quarks = []
        for i in range(len(event)):
                if event[i][0] >= -5 and event[i][0] <= 5 and event[i][1] == 1:
                        quarks.append(event[i])
        return quarks

# Generate the power set for a given set
def pset(order):
        powerset = []
        for i in xrange(2**len(order)):
                # for each i (in binary): include the elements whose
                # corresponding digit is 1 in the subset
                subset = tuple([x for j,x in enumerate(order) if (i >> j) & 1])
                powerset.append(subset)
        return powerset

# Generate all possible buckets with exactly one b-jet in it and at least one other jet
def top_pset(quarks):
        order = range(len(quarks))
        powerset = pset(order)
        top_powerset = []
        for i in range(len(powerset)):
                count = 0
                top_set = []
                for j in range(len(powerset[i])):
                        index = powerset[i][j]
                        if quarks[index][0] == 5 or quarks[index][0] == -5:
                                count += 1
                        top_set.append(quarks[index])
                if count == 1 and len(top_set) > 1:
                        top_powerset.append(top_set)
        return top_powerset

def comb_buckets(event):
        quark = quarks(event)
        pset = top_pset(quark)
        bucket1 = []
        bucket2 = []
	mass1 = 0
	mass2 = 0
        delta = float('inf')
	t0 = []
        for i in range(len(pset)):
                p1 = [0.0,0.0,0.0,0.0]
                p2 = [0.0,0.0,0.0,0.0]
                #temp_event = quark
                b1 = pset[i]
                b2 = [x for x in quark if x not in b1]
                for j in range(len(b1)):
                        for k in range(len(p1)):
                                p1[k] += b1[j][6+k]
                for j in range(len(b2)):
                        for k in range(len(p2)):
                                p2[k] += b2[j][6+k]
                m1 = mass(p1)
                m2 = mass(p2)
                d = 100 * ((m1 - 173.5) ** 2) + (m2 - 173.5) ** 2
                if d < delta:
                        delta = d
                        bucket1 = b1
                        bucket2 = b2
			mass1 = m1
			mass2 = m2
	if mass1 > 155 and mass1 < 200:
		t0.append('not t0')
	else:
		t0.append('t0')
	if mass2 > 155 and mass2 < 200:
		t0.append('not t0')
	else: 
		t0.append('t0')
        return bucket1, bucket2, t0

def label_buckets(event):
	buckets = list(comb_buckets(event)[0:2])
	t0 = comb_buckets(event)[2]
	labels = []
	# Loop over all possible two-jet combinations
	# See if any satisfies the W candidate condition
	for i in range(len(buckets)):
		if t0[i] == 't0':
			labels.append('t0')
		else:
			p_b = [0.0, 0.0, 0.0, 0.0]
			for j in range(len(buckets[i])):
				for k in range(len(p_b)):
					p_b[k] += buckets[i][j][6+k]
			m_b = mass(p_b)
			count = 0
			for j in range(len(buckets[i])-1):
				p1 = buckets[i][j][6:10]
				for k in range(j+1,len(buckets[i])):
					#if count == 1:
					#	break
					p2 = buckets[i][k][6:10]
					p = [sum(x) for x in zip(p1, p2)]
					m_jk = mass(p)
					if math.fabs(m_jk/m_b - 80.4/173.5) < 0.15:
						labels.append('tw')
						count += 1
						break #
				if count == 1:
					break
			if count == 1:
				continue
			labels.append('t_')
	return labels 

def top_buckets(event):
	quark = quarks(event) 
	buckets = list(comb_buckets(event)[0:2])
	labels = label_buckets(event)
	count = 0
	index = 0
	for i in range(len(labels)):
		if labels[i] == "t_":
			count += 1
			index = i
	if count == 0:
		buckets.append([])
	if count == 1:
		delta_bj = float('inf')
		m_bj = 0	
		new_bucket = [[],[]]
		# [1-index] only works when there are two buckets
		quark = [x for x in quark if x not in buckets[1-index]]
		for j in range(len(quark)-1):
			if quark[j][0] == 5 or quark[j][0] == -5:
				p1 = quark[j][6:10]
				for k in range(j+1,len(quark)):
					p2 = quark[k][6:10]
					p = [sum(x) for x in zip(p1, p2)]
					m_b = mass(p) 
					delta = math.fabs(m_b - 145)
					if m_b > 155:
						delta = float('inf')
					if delta < delta_bj:
						delta_bj = delta
						new_bucket[0] = quark[j]
						new_bucket[1] = quark[k]
						buckets[index] = new_bucket
						m_bj = m_b
		if m_bj <= 75 or m_bj >= 155:
			labels[index] = 't0'
                temp = [x for x in quark if x not in buckets[0]]
                t_extra = [x for x in temp if x not in buckets[1]]
		if len(t_extra) == 0:
			buckets.append([])
		else:
                	buckets.append(t_extra)
	if count > 1:
		delta_bj = float('inf')
		new_bucket1 = [[],[]]
		new_bucket2 = [[],[]]
		m_bj1 = 0
		m_bj2 = 0
		for j in range(len(quark)-1):
			if quark[j][0] == 5 or quark[j][0] == -5:
                        	p11 = quark[j][6:10]
                        	for k in range(j+1,len(quark)):
					bucket1 = [quark[j], quark[k]]
                                        p12 = quark[k][6:10]
                    	                p1 = [sum(x) for x in zip(p11, p12)]
                            	        m_b1 = mass(p1)
                              	        delta_bj1 = math.fabs(m_b1 - 145)
					if m_b1 > 155:
                                	        delta_bj1 = float('inf')
					quark_remain = [x for x in quark if x not in bucket1]
					for n in range(len(quark_remain)-1):
						if quark_remain[n][0] == 5 or quark_remain[n][0] == -5:
							p21 = quark_remain[n][6:10]
							for m in range(n+1,len(quark_remain)):
								p22 = quark_remain[m][6:10]
								p2 = [sum(x) for x in zip(p21, p22)]
								m_b2 = mass(p2)
								delta_bj2 = math.fabs(m_b2 - 145)
								if m_b2 > 155:
									delta_bj2 = float('inf')
								delta = delta_bj1 + delta_bj2
								if delta < delta_bj:
									delta_bj = delta
									buckets[0] = bucket1
									new_bucket2[0] = quark_remain[n]
                                              		                new_bucket2[1] = quark_remain[m]
									buckets[1] = new_bucket2 	
									m_bj1 = m_b1
									m_bj2 = m_b2
		print "m1: " + str(m_bj1)
		print "m2: " + str(m_bj2)
		if m_bj1 <= 75 or m_bj1 >= 155:
                        labels[0] = 't0'
		if m_bj2 <= 75 or m_bj2 >= 155:
                        labels[1] = 't0'
		temp = [x for x in quark if x not in buckets[0]]
		t_extra = [x for x in temp if x not in buckets[1]]
                if len(t_extra) == 0:
                        buckets.append([])
                else:
                        buckets.append(t_extra)
	return labels, buckets 
			
