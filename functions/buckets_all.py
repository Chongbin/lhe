import math

def mass(momentum):
  mass_sq = momentum[3]**2-momentum[2]**2-momentum[1]**2-momentum[0]**2
  #if math.fabs(mass_sq) < 0.01:
  if mass_sq < 0:
    mass = 0
  else:
    mass = math.sqrt(mass_sq)
  return mass

# Filter out the quarks in the data
def quarks(event):
  quarks = []
  bottom = 0
  ord_index = 0
  for i in range(len(event)):
    if ((event[i][0] >= -5 and event[i][0] <= 5) or math.fabs(event[i][0]) == 21) and event[i][1] == 1:
      ord_index += 1
      extrainfo = event[i][0:6]*1
      extrainfo.append(ord_index)  #particle order added
      quarks.append([event[i][6:10],extrainfo])
      if (math.fabs(event[i][0]) == 5):
        bottom += 1
  if bottom == 2:
    return quarks
  else:
    return []

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
      if quarks[index][1][0] == 5 or quarks[index][1][0] == -5:
        count += 1
      top_set.append(quarks[index])
    if count == 1 and len(top_set) >= 1: #modified to allow single b bucket too
      top_powerset.append(top_set)
    #print top_powerset 
  return top_powerset

def comb_buckets(event):
  quark = quarks(event)
  if len(quark) > 0:
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
        p1 = [sum(x) for x in zip(p1, b1[j][0])]
      for j in range(len(b2)):
        p2 = [sum(x) for x in zip(p2, b2[j][0])]
      m1 = mass(p1)
      m2 = mass(p2)
      d = 100 * ((m1 - 173.5) ** 2) + (m2 - 173.5) ** 2
      piD1 = [x[-1][-1] for x in b1]
      piD2 = [x[-1][-1] for x in b2]
      #print "del: %.3f"%d, piD1, "\t", piD2
      if d < delta:
        delta = d
        bucket1 = b1
        bucket2 = b2
        mass1 = m1
        mass2 = m2
    #piD1 = [x[-1][0] for x in bucket1]
    #piD2 = [x[-1][0] for x in bucket2]
    #^#print "del: ", delta
    #^#print "B1: ", piD1
    #^#print "B2: ", piD2
    if mass1 > 155 and mass1 < 200:
      t0.append('not t0')
    else:
      t0.append('t0')
    if mass2 > 155 and mass2 < 200:
      t0.append('not t0')
    else: 
      t0.append('t0')
    return bucket1, bucket2, t0
  else:
    return []

def label_buckets(event):
  combs = comb_buckets(event)
  if len(combs) > 0:
    buckets = list(combs[0:2])
    t0 = combs[2]
    labels = []
    masses = []
    # Loop over all possible two-jet combinations
    # See if any satisfies the W candidate condition
    for i in range(len(buckets)):
      p_b = [0.0, 0.0, 0.0, 0.0]
      for j in range(len(buckets[i])):
        p_b = [sum(x) for x in zip(p_b, buckets[i][j][0])]
      m_b = mass(p_b)
      #piD = [x[-1][-1] for x in buckets[i]]
      #if ((m_b < 10) and (m_b > -1)): 
      #  print "massB: ", m_b, "\tbucket: ", piD
      ratio_min = float('inf')
      m_jk_min = 0
      for j in range(len(buckets[i])):
        p1 = buckets[i][j][0]
        for k in range(len(buckets[i])):
          if (k < j): continue
          p2 = buckets[i][k][0]
          p = [sum(x) for x in zip(p1, p2)]
          m_jk = mass(p)
          diff = math.fabs(m_jk/m_b - 80.4/173.5)
	  #^#print "dd: ", diff, "\tratmin: ", ratio_min
          if diff < ratio_min:
            m_jk_min = m_jk
            ratio_min = diff
	    #^#print "temp mass W: ", m_jk_min
      if (t0[i] == "t0"): labels.append("t0")
      else:
        if (ratio_min < 0.15):
          labels.append('tw')
        else:
          labels.append('t_')
      #^#print "W_cand_mass: ", m_jk_min
      #^#print "ratio_min: ", ratio_min
      masses.append([m_jk_min, m_b, ratio_min])
    return labels, masses 
  else:
    return []

def top_buckets(event):
  quark = quarks(event)
  if len(quark) > 0:
    buckets = list(comb_buckets(event)[0:2])
    labels = label_buckets(event)[0]
    count = 0
    index = 0
    for i in range(len(labels)):
      if labels[i] == "t_":
        count += 1
        index = i
    if count == 0:
      buckets.append([])
    elif count == 1:
      delta_bj = float('inf')
      m_bj = 0  
      new_bucket = [[],[]]
      # [1-index] only works when there are two buckets
      quark = [x for x in quark if x not in buckets[1-index]]
      buckets[index] = []
      ##print quark
      ##print buckets[1-index]
      for j in range(len(quark)):
        if quark[j][1][0] == 5 or quark[j][1][0] == -5:
          p1 = quark[j][0]
          for k in range(len(quark)):
            if (k == j): continue
            ##print "}}}}", k, "\t", j
            p2 = quark[k][0]
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
            ##print "~~~~~~~~~~~~~", delta_bj , "\t", quark[k][1]
      ##print "```````````````````", delta_bj, "\t", new_bucket, "\t", m_bj

      if ((m_bj <= 75) or (m_bj >= 155)):
        labels[index] = 't0'
      temp = [x for x in quark if x not in buckets[0]]
      t_extra = [x for x in temp if x not in buckets[1]]
      #~#print "EXTRA", t_extra
      if len(t_extra) == 0:
        buckets.append([])
      else:
        buckets.append(t_extra)
    elif count == 2:
      delta_bj = float('inf')
      new_bucket2 = [[],[]]
      m_bj1 = 0
      m_bj2 = 0
      buckets = [[],[]]
      for j in range(len(quark)):
        if quark[j][1][0] == 5 or quark[j][1][0] == -5:
          p11 = quark[j][0]
          for k in range(len(quark)):
            if (abs(quark[k][1][0]) == 5): continue
            bucket1 = [quark[j], quark[k]]
            p12 = quark[k][0]
            p1 = [sum(x) for x in zip(p11, p12)]
            m_b1 = mass(p1)
            delta_bj1 = math.fabs(m_b1 - 145)
            if m_b1 > 155:
              delta_bj1 = float('inf')
            quark_remain = [x for x in quark if x not in bucket1]
            for n in range(len(quark_remain)):
              if quark_remain[n][1][0] == 5 or quark_remain[n][1][0] == -5:
                p21 = quark_remain[n][0]
                for m in range(len(quark_remain)):
                  if (m == n): continue
                  p22 = quark_remain[m][0]
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
      
      
      labels[0] = 't0' if ((m_bj1 <= 75) or (m_bj1 >= 155)) else 't_'
      labels[1] = 't0' if ((m_bj2 <= 75) or (m_bj2 >= 155)) else 't_'
      temp = [x for x in quark if x not in buckets[0]]
      t_extra = [x for x in temp if x not in buckets[1]]
      #~#print "EXTRA", t_extra
      if len(t_extra) == 0:
        buckets.append([])
      else:
        buckets.append(t_extra)
    return labels, buckets 
  else:
    return []    
