# vim: set filetype=python ts=2 sw=2 sts=2 expandtab: 

import sys,random, argparse
sys.dont_write_bytecode=True

parser = argparse.ArgumentParser(
    description="Apply Scott-Knot test to data read from standard input")

class o():
  def __init__(i,**fields) : i.override(fields)
  def override(i,d): i.__dict__.update(d); return i
  def __repr__(i):
    return'{'+ str([k for k in sorted(i.__dict__.keys()) 
                   if not "_" in k]) + '}'

The=o(cohen=0.3, conf=0.01, cliffs=0.147, b=240)

p = parser.add_argument

p("--demo",default=False, action="store_true")
p("--cohen",   type=float, 
    default=The.cohen, metavar='N',
    help="too small if delta less than N*std of the data)")
p("--cliffs",type=float, default=The.cliffs, metavar="N",
    help="threshold for cliffs test: disable,small,med,large=0,0.147,0.33,0.474") 
p("--conf", type=float, default=The.conf, metavar="N",
    help="bootstrap tests with confidence 1-n")
p("--b", type=float, default=The.b, metavar="N",
    help="repeated samples for bootstrap")

args = parser.parse_args()
The.cohen  = args.cohen + 0
The.cliffs = args.cliffs + 0
The.conf   = args.conf + 0
The.b      = args.b + 0

# ----------------------------------------------
# cliffs delta
def cliffsDelta(lst1, lst2, dull = The.cliffs):
  "Returns true if there are more than 'dull' differences"
  def runs(lst):
    for j,two in enumerate(lst):
      if j == 0: one,i = two,0
      if one!=two:
        yield j - i,one
        i = j
      one=two
    yield j - i + 1,two  
  #---------------------
  m, n = len(lst1), len(lst2)
  lst2 = sorted(lst2)
  j = more = less = 0
  for repeats,x in runs(sorted(lst1)):
    while j <= (n - 1) and lst2[j] <  x: j += 1
    more += j*repeats
    while j <= (n - 1) and lst2[j] == x: j += 1
    less += (n - j)*repeats
  d= (more - less) / (m*n) 
  return abs(d)  > dull

# ----------------------------------------------
# bootstrap
class Num:
  "An Accumulator for numbers"
  def __init__(i,name,inits=None): 
    i.n = i.m2 = i.mu = 0.0
    i.rank = 0
    i.name = name
    i.all  = []
    i.adds(inits or [])
  def adds(i,lst):
    [i.add(x) for x in x in lst]
  def sd(i): 
    return (i.m2/(i.n - 1))**0.5
  def add(i,x):
    i.n   += 1  
    i.all + [x] 
    delta  = x - i.mu
    i.mu  += delta*1.0/i.n
    i.m2  += delta*(x - i.mu)
  def __add__(i,j):
    if not isintance(j,(list,tuple)):
	return i + [j]
    out = Num(i.name + "*")
    [out.adds[n.all] for n in j]
    return out
  def delta(i,j): 
    if i.sd() + j.sd():
      return (j.mu - i.mu)/((i.sd()/i.n + j.sd()/j.n)**0.5)
    else:
      return z.mu - y.mu

def bootstrap(y0,z0,conf=The.conf,b=The.b):
  "From p220 to 223 of Efron's introduction to the boostrap."
  def testStatistic(a,b): return Num("/",a).delta(Num("/",b))
  y, z   = Num("totaly",y0), Num("totalz",z0)
  x      = y + z
  tobs   = testStatistic(y,z)
  yhat   = [y1 - y.mu + x.mu for y1 in y.all]
  zhat   = [z1 - z.mu + x.mu for z1 in z.all]
  bigger = 0
  for i in range(b):
    if testStatistic(total(sampleWithReplacement(yhat)),
                     total(sampleWithReplacement(zhat))) > tobs:
      bigger += 1
  return bigger / b < conf 

# -----------------------------------------------
def rdiv(data,  # a list of class Nums
         all,   # all the data combined into one num
         div,   # function: find the best split
         big,   # function: rejects small splits
         same, # function: rejects similar splits
         epsilon): # small enough to split two parts
  def recurse(parts,all,rank=0):
    "Split, then recurse on each part."
    cut,left,right = maybeIgnore(div(parts,all,big,epsilon),
                                 same,parts)
    if cut: 
      # if cut, rank "right" higher than "left"
      rank = recurse(parts[:cut],left,rank) + 1
      rank = recurse(parts[cut:],right,rank)
    else: 
      # if no cut, then all get same rank
      for part in parts: 
        part.rank = rank
    return rank
  recurse(sorted(data),all)
  return data

def div(parts,all,big,epsilon):
  cut,left,right = None,None,None
  before, mu     =  0, all.mu
  for i,l,r in leftRight(parts,epsilon):
    if big(l.n) and big(r.n):
      n   = all.n * 1.0
      now = l.n/n*(mu- l.mu)**2 + r.n/n*(mu- r.mu)**2  
      if now > before:
        before,cut,left,right = now,i,l,r

def leftRight(parts,epsilon=The.epsilon):
  rights = {}
  n = j = len(parts) - 1
  while j > 0:
    rights[j] = parts[j]
    if j < n: rights[j] += rights[j+1]
    j -=1
  left = parts[0]
  for i,one in enumerate(parts):
    if i> 0:
      if parts[i]._median - parts[i-1]._median > epsilon:
        yield i,left,rights[i]
      left += one

