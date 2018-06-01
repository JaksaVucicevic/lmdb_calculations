try:
  from pytriqs.archive import *
  from pytriqs.gf.local import *
  import math
  import os.path
  import numpy
  import json
  import sys
  from itertools import product
  from copy import deepcopy 
except:
  print "ERROR: not all requirements met."
  quit()

if len(sys.argv)<2:
  print "ERROR: Expected arguments."
  quit()

input_json = sys.argv[1]
try:
  print "raw input_json:",input_json
  input_json_str = str(input_json).replace("\\\"",'\"')
  input_dict = json.loads(input_json_str)
except:
  print "ERROR: unparseable input\n\n"
  print "input_json_str:",input_json_str  
  quit() 

try:
  filename = input_dict['filename']
except:
  print "ERROR: filename not found in input."
  quit()

if not os.path.exists(filename):
  print "ERROR: file not found"
  quit()
try: 
  A = HDFArchive(str(filename),'a')     
except:
  print "ERROR: HDF5 archive corrupted, or library problem"
  quit()

try: 
  lmdb = A['lmdb']
except:
  print "ERROR: HDF5 archive does not contain lmdb dict"
  quit()

try:
  T = lmdb['T']
except:
  print "ERROR: HDF5 archive does not contain Temperature"
  quit()

try: 
  quantity = input_dict['quantity']    
except:
  print "ERROR: input does not specify quantity"
  quit()
   
try: 
  spin_argument = input_dict['spin_argument']
  orbital_argument = input_dict['orbital_argument']    
  spatial_argument = input_dict['spatial_argument']
  w_min = float(input_dict['w_min'])
  w_max = float(input_dict['w_max'])
  n_points = int(input_dict['n_points'])
except:
  print "ERROR: input does not specify all arguments, orbital, spatial and temporal"
  quit()

Q_name = quantity
for arg in [spin_argument, orbital_argument, spatial_argument]:
  if arg!='':
    Q_name+='_'+arg
Q_name+='_iw'

Qpade_name = Q_name[:-3]+'_w'

try:
  Q = A['lmdb'][Q_name]
except:
  print "ERROR: HDF5 archive does not contain the specified quantity",Q_name
  quit()
  
try:
  sha = numpy.shape(Q['data'])
  no_mesh = False
except:
  try: 
    sha = numpy.shape(Q)
    no_mesh = True
  except:
    print "ERROR: Quantity in unrecognized format: ",Q_name
    quit()

dims_per_arg = {
 '': 0,
 's': 1,
 's_s': 2,
 'o_o': 2,
 'n_n': 2,
 'k': 2,
 'r': 2,
 'tau': 1,
 't': 1,
 'iw': 1,
 'w': 1
}
dims = sum([dims_per_arg[arg] for arg in [spin_argument, orbital_argument, spatial_argument, 'iw']])
assert len(sha) == dims, "ERROR: incorrect number of dimensions in data for the given choice of arguments. len(sha)="+str(len(sha))+" dims="+str(dims)
   
Qpade = deepcopy(Q)
if not no_mesh:  
  Qpade['mesh'][-1] = list(numpy.linspace(w_min,w_max,n_points,endpoint=True))
  Qpade['data'] = numpy.zeros(tuple(list(sha[:-1]) + [n_points]),dtype=numpy.complex_)
else:
  Qpade = numpy.zeros(tuple(list(sha) + [n_points]),dtype=numpy.complex_)
index_ranges = []
for i in range(dims-1):
  index_ranges.append(range(sha[i]))
for inds in product(*index_ranges):
  indices = tuple(list(inds) + [slice(sha[-1])])
  qiw = GfImFreq(indices = [0], beta = 1.0/T , n_points = sha[-1]/2, name = "whatever")
  if no_mesh: 
    qiw.data[:,0,0]=Q[indices]
  else:
    qiw.data[:,0,0]=Q['data'][indices]
  qw = GfReFreq(indices = [0], window = (w_min, w_max), n_points = n_points, name = "whatever")
  qw.set_from_pade(qiw, n_points = sha[-1]/2, freq_offset = 0.0)
  indices = tuple(list(inds) + [slice(n_points)])
  if no_mesh: 
    Qpade[indices] = qw.data[:,0,0]
  else:
    Qpade['data'][indices] = qw.data[:,0,0]

lmdb = A['lmdb']
lmdb[Qpade_name] = Qpade
A['lmdb'] = lmdb
del A

print "SUCCESS!!"
       

