try:
  from pytriqs.archive import *
  from pytriqs.gf.local import *
  import math
  import os.path
  import numpy
  import json
  import sys
  from itertools import product
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
  quantity = input_dict['quantity']    
except:
  print "ERROR: input does not specify quantity"
  quit()
   
try: 
  spin_argument = input_dict['spin_argument']
  orbital_argument = input_dict['orbital_argument']    
  spatial_argument = input_dict['spatial_argument']
  temporal_argument = input_dict['temporal_argument']        
except:
  print "ERROR: input does not specify all arguments, orbital, spatial and temporal"
  quit()

Q_name = quantity
for arg in [spin_argument, orbital_argument, spatial_argument, temporal_argument]:
  if arg!='':
    Q_name+='_'+arg

Qfft_name = quantity
for arg in [spin_argument, orbital_argument, spatial_argument, temporal_argument]:
  if arg=='k':
    Qfft_name+='_r'
  elif arg=='r':
    Qfft_name+='_k'
  elif arg!='':    
    Qfft_name+='_'+arg

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
dims = sum([dims_per_arg[arg] for arg in [spin_argument, orbital_argument, spatial_argument, temporal_argument]])
assert len(sha) == dims, "ERROR: incorrect number of dimensions in data for the given choice of arguments. len(sha)="+str(len(sha))+" dims="+str(dims)

first_spatial_index = sum([dims_per_arg[arg] for arg in [spin_argument, orbital_argument]])
spatial_indices = [first_spatial_index,first_spatial_index+1]
assert sha[first_spatial_index] == sha[first_spatial_index+1], "ERROR: it is assumed that the spatial grid is square, NxN points"
    
Qfft = Q.copy()
if not no_mesh:
  if spatial_argument == 'k':
    Qfft['mesh'][first_spatial_index] = range(sha[first_spatial_index])
    Qfft['mesh'][first_spatial_index+1] = range(sha[first_spatial_index])
  elif spatial_argument == 'r':
    Qfft['mesh'][first_spatial_index] = list(numpy.linscace(0,2.0*numpy.pi,len(sha[first_spatial_index]),endpoint=False))
    Qfft['mesh'][first_spatial_index+1] = list(numpy.linscace(0,2.0*numpy.pi,len(sha[first_spatial_index]),endpoint=False))
  else:
    print "ERROR: unknown spatial argument", spatial_argument
    quit()

index_ranges = []
for i in range(dims):
  if i in spatial_indices: continue
  index_ranges.append(range(sha[i]))
for inds in product(*index_ranges):
  indices = tuple(list(inds)[:first_spatial_index] + [slice(sha[first_spatial_index]),slice(sha[first_spatial_index])] + list(inds)[first_spatial_index:])
  if spatial_argument == 'k': XX=numpy.fft.ifft2
  if spatial_argument == 'r': XX=numpy.fft.fft2
  if no_mesh: 
    Qfft[indices]=XX(Q[indices])
  else:
    Qfft['data'][indices]=XX(Q['data'][indices])

lmdb = A['lmdb']
lmdb[Qfft_name] = Qfft
A['lmdb'] = lmdb
del A

print "SUCCESS!!"
       

