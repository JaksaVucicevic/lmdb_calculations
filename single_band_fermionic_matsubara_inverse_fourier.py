try:
  from pytriqs.archive import *
  from pytriqs.gf.local import *
  import math
  import os.path
  import numpy
  import json
  import sys
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
  T = A['lmdb']['T']
except:
  print "ERROR: HDF5 archive does not contain Temperature"
  quit()

try: 
  Q = input_dict['quantity']    
except:
  print "ERROR: input does not specify quantity"
  quit()
   
try: 
  spatial_argument = input_dict['spatial_argument']    
except:
  print "ERROR: input does not specify spatial_argument"
  quit()

try: 
  starting_iw = input_dict['starting_iw']
except:
  print "ERROR: input does not specify starting_iw"
  quit()

try: 
  starting_iw = float(starting_iw)
except:
  print "ERROR: starting_iw cannot be converted to float"
  quit()


Q_sa_iw_name = Q+'_'+spatial_argument+'_iw'
Q_sa_tau_name = Q+'_'+spatial_argument+'_tau'

try:
  Q_sa_iw = A['lmdb'][Q_sa_iw_name]
except:
  print "ERROR: HDF5 archive does not contain the specified quantity",Q_sa_iw_name
  quit()
  
beta = 1.0/T
try:
  nk,nk,niws = numpy.shape(Q_sa_iw['data'])
  Q_sa_iw = Q_sa_iw['data']
except:
  try: 
    nk,nk,niws = numpy.shape(Q_sa_iw)
  except:
    print "ERROR: Quantity in unrecognized format: ",Q_sa_iw_name
    quit()
    
assert niws % 2 ==0, "ERROR: must be fermionic Green's function with both negative and postivie freqs"
n_iw = niws/2
ntau = n_iw*3+1

Q_sa_tau = {}
Q_sa_tau['data'] = numpy.zeros((nk,nk,ntau),dtype=numpy.complex_)
if spatial_argument=='k':
  sa_values = numpy.linspace(0.0,2*numpy.pi,nk,endpoint=False)
elif spatial_argument=='r':
  sa_values = range(nk)
else:
  print "unrecognized spatial argument:", spatial_argument
  quit()
Q_sa_tau['mesh'] = [ 
  sa_values,
  sa_values,
  numpy.linspace(0,beta,ntau,endpoint=True)
]

for kxi in range(nk):
    for kyi in range(nk):
        qiw = GfImFreq(indices = [0], beta = beta, n_points = n_iw, statistic='Fermion')        
        qtau = GfImTime(indices = [0], beta = beta, n_points = ntau, statistic='Fermion') 

        qiw.data[:,0,0] = Q_sa_iw[kxi,kyi,:]
        Nc = 1
        known_coeff = TailGf(Nc,Nc,1,-1)
        known_coeff[-1] = numpy.zeros((Nc,Nc))
        nmax = qiw.mesh.last_index()
        max_order=5
        overwrite_tail=False
        nmin = int(((starting_iw * beta)/math.pi-1.0)/2.0) 
        nmin = max(nmin,1)
        qiw.fit_tail(known_coeff,max_order, nmin,nmax, overwrite_tail)
        tail0 = qiw.tail[0][0,0]  
        qiw[0,0] -= tail0
        qiw.fit_tail(known_coeff,max_order, nmin,nmax, overwrite_tail)
        qtau << InverseFourier(qiw)  
        Q_sa_tau['data'][kxi,kyi,:] = qtau.data[:,0,0]


lmdb = A['lmdb']
lmdb[Q_sa_tau_name] = Q_sa_tau
A['lmdb'] = lmdb
del A

print "SUCCESS!!"
       

