try:
  from pytriqs.archive import *
  from pytriqs.gf.local import GfImFreq
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

input_json = sys_argv[1]
try:
  input_dict = json.loads(input_json[1:-1].replace("'",'"'))
except:
  print "ERROR: unparseable input" 
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
  Q_k_iw = A['lmdb'][Q+'_k_iw']
except:
  print "ERROR: HDF5 archive does not contain the specified quantity",Q+'_k_iw'
  quit()
  
beta = 4.0/T #!!!!!!!!
nk,nk,niws = numpy.shape(Q_k_iw)
assert niws % 2 ==0, "ERROR: must be fermionic Green's function with both negative and postivie freqs"
n_iw = niws/2
ntau = n_iw*3+1

Q_k_tau = {}
Q_k_tau['data'] = numpy.zeros((nk,nk,ntau),dtype=numpy.complex_)
Q_k_tau['mesh'] = [ 
    numpy.linspace(0.0,2*numpy.pi,nk,endpoint=False),
    numpy.linspace(0.0,2*numpy.pi,nk,endpoint=False),
    numpy.linspace(0,beta,ntau,endpoint=True)
]

Q_F_k = {}
Q_F_k['data'] = numpy.zeros((nk,nk),dtype=numpy.complex_)
Q_F_k['mesh'] = [ 
    numpy.linspace(0.0,2*numpy.pi,nk,endpoint=False),
    numpy.linspace(0.0,2*numpy.pi,nk,endpoint=False)
]

for kxi in range(nk):
    for kyi in range(nk):
        qiw = GfImFreq(indices = [0], beta = beta, n_points = n_iw, statistic='Fermion')        
        qtau = GfImTime(indices = [0], beta = beta, n_points = ntau, statistic='Fermion') 

        qiw.data[:,0,0] = Q_k_iw[kxi,kyi,:]
        Nc = 1
        known_coeff = TailGf(Nc,Nc,1,-1)
        known_coeff[-1] = numpy.zeros((Nc,Nc))
        nmax = qiw.mesh.last_index()
        starting_iw=14.0 
        max_order=5
        overwrite_tail=False
        nmin = int(((starting_iw * beta)/math.pi-1.0)/2.0) 
        nmin = max(nmin,1)
        qiw.fit_tail(known_coeff,max_order, nmin,nmax, overwrite_tail)
       
        qtau << InverseFourier(qiw)  
        Q_k_tau['data'][kxi,kyi,:] = qtau.data[:,0,0]
        Q_F_k['data'][kxi,kyi] =  - (beta/numpy.pi) * qtau.data[ntau/2,0,0].real

A[Q+'_k_tau'] = Q_k_tau
A[Q+'_F_k'] = Q_F_k
del A

print "SUCCESS!!"
       

