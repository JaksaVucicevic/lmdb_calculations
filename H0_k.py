try:
  from pytriqs.archive import *
  from pytriqs.gf.local import *
  from py_expression_eval import Parser
  import math
  import os.path
  import numpy
  numpy.warnings.filterwarnings('ignore')
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


for q in ['filename','hamiltonian_terms','hamiltonian_params','model_id','params','nk']:
  try:
    globals()[q] = input_dict[q]
  except:
    print "ERROR: ",q," not found in input."
    quit()
try:
  nk = int(nk)
except:
  print 'nk cannot be converted to integer'

try: 
  A = HDFArchive(str(filename),'a')     
except:
  print "ERROR: HDF5 archive corrupted, or library problem"
  quit()

try: 
  lmdb = A['lmdb']
except:
  lmdb = {}

#################################################################3
def evaluate_term_coefficient(coeff, params):    
  parser = Parser()
  expr = parser.parse(str(coeff))
  vrs = expr.variables()
  has_j = ( 'j' in vrs )
  params_cpy = {'j': 1}
  for k,v in params.iteritems():
    try: 
      params_cpy[str(k)] = float(str(v)) 
    except:
      params_cpy[str(k)] = str(v) 
  try: 
    val = expr.evaluate( params_cpy )
  except:      
    print "Cannot evaluate term coefficient!!", coeff, "\n params:",params_cpy
    raise
    quit()       
  if ((len(vrs)>2) and has_j):
    print "Terms with complex coefficients can only have one variable, sorry. If additional constants are present, the code will give incorrect results!"
    quit()    
  return str(val)+("j" if has_j else "")

for param in hamiltonian_params:
  lmdb[param] = float(evaluate_term_coefficient(param, params))
lmdb['model'] = model_id
lmdb['method'] = "DMFT(IAIPT)"

terms = [ term for term in hamiltonian_terms if len(term[1])==2 ]

#term proto
# [ coeff, [ [ [0,0], 'c', 'up', True ], [ [0,0], 'f', 'up', False ],...   ]

field_ids = []
for term in terms:
  field_ids.extend([term[1][0][1],term[1][1][1]])
field_ids = list(set(field_ids))

nfields = len(field_ids)

nambu = False
for term in terms:
    dag1 = term[1][0][3]
    dag2 = term[1][1][3] 
    if (dag1 and dag2) or ((not dag1) and (not dag2)):
        nambu = True
        break

#we are assuming su2 symmetry!!!!        

ks = numpy.linspace(0,2.0*numpy.pi,nk, endpoint=False)

H0k = {}
if (1+nambu)*nfields==1:
  H0k['data'] = numpy.zeros(( nk, nk ), dtype=numpy.complex_)
  H0k['mesh'] = [ks, ks]
else:
  H0k['data'] = numpy.zeros(( (1+nambu)*nfields, (1+nambu)*nfields, nk, nk ), dtype=numpy.complex_)
  if nambu:
    H0k['mesh'] = [range(2*nfields),range(2*nfields),ks, ks]
  else:  
    H0k['mesh'] = [field_ids,field_ids,ks, ks]

for term in terms:
  if term[1][0][1]=='dn': continue 
  coeff = term[0]
  t = complex(evaluate_term_coefficient(coeff, params))

  r = numpy.array(term[1][1][0])-numpy.array(term[1][0][0])
  dag1 = term[1][0][3]
  dag2 = term[1][1][3] 
  if (1+nambu)*nfields==1: h0k = H0k['data'][:,:]
  else: 
    a = field_ids.index(term[1][0][1])
    b = field_ids.index(term[1][1][1])
    h0k = H0k['data'][a,b,:,:]

  for kxi, kx in enumerate(ks):
    for kyi, ky in enumerate(ks):
      if r[0]==0 and r[1]==0: val = t
      else: val = t*numpy.exp(1j*numpy.dot([kx,ky],r))
      if dag1 and not dag2:
        h0k[kxi, kyi] += val
        if nambu: H0k['data'][ nfields+a, nfields+b, kxi, kyi] -= val                    
      elif dag1 and dag2:
        H0k['data'][ a, b+nfields, kxi, kyi] += val
      elif ((not dag1) and (not dag2)):
        H0k['data'][ a+nfields, b, kxi, kyi] += val

counter = 0
for kxi, kx in enumerate(ks):
  for kyi, ky in enumerate(ks):
    if (1+nambu)*nfields==1:
      if H0k['data'][kxi, kyi].real < 0: counter+=1
    else:
      vals, vecs = numpy.linalg.eig(H0k['data'][ :, :, kxi, kyi])
      if numpy.amax(numpy.abs(vals.imag))>1e-10:
        print "inaginary eigenvalues!!! non-Hermitan Hamiltonian!!!"
        quit()
      for val in vals.real: 
        if val.real < 0.0:
          counter+=1
lmdb['n'] = counter/ ( nk**2.0 * (1+nambu) )

if nfields>1:
  if nambu:
    Q_name = "H0_n_n_k"
  else:
    Q_name = "H0_o_o_k"
else:
  Q_name = "H0_k"
lmdb[Q_name] = H0k
A['lmdb'] = lmdb

print "SUCCESS!"
