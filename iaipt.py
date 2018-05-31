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
  print "about to import smart_scripts...",
  sys.path.insert(0,'/home/jaksa/TRIQS/source/smart_scripts')
  from smart_scripts import *
  print "done"
except:
  print "ERROR: not all requirements met."
  raise

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


for q in ['filename','hamiltonian_terms','hamiltonian_params','model_id','params','nk','starting_iw','max_its','min_its','accr','initial_mu']:
  try:
    globals()[q] = input_dict[q]
  except:
    print "ERROR: ",q," not found in input."
    quit()
try:
  filename = str(filename)
  nk = int(nk)
  starting_iw = float(starting_iw)
  max_its = int(max_its)
  min_its = int(min_its)
  accr = float(accr)
  initial_mu = float(initial_mu)
except:
  print 'nk cannot be converted to integer'

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
  assert not ((len(vrs)>2) and has_j), "Terms with complex coefficients can only have one variable, sorry. If additional constants are present, the code will give incorrect results!"
  return str(val)+("j" if has_j else "")

def get_evaluated_params(params):
  evaluated_params = {}
  for k,v in params.iteritems():  
    try:
      evaluated_params[k] = evaluate_term_coefficient(v, params)
    except:
      evaluated_params[k] = v
  return evaluated_params

evaluated_params = get_evaluated_params(params)

terms = [ term for term in hamiltonian_terms if len(term[1])==2 ]

#term proto
# [ coeff, [ [ [0,0], 'c', 'up', True ], [ [0,0], 'f', 'up', False ],...   ]

field_ids = []
for term in terms:
  field_ids.extend([term[1][0][1],term[1][1][1]])
field_ids = list(set(field_ids))

nfields = len(field_ids)
for term in terms:
    s1 = term[1][0][2]
    s2 = term[1][1][2]
    assert s1==s2, "no spin-orbit allowed!!!"
    dag1 = term[1][0][3]
    dag2 = term[1][1][3] 
    assert dag1 and not dag2, "no sc_field allowed in the Hamiltonian, and normal order is assumed!!!"
for term in terms:
    co11, f11, s1, _ = term[1][0]
    co21, f21, s1, _ = term[1][1]
    counter = 0
    if s1=='up':
      for term2 in terms:
        co12, f12, s2, _ = term2[1][0]
        co22, f22, s2, _ = term2[1][1]        
        if s2=='dn' and co11==co12 and co21==co22 and f11==f12 and f21==f22: counter+=1
      if counter!=1: assert "Hamiltonian is not su2 symmetric or is not sanitized (has identical terms)!"
      
ks = numpy.linspace(0,2.0*numpy.pi,nk, endpoint=False)
ss = ['up','dn']

def fill_H0k(mu, terms):
  H0k = numpy.zeros(( nk, nk, nfields, nfields ), dtype=numpy.complex_)
  params_copy = evaluated_params.copy()
  params_copy['mu'] = mu
  for coeff, [[co1,f1,s,_],[co2,f2,s,_]] in terms:
    if s=='dn': continue
    t = complex(evaluate_term_coefficient(coeff, params_copy))  
    r = numpy.array(co2)-numpy.array(co1)
    a = field_ids.index(f1)
    b = field_ids.index(f2)
    for kxi, kx in enumerate(ks):
      for kyi, ky in enumerate(ks):
        if r[0]==0 and r[1]==0: val = t
        else: val = t*numpy.exp(1j*numpy.dot([kx,ky],r))
        H0k[kxi, kyi, a, b] += val
  return H0k
 
get_H0k = partial(fill_H0k, terms=terms)

Us = {}
terms = [ term for term in hamiltonian_terms if len(term[1])==4 ]
for term in terms:  
  assert len(set([t[1] for t in term[1]]))==1, 'no interorbital interactions allowed!!'
  for (ti,t),flavor,dagger in zip(enumerate(term[1]),['up','dn','dn','up'],[True,True,False,False]):
    assert (t[2]==flavor) and (t[3]==dagger), 'interaction is not of required type!!'
  co = term[1][0][0]
  for t in term[1][1:]:
    assert t[0]==co, 'interaction must be local!!'
  field_id = term[1][0][1]  
  coef = float(evaluate_term_coefficient(term[0], params))
  if str(field_id) in Us.keys(): Us[str(field_id)]+=coef
  else: Us[str(field_id)]=coef

assert 'n' in params.keys(), "total n must be specified for this calculation"
assert 'T' in params.keys(), "T must be specified for this calculation"
n = float(evaluate_term_coefficient(params['n'], params))
T = float(evaluate_term_coefficient(params['T'], params))

for sol in [ "metal", "insulator" ]:
  globals()["dt_"+sol] , globals()["convergers_"+sol]  = iaipt_launcher( 
    [str(fid) for fid in field_ids],
    get_H0k,
    Us,
    n, 
    T, 
    starting_iw,
    ks,
    sol, #this is initial guess
    initial_mu,
    max_its,
    min_its, 
    accr,
    filename[:-3]+".%s.h5"%sol
  )

#=========================== when done prepare lmdb dictionary ================#

for sol in [ "metal", "insulator" ]:    
  dt = globals()["dt_"+sol] 
  cs = globals()["convergers_"+sol] 
  if any(c.diffs[-1]>accr for c in cs):
    print "calculation from metal failed to converge!!!"
    continue
  lmdb = {}
  lmdb['model'] = model_id
  lmdb['method'] = "DMFT(IAIPT)"
  lmdb['n'] = dt.get_n()
  invalid = False
  for p in hamiltonian_params:
    try:
      print "dt.",p,"=",vars(dt)[p] 
      lmdb[p] = vars(dt)[p]
    except:
      print "dt.",p,"taking from parameters"
      try:
        lmdb[p] = float(evaluate_term_coefficient(p, params))
      except:
        print "impossible to determine all model parameters"
        invalid=True
        continue  
  if invalid: continue
  H0k = {}
  if nfields==1:
    H0k['data'] = numpy.zeros(( nk, nk ), dtype=numpy.complex_)
    H0k['mesh'] = [ks, ks]
  else:
    H0k['data'] = numpy.zeros(( nfields, nfields, nk, nk ), dtype=numpy.complex_)
    H0k['mesh'] = [field_ids,field_ids, ks, ks]
  if nfields>1: 
    for a in range(nfields):
      for b in range(nfields):
        H0k['data'][a,b,:,:] = dt.H0k[:,:,a,b]
  else: H0k['data'][:,:] = dt.H0k[:,:,0,0]
  lmdb["H0_o_o_k" if nfields>1 else "H0_k"] = H0k

  Gk = {}
  if nfields==1:
    Gk['data'] = numpy.zeros(( nk, nk, 2*dt.niw ), dtype=numpy.complex_)
    Gk['mesh'] = [ks, ks, dt.iws.imag ]
  else:
    Gk['data'] = numpy.zeros(( nfields, nfields, nk, nk, 2*dt.niw ), dtype=numpy.complex_)
    Gk['mesh'] = [field_ids,field_ids, ks, ks, dt.iws.imag]
  if nfields>1:
    for a in range(nfields):
      for b in range(nfields):
        for iwi,iw in enumerate(dt.iws):
          Gk['data'][a,b,:,:,iwi] = dt.G_ab_k_iw[iwi,:,:,a,b]
  else: 
    for iwi,iw in enumerate(dt.iws):
      Gk['data'][:,:,iwi] = dt.G_ab_k_iw[iwi,:,:,0,0]
  lmdb["G_o_o_k_iw" if nfields>1 else "G_k_iw"] = Gk

  Sigmak = {}
  if nfields==1:
    Sigmak['data'] = numpy.zeros(( nk, nk, 2*dt.niw ), dtype=numpy.complex_)
    Sigmak['mesh'] = [ks, ks, dt.iws.imag ]
  else:
    Sigmak['data'] = numpy.zeros(( nfields, nfields, nk, nk, 2*dt.niw ), dtype=numpy.complex_)
    Sigmak['mesh'] = [field_ids,field_ids, ks, ks, dt.iws.imag]
  if nfields>1:
    for a in range(nfields):
      Sigmak['data'][a,a,:,:,:] = dt.Sigma_imp_iw[field_ids[a]].data[:,0,0]
  else: Sigmak['data'][:,:,:] = dt.Sigma_imp_iw[field_ids[0]].data[:,0,0]
  lmdb["Sigma_o_o_k_iw" if nfields>1 else "Sigma_k_iw"] = Sigmak


  A=HDFArchive(  filename[:-3]+".%s.h5"%sol, "a" )
  A['lmdb'] = lmdb
  del A 

#============================= if there is a single solution, copy the file for aoutomatic upload to lmdb, otherwise add hysteresis tags, and leave it to the user to upload ==================#

if all(c.diffs[-1]<accr for c in convergers_metal) and all(c.diffs[-1]<accr for c in convergers_insulator):
  if numpy.amax(numpy.abs(dt_metal.G_loc_iw.data-dt_insulator.G_loc_iw.data)):
    os.system("cp %s %s"%(filename[:-3]+".metal.h5", filename))
  else:
    for sol in [ "metal", "insulator" ]:
      A=HDFArchive(  filename[:-3]+".%s.h5"%sol, "a" )
      lmdb = A['lmdb']
      lmdb['hysteresis_tag'] = sol
      A['lmdb'] = lmdb
      del A 

print "SUCCESS!"
