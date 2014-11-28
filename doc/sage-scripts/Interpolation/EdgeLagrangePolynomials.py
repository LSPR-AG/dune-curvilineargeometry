# This file was *autogenerated* from the file EdgeLagrangePolynomials.sage.
from sage.all_cmdline import *   # import sage library
_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_6 = Integer(6); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4); _sage_const_0p00001 = RealNumber('0.00001')# This Sage code computes the Lagrange basis functions for triangles

u = var('u')
v = var('v')

vars = [u]
nDofs = [_sage_const_2 , _sage_const_3 , _sage_const_4 , _sage_const_5 , _sage_const_6 ]

# Produce Enumerated Edge
def enumeratedEdge(n) :
  PT = []
  for i in range(_sage_const_0 , n + _sage_const_1 ) : 
    PT.append(i)
  return PT;

# Produces sets of tetrahedral points points  
def edgePoints(n) :
  PT = enumeratedEdge(n);
  return [ PT[i] / n    for i in range(len(PT)) ]

# Returns a set of all standard 3D polynomials of all orders up to and including n
def dofEdgeSet(n) :
  PT = enumeratedEdge(n);
  return [u**PT[i]          for i in range(len(PT))]  

# Produces set of Lagrange interpolatory polynomials of order n
def Lagrange1DPolynomials(n) : 
  dofs = dofEdgeSet(n)
  dofNum = nDofs[n-_sage_const_1 ]
  p = edgePoints(n)
  V = matrix(dofNum, dofNum, _sage_const_1 /_sage_const_2 )

  for i in range(_sage_const_0 , len(dofs)) :
    for j in range(_sage_const_0 , len(dofs)) :
      V[i,j] = dofs[i].subs(u == p[j])
      
      

  return V.inverse() * matrix(dofs).transpose()
  

  
# Print all interpolatory polynomials
def printList(L) :
  for i in range(_sage_const_0 , L.nrows()) :
    print(L)[i]
  print("----------------------------------------")

printList(Lagrange1DPolynomials(_sage_const_1 ))
printList(Lagrange1DPolynomials(_sage_const_2 ))
printList(Lagrange1DPolynomials(_sage_const_3 ))
printList(Lagrange1DPolynomials(_sage_const_4 ))
printList(Lagrange1DPolynomials(_sage_const_5 ))





# Self-test all lagrange polynomials by plugging in the points
def polynomialSelfTest(n) : 
  dofs = dofEdgeSet(n)
  dofNum = nDofs[n-_sage_const_1 ]
  p = edgePoints(n)
  
  LPols = Lagrange1DPolynomials(n)
  
  return [LPols.subs(u == p[j]) for j in range(_sage_const_0 , dofNum)]

print("Self-test started")
print("Self-test passes if for each point there is exactly 1 polynomial that evaluates to 1.0 and all others evaluate to 0. Further, for each point, a different polynomial must evaluate to 1");
print("*******************************************************")
for n in range(_sage_const_0 , _sage_const_5 ) : 
  print("testing, ", n + _sage_const_1 , " which has DoF=", nDofs[n])
  PST = polynomialSelfTest(n+_sage_const_1 )
  
  for i in range(_sage_const_0 , nDofs[n]) :
    test_zeros = _sage_const_0 ;
    test_ones = _sage_const_0 ;
    ones_pos = _sage_const_0 ;
    for j in range(_sage_const_0 , nDofs[n]) :
      if (abs(PST[i][j][_sage_const_0 ]) < _sage_const_0p00001 )	: test_zeros += _sage_const_1 ;
      if (abs(PST[i][j][_sage_const_0 ] - _sage_const_1 ) < _sage_const_0p00001 )	:
	test_ones += _sage_const_1 ;
	ones_pos = j;
    print("-- test result for point ", i, " is ", test_zeros, " zeros and ", test_ones, " ones, with one at polynomial ", ones_pos)
print("*******************************************************")      


