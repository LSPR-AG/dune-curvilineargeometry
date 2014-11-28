# This Sage code computes the Lagrange basis functions for triangles

u = var('u')
v = var('v')

vars = [u]
nDofs = [2, 3, 4, 5, 6]

# Produce Enumerated Edge
def enumeratedEdge(n) :
  PT = []
  for i in range(0, n + 1) : 
    PT.append(i)
  return PT;

# Produces sets of tetrahedral points points  
def edgePoints(n) :
  PT = enumeratedEdge(n);
  return [ PT[i] / n    for i in range(len(PT)) ]

# Returns a set of all standard 3D polynomials of all orders up to and including n
def dofEdgeSet(n) :
  PT = enumeratedEdge(n);
  return [u^PT[i]          for i in range(len(PT))]  

# Produces set of Lagrange interpolatory polynomials of order n
def Lagrange1DPolynomials(n) : 
  dofs = dofEdgeSet(n)
  dofNum = nDofs[n-1]
  p = edgePoints(n)
  V = matrix(dofNum, dofNum, 1/2)

  for i in range(0, len(dofs)) :
    for j in range(0, len(dofs)) :
      V[i,j] = dofs[i].subs(u == p[j])
      
      

  return V.inverse() * matrix(dofs).transpose()
  

  
# Print all interpolatory polynomials
def printList(L) :
  for i in range(0, L.nrows()) :
    print(L)[i]
  print("----------------------------------------")

printList(Lagrange1DPolynomials(1))
printList(Lagrange1DPolynomials(2))
printList(Lagrange1DPolynomials(3))
printList(Lagrange1DPolynomials(4))
printList(Lagrange1DPolynomials(5))





# Self-test all lagrange polynomials by plugging in the points
def polynomialSelfTest(n) : 
  dofs = dofEdgeSet(n)
  dofNum = nDofs[n-1]
  p = edgePoints(n)
  
  LPols = Lagrange1DPolynomials(n)
  
  return [LPols.subs(u == p[j]) for j in range(0, dofNum)]

print("Self-test started")
print("Self-test passes if for each point there is exactly 1 polynomial that evaluates to 1.0 and all others evaluate to 0. Further, for each point, a different polynomial must evaluate to 1");
print("*******************************************************")
for n in range(0, 5) : 
  print("testing, ", n + 1, " which has DoF=", nDofs[n])
  PST = polynomialSelfTest(n+1)
  
  for i in range(0, nDofs[n]) :
    test_zeros = 0;
    test_ones = 0;
    ones_pos = 0;
    for j in range(0, nDofs[n]) :
      if (abs(PST[i][j][0]) < 0.00001)	: test_zeros += 1;
      if (abs(PST[i][j][0] - 1) < 0.00001)	:
	test_ones += 1;
	ones_pos = j;
    print("-- test result for point ", i, " is ", test_zeros, " zeros and ", test_ones, " ones, with one at polynomial ", ones_pos)
print("*******************************************************")      


