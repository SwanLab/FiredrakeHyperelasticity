from firedrake import *

mesh = Mesh('openHole.msh')

#mesh = BoxMesh(15,5,5,5,1,1)


left = 1
V = VectorFunctionSpace(mesh,'CG',1)

bcs = DirichletBC(V,0,[12])
g = Constant((1e7,0))

u = Function(V)


(mu,lam) = (27.4e9,64.0e9)
strain = sym(grad(u))

internalEnergy = 0.5*(2*mu*inner(strain,strain) + lam*inner(div(u),div(u)))*dx
#externalWork = inner(g,u)*ds
externalWork = inner(g,u('+'))*dS(15)
#externalWork = g*u('+')*dS(12)

J = 0.5*internalEnergy - externalWork

du = TestFunction(V)
dJ = derivative(J,u,du)

solve(dJ==0,u,bcs)
#u.dat.data
File('output.pvd').write(u)
