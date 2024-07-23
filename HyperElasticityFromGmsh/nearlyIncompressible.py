from firedrake import *
mesh = BoxMesh(20,5,5,20,5,5)
left = 1
quad_degree = 5
V = VectorFunctionSpace(mesh,'CG',2)
Q = FunctionSpace(mesh,'CG',1)
Z = V*Q

bcs = DirichletBC(Z.sub(0),0,left)
g = Constant((0,0,-0.01))

z = Function(Z)
#for zk,name in zip(z.split(),['displacement','pressure']):
#    zk.rename(name)

z.split()[0].rename('displacement')
z.split()[1].rename('pressure')

(u,p) = split(z)
mu = Constant(1)
lam = Constant(100)
strain = sym(grad(u))

n = mesh.geometric_dimension()
F = grad(u)+Identity(n)
J = det(F)
logJ = ln(J**2)/2
thetaJ = logJ
dX = dx(degree=quad_degree)
internalEnergy = (mu*inner(F,F) - 2*mu*logJ - (0.5/lam)*inner(p,p) + inner(p,thetaJ))*dX

#internalEnergy = (mu*inner(strain,strain) - 0.5/lam*inner(p,p) + inner(p,div(u)))*dX


#l = inner(u-f,u-f)*dx
externalWork = inner(g,u)*ds
J = internalEnergy - externalWork

test = TestFunction(Z)
dJ = derivative(J,z,test)



#assemble(internalEnergy)
#assemble(dJ)
#assemble(derivative(dJ,z,TrialFunction(Z)))
solve(dJ==0,z,bcs,solver_parameters={
    'snes_monitor':None,
    'ksp_monitor':None,
})
File('output.pvd').write(*z.split())
