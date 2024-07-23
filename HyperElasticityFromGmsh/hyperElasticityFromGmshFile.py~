from firedrake import *
#from defcon import *

class HyperElasticProblem():

    def __init__(self):
        self.mu = Constant(1)#27.4e9
        self.lam = Constant(100)#64.0e9 

    def computeMesh(self):
        return Mesh('metaMaterial.msh')

    def computeInternalEnergy(self,z):
        quad_degree = 5
        (u,p) = split(z)
        n = self.n
        mu  = self.mu
        lam = self.lam        
        F = grad(u)+Identity(n)
        J = det(F)
        logJ = ln(J**2)/2
        thetaJ = logJ
        dX = dx(degree=quad_degree)
        e = (mu*inner(F,F) - 2*mu*logJ - (0.5/lam)*inner(p,p) + inner(p,thetaJ))*dX
        return e

    def computeExternalWork(self,z):
        (u,p) = split(z)
        g = Constant((-0.01,0))
        W = inner(g,u)*ds(12)
        return W
    
    def computeTotalEnergy(self,z):
        E = self.computeInternalEnergy(z)
        W = self.computeExternalWork(z)
        J = 0.5*E - W
        return J           

    def computeResidual(self,u,du):
        J  = self.computeTotalEnergy(u)
        dJ = derivative(J,u,du)
        return dJ
            
    def solve(self):

        mesh = self.computeMesh()
        self.n = mesh.geometric_dimension()
        
        V = VectorFunctionSpace(mesh,'CG',2)
        Q = FunctionSpace(mesh,'CG',1)
        Z = V*Q
        z = Function(Z)

        z.split()[0].rename('displacement')
        z.split()[1].rename('pressure')
        
        bcs = DirichletBC(Z.sub(0),0,[11])
        
        dz = TestFunction(Z)
        r  = self.computeResidual(z,dz)

        solve(r==0,z,bcs,solver_parameters={'snes_monitor':None,'ksp_monitor':None,})        
        File('output.pvd').write(*z.split())
        return{}

e = HyperElasticProblem()
e.solve()        
        
