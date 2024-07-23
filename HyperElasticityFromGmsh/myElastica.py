from firedrake import *
from defcon import *
from   math import sin
from   math import cos
from   math import floor


class ElasticaProblem(BifurcationProblem):

    def __init__(self):
        self.mu = 27.4e9
        self.lam = 64.0e9       

    def boundary_conditions(self, V):
        return DirichletBC(V, 0.0, "on_boundary")
        
    def computeMesh(self):
        return IntervalMesh(1000, 0, 1)

    
    def computeResidual(self, theta, v):
        mu  = self.mu
        lam = self.lam

        F = inner(grad(theta), grad(v))*dx - lam*2*sin(theta)*v*dx + mu*cos(theta)*v*dx

        return F

    def solve(self):
        mesh = self.computeMesh()
        V = VectorFunctionSpace(mesh,'CG',1)
        bcs = self.boundary_conditions(V)
        
        u  = Function(V)
        du = TestFunction(V)        
        r  = self.computeResidual(u,du)

        solve(r==0,u,bcs)
        File('output.pvd').write(u)
        return{}

e = ElasticaProblem()        
e.solve()