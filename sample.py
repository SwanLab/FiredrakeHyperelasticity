from firedrake import *

mesh = UnitSquareMesh(1, 1)
V = FunctionSpace(mesh, "DG", 0)
f = Function(V)
f.interpolate(sin(SpatialCoordinate(mesh)[0]))

outfile = VTKFile("output_sample.pvd")
outfile.write(f)
