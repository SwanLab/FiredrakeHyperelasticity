gmsh -2 metaMaterial.geo -format msh2
python elasticityFromGmshFile.py
paraview output.pvd

#gmsh -2 immersed_domain.geo -format msh2
#python example.py
#paraview output.pvd

