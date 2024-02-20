import gmsh
import numpy as np

## Parameters
r = 6e3  # Domain radius
z = 6e3  # Domain depth

sr = 1e3  # Pressure source radius
sz = 1.25e3  # Pressure source depth

res = 500

# GMSH boilerplate
gmsh.initialize()
gmsh.clear()

# Model name
gmsh.model.add("mctigue")

# Add sphere
v0 = gmsh.model.occ.addSphere(0, 0, -sz, sr)

# Add cylinder
v1 = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, -z, r)

# Cut out sphere from cylinder
v2 = gmsh.model.occ.cut([(3, v1)], [(3, v0)], removeTool=False, removeObject=False)[0][0][1]

# Sync occ objects
gmsh.model.occ.synchronize()

# Remove original sphere and cylinder volumes
# occ.cut will do this by default with remove...=True but that also removes
# object boundaries which are needed for boundary conditions. This just removes
# the volumes
gmsh.model.removeEntities([(3, v0), (3, v1)])

# Make physical surfaces for boundary conditons
ents = gmsh.model.getEntities(dim=2)

# First one is the sphere surface (found this out through trial and error)
gmsh.model.addPhysicalGroup(
    2, [ents[0][1]], tag=1
)

# Second two are the side and top of the cylinder
gmsh.model.addPhysicalGroup(
    2, [ents[1][1], ents[2][1]], tag=2
)

# Add the whole volume too, or it doesn't mesh for some reason
gmsh.model.addPhysicalGroup(
    3, [v2], tag=3
)

# Set all point resolutions to res
ents = gmsh.model.getEntities(dim=0)
gmsh.model.mesh.setSize(ents, res)

# Set points around sphere to res/4
ents = gmsh.model.getEntitiesInBoundingBox(-sr-1, -sr-1, -sz-sr-1, sr+1, sr+1, -sz+sr+1, dim=0)
gmsh.model.mesh.setSize(ents, res/4)

# Set points at top of domain to res/4
ents = gmsh.model.getEntitiesInBoundingBox(-r-1, -r-1, -1, r+1, r+1, 1, dim=0)
gmsh.model.mesh.setSize(ents, res/4)

gmsh.model.mesh.generate(3)

gmsh.write("mctigue.msh")

gmsh.finalize()