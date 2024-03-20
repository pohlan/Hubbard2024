import fiona
import numpy as np
import matplotlib.pyplot as plt

shp = fiona.open("frenchjoe.gpkg", mode="r")
coords = np.array(list(shp.values())[1]['geometry']['coordinates']).squeeze()
shp.close()

# Resolution
lc = 200

# Trim last point (is identical to first)
coords = coords[:-1, :]

# Scale coords
coords /= 1000
lc /= 1000

with open("frenchjoe.geo", mode="w") as fd:
    fd.write("// This is a geo file created using the python script makegeo.py\n")
    fd.write("Mesh.Algorithm=5;\n")
    fd.write("// Element size\n")
    fd.write("lc = %.6f;\n" % lc)
    for i in range(coords.shape[0]):
        fd.write("Point(%d) = { %.3f, %.3f, 0.0, lc};\n" % (i+1, coords[i,0], coords[i,1]))

    for i in range(coords.shape[0]-1):
        fd.write("Line(%d) = {%d, %d};\n" % (i+1, i+1, i+2))
    fd.write("Line(%d) = {%d, %d};\n" % (coords.shape[0], coords.shape[0], 1))

    fd.write("Curve Loop(1) = {")
    for i in range(coords.shape[0]-1):
        fd.write("%d," % (i+1))
    fd.write("%d};\n" % (coords.shape[0]))

    fd.write("Plane Surface(2) = {1};\n")
