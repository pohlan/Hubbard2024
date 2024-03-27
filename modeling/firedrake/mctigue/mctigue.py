import os

os.environ["OMP_NUM_THREADS"] = "1"

import matplotlib.pyplot as plt
from firedrake import *

# Load mesh
mesh = Mesh("mctigue.msh")

V = VectorFunctionSpace(mesh, "CG", 1)

# No displacement on bottom and sides (physical surface 2)
bcs = [DirichletBC(V, Constant([0, 0, 0]), 2)]

# Elastic parameters
G = 10e9
v = 0.25
E = G*2*(1+v)

# Pressure in chamber
P = Constant(10e6)

# Lame parameters
lam = Constant((E*v)/((1+v)*(1-2*v)))
mu = Constant(G)
Id = Identity(mesh.geometric_dimension())  # 2x2 Identity tensor


def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)


def sigma(u):
    return lam * div(u) * Id + 2 * mu * epsilon(u)


u = TrialFunction(V)
v = TestFunction(V)
a = inner(sigma(u), epsilon(v)) * dx

# Apply pressure boundary condition inside of chamber (physical surface 1)
L = dot(-P*FacetNormal(mesh), v) * ds(1)

uh = Function(V, name="u")
solve(a == L, uh, bcs=bcs)

File("u.pvd").write(uh)


# Plot surface deformation
if True:
    xs = np.linspace(0, 6e3, 100)
    ys = np.zeros_like(xs)
    zs = np.zeros_like(xs)

    dx, dy, dz = zip(*uh.at(list(zip(xs, ys, zs))))
    plt.plot(xs, dx, "k--", label="Horiz. disp.")
    plt.plot(xs, dz, "k-", label="Vert. disp.")
    plt.xlabel("Radial distance (m)")
    plt.ylabel("Displacement (m)")
    plt.legend()
    plt.show()
