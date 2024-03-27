import os
os.environ["OMP_NUM_THREADS"] = "1"

from firedrake import *

# Set dimension:
dim = 1  # 1 or 2 or 3

# Make mesh and initial condition
if(dim == 1):
    mesh = UnitIntervalMesh(20)
    x = SpatialCoordinate(mesh)[0]
    y = .5
    z = .5
elif(dim == 2):
    mesh = UnitSquareMesh(20, 20)
    x = SpatialCoordinate(mesh)[0]
    y = SpatialCoordinate(mesh)[1]
    z = .5
elif(dim == 3):
    mesh = UnitCubeMesh(20, 20, 20)
    x = SpatialCoordinate(mesh)[0]
    y = SpatialCoordinate(mesh)[1]
    z = SpatialCoordinate(mesh)[2]

u_init = exp(-100*sqrt((x-.5)**2 + (y-.5)**2 + (z-.5)**2)**2)

# Set up function space and functions
# "CG", 1 means first degree Lagrange elements 
V = FunctionSpace(mesh, "CG", 1)

# Temperature at times n and n+1
u_n = Function(V, name="u^{n}")
u_n1 = Function(V, name="u^{n+1}")

v = TestFunction(V)

# Apply initial condition
u_n.interpolate(u_init)

# diffusivity
alpha = .2

# Time stepping
dt = 1.0/1000

# Weak form of heat equation
# backward euler time step
F = (inner((u_n1 - u_n)/dt, v) +
     alpha*inner(grad(u_n1), grad(v)))*dx

# Set boundary conditions
# Sides are numbered by the Unit...Mesh function
if(dim == 1):
    bcs = [DirichletBC(V, (0.0), 1),
           DirichletBC(V, (0.0), 2)]
elif(dim == 2):
    bcs = [DirichletBC(V, (0.0), 1),
           DirichletBC(V, (0.0), 2),
           DirichletBC(V, (0.0), 3),
           DirichletBC(V, (0.0), 4)]
elif(dim == 3):
    bcs = [DirichletBC(V, (0.0), 1),
           DirichletBC(V, (0.0), 2),
           DirichletBC(V, (0.0), 3),
           DirichletBC(V, (0.0), 4),
           DirichletBC(V, (0.0), 5),
           DirichletBC(V, (0.0), 6)]

# Set up output
u_file = File("u.pvd")
u_file.write(u_n,time=0.0)

# Stop time and solving iterations
t_end = 0.25
for t in ProgressBar("Time step").iter(np.linspace(0.0, t_end, int(t_end/dt))):
    solve(F == 0, u_n1, bcs=bcs)
    u_n.assign(u_n1)
    u_file.write(u_n,time=t)
