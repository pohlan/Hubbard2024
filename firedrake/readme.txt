Fiddling around with firedrake.

heat.py:
    Time-steps the heat equation in 1, 2, or 3 dimensions, starting from a gaussian bump initial condition.
    Dirichlet boundary conditions in all scenarios.
    Generates paraview output.

mctigueX:
    mctigue_mesh.py - Makes a cylindrical mesh with a spherical void at the center. This is a common model for thinking about volcano deformation (the cylindrial void represents a magma chamber). An analytical approximation of the deformation caused in this scenario was derived by a person named McTigue.
    mctigue.py - This runs a linear-elastic deformation model where the inside of the buried sphere is pressurized, the top of the cylinder is a free surface, and the bottom and sides of the cylinder are held at zero deformation.
    Generates a paraview output and a plot.
