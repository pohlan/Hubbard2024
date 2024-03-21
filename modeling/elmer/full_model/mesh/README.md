Read a contour of the area around the Hubbard glacier and build the `geo`input file that will be used to make the mesh
`python Contour2geo.py -r 150.0 -i hubbard_contour.dat -o hubbard_mesh.geo`

You can change the resolution of the mesh directly by modifying the **lc** variable in hubbard_mesh.geo

Then make the Elmer mesh:
`gmsh -1 -2 hubbard_mesh.geo hubbard_mesh.msh`
`ElmerGrid 14 2 hubbard_mesh.msh -autoclean`

Create a vtu file to visualise the mesh:
`ElmerGrid 14 5 hubbard_mesh.msh -autoclean`

Now move this folder one directory up, where the sif files are

------------------------------
If going parallel:
`ElmerGrid 14 2 hubbard_mesh.msh -autoclean -metis 4 0`

Create a vtu file to visualise the mesh:
`ElmerGrid 14 5 hubbard_mesh.msh -autoclean -metis 4 0`