# convex_contact_patch
1. This package is used as a tutorial for our complementarity-based motion prediction algorithm. We have examples of "point", "cylinder", "cube" and "ellipsoid". In  this package, there are folders of (1) "convex_contact_patches_interface", which contains the lectures and the main scripts. (2) The "funjac" folder, which contains the function and jacobian evaluation file for each scenario. (3) the "symbolic_computation" folder, (4) the "utility" folder, (5) the "pathmexmaci64" folder and the (6) visualize folder. In the following, we will introduce the "convex_contact_patches_interface" folder in detail.

2. The "convex_contact_patches_interface" folder contains 

   (1) "lecture.mlx" and "lecture2.mlx" live scripts files. They will teach you the basics of our motion prediction algorithm and examples of main script files.            Please go through the files first.

   (2) Two main script files "interface_cuboid.m" and "interface_cylinder.m". Before running the main files, you need to provide the information of object shape,          mass, dimensions, friction parameters, time step length, initial state and configuration, and the unit.

   (3) The "planner_xxx.m" files provide the applied impulses to the object at each time step.

   (4) The "path.opt" file contains the options and parameters for the PATH solver.
