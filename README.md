## EasyDD.m (11/11/2020)

3D Discrete Dislocation Plasticity. Is currently capable of simulating
nanoindentation, micropillar compression and cantilever bending.
Mixed-language model C, CUDA and Matlab so requires at least a C/C++
compatible compiler, CUDA computation is optional. Explicitly calculates
dislocation-dislocation interactions O(N^2).

## Requirements

- C/C++ compiler,
- CUDA compiler (optional),
- MATLAB (2019b+ preferred).

## Capabilities

- Cantilever bending,
- Pillar compression and extension,
- Nanoindentation,
- Inclusion-dislocation interaction.

## Folder structure

- input: contains inputs,
- output: contains outputs,
- scr: source files,
	- wip: contains folders for work in progress of the researches using the code,
- test: place for unit and integral tests,
- shadow_realm: compressed folder with old code kept around for posterity,
- publications: each publication is a folder containing the commit and any auxiliary code used to generate the publication.

## Data structure

- rn: (:,4) array of nodal positions (x, y, z, flag)
	- flag == 7, fixed nodes. They are not allowed to move.
	- flag == 6, are only allowed to move on the surface they live at.
	- flag == 67, are virtual nodes that are not allowed to move.
- vn: (:,3) array of nodal velocities (vx, vy, vz)
- fn: (:,3) array of nodal forces (fx, fy, fz)
- links: (:,8) array of links (idx1, idx2, bx, by, bz, nx, ny, nz)

## TODO

- Make vertices and faces an argument, if they are not defined by the input provide a default.
- In remesh_surf, make it so dislocations do not leave the domain via the fixed end depending on the simulation type.

### Known issues and improvement wish list

- Memory model is highly inefficient. There is a lot of naive dynamic resizing of arrays.
- Some matrices variables change in their number of columns are hard to track/debug (rn, connectivity).
- Some variables hang around unused.
- Collision, Separation and Remeshing can be antagonistic at higher dislocation densities and may end up repeating processes until the simulation progresses enough.
- FEM boundary conditions need reworking so they can be passed in as inputs rather than soft-coded into the FEM mesh builder.
- There is no sanity check on input parameters.
- Some performance-critical functions have no Matlab equivalent so a C/C++ compiler is strictly required.
- There are potential opportunities for further parallelisation in collision, separation and integration processes.
- Better modularisation would help increase the usability of the code.
- Slip system in links(:, 6:8) may be inaccurate for some dislocations.