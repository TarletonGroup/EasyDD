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