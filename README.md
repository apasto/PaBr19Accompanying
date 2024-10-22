# Accompanying data to: "A geothermal application for GOCE satellite gravity data: modelling the crustal heat production and lithospheric temperature field in Central Europe"
by [Alberto Pastorutti](https://orcid.org/0000-0002-0279-5762) and [Carla Braitenberg](https://orcid.org/0000-0001-7277-816X) (2019)

[![GJI paper DOI](https://img.shields.io/static/v1?label=DOI&message=10.1093/gji/ggz344&color=brightgreen)](https://doi.org/10.1093/gji/ggz344)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3358902.svg)](https://doi.org/10.5281/zenodo.3358902)

This repository contains accompanying data sets to the paper
>Pastorutti, A., & Braitenberg, C., 2019. A geothermal application for GOCE satellite gravity data: modelling the crustal heat production and lithospheric temperature field in Central Europe, Geophysical Journal International, doi:10.1093/gji/ggz344
which has been accepted for publication in *Geophysical Journal International* following peer review.

**The version of record is available online at: [doi.org/10.1093/gji/ggz344](https://doi.org/10.1093/gji/ggz344).**

This accompanying data repository has the following DOI: [doi.org/10.5281/zenodo.3358902](https://doi.org/10.5281/zenodo.3358902)

(note that the paper itself is accompanyed by a supplementary material / supporting information pdf, not to be mistaken with this data set)

## `grav` directory

We are providing the grids of gravity data and reductions as plain ASCII files, in (x, y, z) format, comma delimited. The grids cover the extents of area "B", Fig. 2.
All the data is provided in mGal (10<sup>-5</sup> m s<sup>-2</sup>).

* `GGM.xyz` gravity disturbance from `GO_CONS_GCF_2_TIM_r5` GGM, up to d/o 280.
*-* `TOPO.xyz` topographic effect, as gravity disturbance from `RET2014`, up to d/o 280.
* `ISOS.xyz` far-field effect of crustal roots (global), radial component of _g_, synthesised up to d/o 280.
* `SEDS.xyz` sediment infill effect, contrast against 2670 kg/m³. Radial component of _g_, synthesised up to d/o 280.
* `GGMr.xyz` gravity disturbance from the GGM after applying the data reductions.

## `inv` directory

Inverted Moho depth and residuals, provided in the same format as gravity data. For plotting purposes, these grids were un-projected from a UTM35N, and a minimum rectangular bounding box was kept. Due to convergence of meridians, they are slighty smaller than area "B" and area "C".
The thermal model includes the full UTM35N extent of the inverted Moho depth (area "C").

* `MOHO.xyz` Moho depth, km.
* `MOHOres.xyz` Gravity residuals of inverted Moho, mGal.

## `thermal` directory

Input and output of the thermal solver and iterative RHP fitting procedure. Data is provided as .mat files. They can be accessed natively in Matlab or Octave, or in Python through [scipy.io.loadmat](https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.loadmat.html). Default behaviour of loadmat will output Matlab structs as numpy record arrays.

* `Tgrid.mat`, grid coordinates for the thermal model. It contains a struct with the following fields:

  * `x` - vector of UTM35N easting of nodes, step is constant
  * `y` - vector of UTM35N northing of nodes, step is constant 
  * `z` - vector of z coordinates of nodes, positive downwards (this means that topography above msl is negative). Step is coarsening with depth.
  * `UTMstruct` - used map projection (UTM35N), in the form of a map projection structure (see Matlab `defaultm` [documentation](https://mathworks.com/help/map/ref/defaultm.html))
  * `Extents` - vertices of a (x, y) rectangle defining the thermal model area. Vertices are given starting from SW vertex, clockwise.

  The three (x, y, z) coordinate vectors define a rectangular 3D volume, with constant sampling along the horizontal direction and variable sampling along depth. When the volume is visualised this must be taken into account accordingly (by providing z coordinates or by resampling to a constant step).

* `TOUT_VOL.mat`, model volume at the last iteration. Nodes are arranged in a (z, x, y) order. It is a struct named `TOUT_VOL` containing the following 3D arrays, as fields:

  * `L` - layer ID (1 "air" above topography, 2 sediments, 3 upper crust, 4 lower crust, 5 SCLM, 6 "asthenosphere" below SCLM)
  * `A` - radioactive heat production, per unit of volume, in 10<sup>-3</sup> W m<sup>-1</sup>. A is set to 99 under the LAB and 
  * `k` - thermal conductivity, in W m<sup>-1</sup> K<sup>-1</sup>
  * `P` - lithostatic pressure, in Pa
  * `Rho` - density, in kg m<sup>-3</sup>

* `TOUT_iter.mat`, model derived quantities trough the RHP fitting iterations. It is a struct named `TOUT_iter` containing the following fields:

  * `Q0` - surface heat flow, in 10<sup>-3</sup> W m<sup>-1</sup>
  * `Qm` - basal heat flow (at Moho), in 10<sup>-3</sup> W m<sup>-1</sup>
  * `Qs` - heat flow at the sediment to crystalline crust basement, in 10<sup>-3</sup> W m<sup>-1</sup>
  * `misfit` - surface heat flow difference between forward model (Q0) and measurements (where available, NaN everywhere else)
  * `UCA` - Upper Crust RHP used in this iteration, in 10<sup>-3</sup> W m<sup>-1</sup>
  * `LCA` - Lower Crust RHP used in this iteration, in 10<sup>-3</sup> W m<sup>-1</sup>
  * `kEQ` - Series thermal conductivity along each LAB-to-surface column, in W m<sup>-1</sup> K<sup>-1</sup>

  All the `TOUT_iter.mat` fields are 1-by-7 cell arrays. Each cell therein represents one RHP fitting iteration (including the starting conditions, in the first cell). Each cell is a (x, y) array.

* `ProcessedHF.mat`, 2D array (x, y) of filtered surface heat flow measurements on the model grid. NaN where no data is available.

* `TIN_layers.mat`, input data used to build the thermal solver 3D volume. A struct containing:
  * `AIR`, `SEDS`, `UCRUST`, `LCRUST`, `LID`, `AST` - a field for each layer (AIR and AST define the top and bottom solver boundaries)
  * `DefGrid` - grid parameters (coordinates, steps) for the aforementioned layers. This grid, compared to the one provided in `Tgrid`, includes a two-cell wide edge padding, on each side.

## `src` directory

We provide two Matlab/Octave functions to ease the exploration of the thermal volume:

* `TOUT_section.m`
* `TOUT_sectionCall.m`

The `TOUT_sectionCall` function is a tool to interactively plot sections (vertical slices) of the volumes, and of two columns in each section (i.e. depth-wise plots of temperature, thermal conductivity, heat production, heat flow).
It acts as a wrapper for `TOUT_section`, which can be called directly.

`TOUT_sectionCall` must be called without arguments.
In depth usage documentation is provided in each .m file.

First figure: map view. The figure waits for two clicks, defining the section path.

![Map view in figure 1](./images/sectionCall_map.png)

Second figure: section. Same format as section provided in the manuscript. The figure waits for two clicks, defining the along-section position of columns.

![Section view in figure 2](./images/sectionCall_section.png)

Third figure: columns.

![Columns view in figure 3](./images/sectionCall_columns.png)

### Dependencies
Due to a call to [`minvtran`](https://mathworks.com/help/map/ref/minvtran.html), the MATLAB Mapping Toolbox is required.
`minvtran` has not been implemented yet in OCTAVE (see [here](https://wiki.octave.org/Mapping_package#Missing_functions)).
To overcome this limitation, `TOUT_section_NoTopo.m` is provided. It is a version of `TOUT_section.m` stripped of the topography profile extraction; therefore, no calls to `minvtran` are performed.

## `topo` directory

`ETOPO1_005d_crop.mat` is a crop of ETOPO1 ([doi:10.7289/V5C8276M](http://dx.doi.org/10.7289/V5C8276M)) in the study area. It is used in `TOUT_sectionCall` to provide a geographical background for the area map.
