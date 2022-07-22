# TopoToolbox - a set of Matlab functions for topographic analysis

<img src="https://github.com/wschwanghart/topotoolbox/blob/master/topotoolbox.jpg" align="center" height="100">

[![View TopoToolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://de.mathworks.com/matlabcentral/fileexchange/50124-topotoolbox)

[TopoToolbox](http://topotoolbox.wordpress.com) provides a set of Matlab functions that support the analysis
of relief and flow pathways in digital elevation models. The major 
aim of TopoToolbox is to offer helpful analytical GIS utilities in a 
non-GIS environment in order to support the simultaneous application
of GIS-specific and other quantitative methods.

If you have any questions or remarks, please contact the authors:

[Wolfgang Schwanghart](https://www.uni-potsdam.de/de/umwelt/institut/alle-mitarbeiterinnen/schwanghart-wolfgang.html)
w.schwanghart[at]geo.uni-potsdam.de

[Dirk Scherler](https://sites.google.com/site/scherlerdirk/home)
scherler[at]gfz-potsdam.de

## Requirements

TopoToolbox is plat-form independent and requires
Matlab 2016b or higher and the Image Processing Toolbox. The Mapping
Toolbox is not mandatory, but good to have to facilitate easy data
exchange with GIS software. Some functions support parallelisation using the 
Parallel Toolbox. Few functions require the Optimization or Statistics
and Machine Learning Toolbox.

## References

When you use TopoToolbox in your work, please reference following 
publication:

- Schwanghart, W., Scherler, D. (2014): TopoToolbox 2 – MATLAB-based 
software for topographic analysis and modeling in Earth surface sciences. 
Earth Surface Dynamics, 2, 1-7. DOI: [10.5194/esurf-2-1-2014](http://dx.doi.org/10.5194/esurf-2-1-2014)

If you are using version 1, then please refer to this publication:

- Schwanghart, W., Kuhn, N.J. (2010): TopoToolbox: a set of MATLAB 
functions for topographic analysis. Environmental Modelling & Software, 
25, 770-781. DOI: [10.1016/j.envsoft.2009.12.002](http://dx.doi.org/10.1016/j.envsoft.2009.12.002)

In addition, various models and algorithms used in TopoToolbox have been published in the following articles. 

### DEM preprocessing and carving

- Schwanghart, W., Groom, G.B., Kuhn, N.J., Heckrath, G., 2013: Flow network derivation from a high 
resolution DEM in a low relief, agrarian landscape. Earth Surface Processes and Landforms, 38, 
1576-1586. DOI: [10.1002/esp.3452](http://dx.doi.org/10.1002/esp.3452)

- Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: uncertainty assessment and smoothing 
using quantile regression techniques. Earth Surface Dynamics, 5, 821-839. DOI: [10.5194/esurf-5-821-2017](https://doi.org/10.5194/esurf-5-821-2017)

### TopoToolbox Landscape Evolution Modelling (TTLEM) and HyLands

- Campforts, B., Schwanghart, W., Govers, G. (2017): Accurate simulation of transient 
landscape evolution by eliminating numerical diffusion: the TTLEM 1.0 model. 
Earth Surface Dynamics, 5, 47-66. DOI: [10.5194/esurf-5-47-2017](http://dx.doi.org/10.5194/esurf-5-47-2017)

- HyLands: Campforts B., Shobe M.C., et al. (2020): HyLands 1.0: a Hybrid Landscape 
evolution model to simulate the impact of landslides and landslide-derived sediment on landscape evolution. 
Geosci. Model Dev., 13, 3863–3886. DOI: [10.5194/gmd-13-3863-2020](http://dx.doi.org/10.5194/gmd-13-3863-2020)

### Excess topography

- Blöthe, J.H., Korup, O., Schwanghart, W., 2015: Large landslides lie low: Excess topography in the Himalaya-Karakorum ranges. 
  Geology, 43, 523-526. [DOI: 10.1130/G36527.1](https://doi.org/10.1130/G36527.1)

### Knickpointfinder

- Stolle, A., Schwanghart, W., Andermann, C., Bernhardt, A., Fort, M., Jansen, J.D., Wittmann, H., Merchel, S., 
  Rugel, G., Adhikari, B.R., Korup, O., 2019. Protracted river response to medieval earthquakes. Earth Surface Processes 
  and Landforms, 44, 331-341. DOI: [10.1002/esp.4517](https://doi.org/10.1002/esp.4517) (The description here is very terse, yet)

### Divide functions

- Scherler, D., Schwanghart, W., 2020. Drainage divide networks – Part 1: Identification and ordering in digital elevation models. 
Earth Surface Dynamics, 8, 245–259. [DOI: 10.5194/esurf-8-245-2020](http://dx.doi.org/10.5194/esurf-8-245-2020)

- Scherler, D., Schwanghart, W., 2020. Drainage divide networks – Part 2: Response to perturbations. 
Earth Surface Dynamics, 8, 261-274. [DOI: 10.5194/esurf-8-261-2020](http://dx.doi.org/10.5194/esurf-8-261-2020)

### Point patterns on stream networks

- Schwanghart, W., Molkenthin, C., & Scherler, D. (2020). A systematic approach and software for the analysis 
of point patterns on river networks. Earth Surface Processes and Landforms, accepted. [DOI: 10.1002/esp.5127](http://dx.doi.org/10.1002/esp.5127)

## Getting started

Before working with TopoToolbox the directories and functions must be on the search 
path of Matlab. Enter following code into the command line:

        addpath(genpath('C:\path\to\wherever\you\installed\this\TopoToolbox-2'))
		
To remove .git-folders from the path, run

        rmpath(genpath('C:\path\to\wherever\you\installed\this\TopoToolbox-2\.git'));
		
To make these paths permanent, use the command

		savepath

Type `doc` in the command line to open the main documentation page. You'll 
find the TopoToolbox documentation in the section Supplemental Software.  
The documentation contains several user's guides that will help you 
getting started. In addition, TopoToolbox functions have extensive help 
sections (e.g. `help gradient8` or `help STREAMobj/modify`. An additional 
resource for code and examples is the [TopoToolbox blog](http://topotoolbox.wordpress.com).


***
## Version History

### pre 2.5

- new function: GRIDobj/randomsample - Create a spatially uniform sample based on a GRIDobj
- new function: tthelp - Search for keywords on the TopoToolbox blog

### 2.4

- HyLands 1.0 added: see paper: Campforts, B. et al. (2020): HyLands 1.0: a Hybrid Landscape evolution 
  model to simulate the impact of landslides and landslide-derived sediment on landscape evolution  
  Geoscientific Model Development. [DOI:10.5194/gmd-13-3863-2020](https://doi.org/10.5194/gmd-13-3863-2020)
- new class: DIVIDEobj [Paper 1 DOI: 10.5194/esurf-8-245-2020](https://doi.org/10.5194/esurf-8-245-2020)
  [Paper 2 DOI: 10.5194/esurf-8-261-2020](https://doi.org/10.5194/esurf-8-261-2020)
- new class: PPS [DOI: 10.1002/esp.5127](http://dx.doi.org/10.1002/esp.5127)
- modification: update to ttscm to Scientific Colormaps 7.0
  see [Fabio Crameri's website](http://www.fabiocrameri.ch/colourmaps.php)
- new function: FLOWobj/plotdbfringe
- new function: FLOWobj/tfactor
- new function: private in FLOWobj graydistparallel
- new function: GRIDobj/diffusion
- new function: GRIDobj/histogram
- new function: GRIDobj/rand
- new function: PPS/extractvaluesaroundpoints
- new function: STREAMobj/extend2divide
- new function: STREAMobj/binarize
- modification: STREAMobj/modify 
- modification: STREAMobj/wmplot
- new function: STREAMobj/loessksn
- new function: STREAMobj/STREAMobj2shape
- new function: STREAMobj/isequal
- new function: STREAMobj/isempty
- new function: STREAMobj/STREAMobj2kml
- new function: STREAMobj/getlocation
- new function: STREAMobj/fastscape
- modification: STREAMobj/smooth
- new function: STREAMobj/mnoptimvar
- modification: STREAMobj/STREAMobj2cell
- modification: STREAMobj/netdist and PPS/netdist
- new function: STREAMobj/istrunk
- modification: STREAMobj/trunk - added second output argument
- modification: GRIDobj/demprofile
- new function: ScaleBar
- new function: xlinerel and ylinerel 
- new function: gif by Chad Greene [see here](https://de.mathworks.com/matlabcentral/fileexchange/63239-gif)
- topoapp was removed from this version 
- updates to readopentopo
- updates to utilities: setextent and getextent 
- new function: padextent
- new function: ukrainecolor


### 2.3

- Documentation in the documentation browser
- new function: ttcmap for creating nice colormaps for DEMs, particularly
  if DEMs include topography and bathymetry
- FLOWobj2gradient renamed to gradient
- new function: ttscm for access to scientific colormaps;
  see [Fabio Crameri's website](http://www.fabiocrameri.ch/colourmaps.php)
- new function: mappingapp (lightweighed GUI for mapping points simultaneously in
                planform and profile view) still beta!
- new function: FLOWobj/mapfromnal
- enhancement:  FLOWobj/multi2single allows area thresholding
- new function: STREAMobj/getvalue
- new function: STREAMobj/hillslopearea
- new function: STREAMobj/zerobaselevel
- new function: STREAMobj/knickpointfinder
- new function: STREAMobj/stackedplotdz
- new function: STREAMobj/sinuosity
- new function: STREAMobj/clean
- new function: STREAMobj/nal2nal
- new function: STREAMobj/netdist
- new function: STREAMobj/tribdir 
- new function: GRIDobj/clip
- new function: GRIDobj/GRIDobj2im
- new function: GRIDobj/getextent
- new function: IOtools/readexample
- renamed GRIDobj/project2GRIDobj to GRIDobj/project. In addition, the function has 
  a number of new functionalities.
- modification: STREAMobj/modify has new options for interactively modifying stream networks
  and to extract streams that confluence from a specified direction.
- changes to readopentopo, getcoordinates, getoutline, polygon2GRIDobj, line2GRIDobj
- new function: hydrosheds2FLOWobj (see in IOtools)
- new function: egm96heights including the grid ww15mgh.grd that contains global geoid 
                undulations based on the EGM96 geoid.


### 2.2

- TTLEM is part of TopoToolbox;
  see our paper: Campforts, B., Schwanghart, W., Govers, G. (2017): Accurate simulation 
  of transient landscape evolution by eliminating numerical diffusion: the TTLEM 1.0 model. 
  Earth Surface Dynamics, 5, 47-66. [DOI: 10.5194/esurf-5-47-2017](http://dx.doi.org/10.5194/esurf-5-47-2017)
- new functions for smoothing and hydrological correction: 
  STREAMobj/crs, STREAMobj/crsapp, STREAMobj/smooth, STREAMobj/crslin, 
  STREAMobj/quantcarve, FLOWobj/quantcarve;
  see our paper: Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: uncertainty 
  assessment and smoothing using quantile regression techniques. Earth Surface Dynamics, 5, 
  821-839. [DOI: 10.5194/esurf-5-821-2017](http://dx.doi.org/10.5194/esurf-5-821-2017)
- modification: GRIDobj way to store referencing information was changed
- modification: FLOWobj now supports multiple flow directions and Dinf.
- modification: several new options for imageschs
- new function: FLOWobj/dbentropy 
- new function: FLOWobj/updatetoposort
- new function: GRIDobj/aggregate
- new function: GRIDobj/createmask
- new function: GRIDobj/dist2line
- new function: GRIDobj/dist2curve
- new function: GRIDobj/line2GRIDobj
- new function: GRIDobj/GRIDobj2pm
- new function: GRIDobj/minmaxnorm
- new function: GRIDobj/reclabel
- new function: GRIDobj/zscore
- new function: GRIDobj/pad
- new function: STREAMobj/aggregate
- new function: STREAMobj/labelreach
- new function: STREAMobj/distance
- new function: STREAMobj/drainagedensity
- new function: STREAMobj/plotc
- new function: STREAMobj/plotdzshaded
- new function: STREAMobj/meanupstream
- new function: STREAMobj/plot3
- new function: STREAMobj/chitransform
- new function: STREAMobj/cumtrapz
- modification: STREAMobj/modify includes option rmconncomps
- new function: STREAMobj/mchi
- new function: STREAMobj/conncomps
- modification: STREAMobj/extractconncomps (new GUI and behavior)
- modification: STREAMobj/union (new syntax)
- modification: STREAMobj/intersect (new syntax)
- new function: STREAMobj/isnal
- new function: STREAMobj/info
- new function: STREAMobj/orientation
- new function: STREAMobj/plotstreamorder
- new function: STREAMobj/removeedgeeffects
- new function: STREAMobj/split
- new function: STREAMobj/streamproj
- new function: STREAMobj/networksegment
- new function: STREAMobj/maplateral
- new function: STREAMobj/plotsegmentgeometry
- new function: STREAMobj/randlocs
- new function: STREAMobj/zerobaselevel
- modification: STREAMobj/streamorder plotting option removed
- modification: STREAMobj/plotdz includes custom distance option
- modification: STREAMobj/distance includes option to derive distance from different
  STREAMobj
- modification: STREAMobj/STREAMobj2cell
- modification: STREAMobj/STREAMobj2mapstruct
- new function: STREAMobj/transformcoords
- new function: FLOWobj/FLOWobj2cell
- update to several FLOWobj methods to avoid speed loss for MATLAB versions newer
  than R2015b
- removed bug in GRIDobj/curvature

### 2.1

- new function: GRIDobj/excesstopography; see our paper Blöthe, J.H., Korup, O., 
  Schwanghart, W. (2015): Large landslides lie low: Excess topography in the 
  Himalaya-Karakorum ranges. Geology, 43, 523-526. [DOI: 10.1130/G36527.1](http://dx.doi.org/10.1130/G36527.1)
- new function: GRIDobj/GRIDobj2polygon
- new function: STREAMobj/getnal
- new function: STREAMobj/sidebranching
- new function: STREAMobj/mincosthydrocon
- new function: STREAMobj/intersectlocs
- new function: STREAMobj/densify
- new function: STREAMobj/plot3d
- new function: STREAMobj/widenstream
- new function: demo_modifystreamnet
- modification of the function slopearea
- better performance of FLOWobj when converting from flow direction matrix 
  by using dmperm to perform topological sort  
- new function: GRIDobj/toposhielding
- new function: demarea (incorporation of Juernjakob Dugge's function: 
  http://www.mathworks.com/matlabcentral/fileexchange/42204-calculate-dem-surface-area )
- new function: GRIDobj/getoutline
- removed a bug when some functions such log, log10 were called with integer 
  GRIDobjs
- additional overloading of built-in functions for GRIDobjs. We added matrix
  arithmetics, which, however, perform element-wise operations (e.g. mtimes can
  be used with GRIDobj now, but performs times)
- the scope of the function GRIDobj/localtopography was enhanced (min, max, 
  percentiles, etc in a disk-shaped neighborhood)
- FLOWobj/streampoi and STREAMobj/streampoi now return 'bconfluences', 
  e.g. stream pixels that in downstream direction are located immediately 
  before confluences.
- new function: STREAMobj/imposemin - limits downstream minima imposition
  to the stream network
- several bug fixes
- demo_modifystreamnet.m
- preprocessapp was removed

### V2.0.1

- removed bug in GRIDobj
- removed bug with case-sensitivity in some functions
- removed bug with internal drainage option in FLOWobj

### V2.0

- new functions STREAMobj/intersect, STREAMobj/union
- new interactive tools in STREAMobj/modify
- new interactive tool GRIDobj/measure


### V2.0beta 

V2.0 introduces an object oriented approach towards grid representation, flow direction, stream 
networks and swath objects. The performance of various, inparticular flow related, functions was 
increased. Mex-files have been written to increase the speed of some functions and are delivered as 
64bit Windows binaries. They have been compiled on Windows 7 with an Intel processor, so they should be 
compiled before using them, if your system differs. However, compiling is not mandatory, since m 
versions are available, too, which are a little slower. Please refer to the Contents.m file for a 
complete list of functions.

### V1.06 --- 11. November 2011

- new function: acv
- new function: cropmat
- new function: dbentropy
- new function: deminpaint
- new function: exaggerate
- new function: label2poly
- new function: routegeodesic (as optimal method for routing through flats and 
  pits). Requires Matlab 2011b and will be made the default routing algorithm
  in future releases.
- new function: upslopesidelength
- new function: upslopestats
- new function: getextent and setextent
- function enhancements: rasterread and rasterwrite. 
- function enhancements: roughness
- new users guide on processing flat areas
- the baranja_hill.mat dataset was added. It was obtained from 
  here: http://geomorphometry.org/content/baranja-hill

### V1.05 --- 15. September 2010

- some of the functions now employ the function validateattributes to check
  input arguments. Note that this might return a warning on older versions than
  2009a.
- a bug in routeflats was removed
- new function: M2UV
- new function: twi
- new function: aspect
- new function: ismulti
- new function: postprocflats
- new function: demsobel
- new function: hypscurve
- new function: roughness
- removed bug in flowacc. flowacc returned an error when called with three
  input arguments
- gradient8 allows you to output different angular units, major speed 
  increase when called with gradient as sole output argument.
- functional enhancements to flowdistance (see help flowdistance)
- functional enhancements to identifyflats
- rasterread and rasterwrite now feature dialog boxes for reading and
  saving files if no filenames are supplied to the function
- flowdistanceds can now calculate the maximum downward flowpath distance 
  for each cell at one step. 

### V1.04 --- 5. January 2010 - *first public release*

- a bug in flowacc_lm was removed. When a weight matrix W0
  was supplied as additional input argument, W0 was set to 
  dem.
- new function: imageschs
- minor changes to hillshade were made. The algorithm is now
  based on method proposed by Katzil and Doytsher, 2003.
- flowacc allows for another input argument (runoff ratio). 

### V1.03 --- 5. November 2009

- sub-basin analysis has been added as new functionality 
  (see sbstruct, sbplot, sbprops)
- new function: adjustgauges

### V1.02 --- 30. October 2009  

- major speed enhancement for fillsinks with maximum
  fill depth

### V1.01 --- 16. September 2009  

- hillshade plots the hillshade matrix when no output
  arguments are defined.
- flowacc_lm was optimized, so that large, flat areas
  are handled more memory efficiently.
- a bug in routeflats was removed.
- crossflats was updated to run most efficient on 
  Matlab R2009a.
- drainagebasins functionality was enhanced to allow
  for the delineation of drainage basins of specified 
  order. The former support of multiple flow direction
  has been removed.
- new function: drainagedensity
- new function: shufflelabel
- influencemap got a second output Mstreams

### V1.0 --- 15. March 2009 

- release
