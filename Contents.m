% TopoToolbox
% Version 2.4  14-June-2022
%
% TopoToolbox provides a set of Matlab functions that support the analysis
% of relief and flow pathways in digital elevation models. The major 
% aim of TopoToolbox is to offer helpful analytical GIS utilities in a 
% non-GIS environment in order to support the simultaneous application
% of GIS-specific and other quantitative methods.
%
% If you have any questions or remarks, please contact the authors:
%
% Wolfgang Schwanghart
% w.schwanghart[at]geo.uni-potsdam.de
%
% Dirk Scherler
% scherler[at]gfz-potsdam.de
% 
% When you use TopoToolbox in your work, please reference one of these 
% publications:
% 
% Schwanghart, W., Scherler, D. (2014): TopoToolbox 2 – MATLAB-based 
% software for topographic analysis and modeling in Earth surface sciences. 
% Earth Surface Dynamics, 2, 1-7. [DOI: 10.5194/esurf-2-1-2014]
% 
% Schwanghart, W., Kuhn, N. J. (2010): TopoToolbox: a set of Matlab 
% functions for topographic analysis. Environmental Modelling & Software, 
% 25, 770-781. [DOI: 10.1016/j.envsoft.2009.12.002]
%
%
% Objects
%
%     GRIDobj         - object for gridded, geospatial data
%     FLOWobj         - object for flow direction
%     STREAMobj       - object for stream (channel) networks
%     DIVIDEobj       - object for drainage divide networks
%     SWATHobj        - object for swath profiles
%     PPS             - object for point patterns on stream networks
%
% Graphical user interfaces
%
%     flowpathapp     - Map, visualize and export flowpaths that start at manually set channelheads                     
%     slopeareatool   - Interactively create slope area plots and fit power laws                    
%     mappingapp      - Map points combining river planform and profile view
%
% GRIDobj methods
%
%     GRIDOBJ           : Create instance of a GRIDobj
%     GRIDOBJ2ASCII     : write/export GRIDobj to ESRI ArcGIS ASCII file
%     GRIDOBJ2GEOTIFF   : Exports an instance of GRIDobj to a geotiff file
%     GRIDOBJ2IM        : GRIDOBJ2IM Create image from GRIDobj
%     GRIDOBJ2MAT       : convert GRIDobj to matrix and coordinate vectors
%     GRIDOBJ2PM        : combine several GRIDobj into a predictor matrix
%     GRIDOBJ2POLYGON   : Conversion from drainage basin grid to polygon or polyline
%     ACV               : Anisotropic coefficient of variation (ACV) 
%     AGGREGATE         : resampling a GRIDobj using aggregation
%     ARCSLOPE          : mean gradient from a digital elevation model sensu ArcGIS
%     ASPECT            : angle of exposition from a digital elevation model (GRIDobj)
%     CASTSHADOW        : cast shadow
%     CELLAREA          : calculate cell areas of a GRIDobj in geographic coordinate system
%     CLIP              : clip a GRIDobj with a polygon or another GRIDobj
%     CONTOUR           : contour plot of an instance of GRIDobj
%     COORD2IND         : convert x and y coordinates to linear index
%     COORD2SUB         : convert x and y coordinates to subscripts into a GRIDobj
%     CREATEMASK        : create a binary mask using polygon mapping
%     CROP              : crop an instance of GRIDobj with axis-aligned minimum bounding box
%     CURVATURE         : 8-connected neighborhood curvature of a digital elevation model 
%     DEMAREA           : Calculate the corrected surface area of a DEM
%     DEMPROFILE        : get profile along path
%     DILATE            : morphological dilation
%     DIST2CURVE        : labels pixels in a GRIDobj by their directed distance to a curved line
%     DIST2LINE         : labels pixels in a GRIDobj by their distance to a straight line
%     DISTANCE          : distance transform
%     ELEVATEMINIMA     : elevate regional minima in a DEM to their lowest neighbor
%     ERODE             : morphological erosion
%     EXCESSTOPOGRAPHY  : reconstruct surface with threshold-slope surface
%     FILLSINKS         : fill/remove pits, sinks or topographic depressions
%     FILTER            : 2D-filtering of DEMs with different kernels 
%     FIND              : Find indices of nonzero elements in GRIDobj
%     FINDCOORD         : Find coordinates of nonzero elements in GRIDobj
%     GETCOORDINATES    : get coordinate vectors of an instance of GRIDobj
%     GETEXTENT         : return extent of a GRIDobj
%     GETOUTLINE        : get or plot extent of GRIDobj
%     GRADIENT8         : 8-connected neighborhood gradient of a digital elevation model
%     GRIDDEDCONTOUR    : plot contours on grid
%     HILLSHADE         : create hillshading from a digital elevation model (GRIDobj)
%     HYPSCURVE         : plot hypsometric curve of a digital elevation model
%     IDENTIFYFLATS     : identify flat terrain in a digital elevation model
%     IMAGESC           : Scale data in GRIDobj and display as image object
%     IMAGESCHS         : plot hillshade image with overlay
%     IND2COORD         : convert linear index to x and y coordinates
%     INFO              : detailed information on GRIDobj instance
%     INPAINTNANS       : interpolate missing values in a grid (GRIDobj)
%     INTERP            : interpolate to query locations
%     INTERP2GRIDOBJ    : Interpolate scattered data to GRIDobj
%     ISNAN             : returns array elements that are NaNs as logical grid
%     KSDENSITY         : kernel density estimator for GRIDobj
%     LATLON            : returns grid lines of latitudes and longitudes
%     LINE2GRIDOBJ      : convert line to a grid
%     LOCALTOPOGRAPHY   : Local topography
%     MEASURE           : take interactive measurements along a polyline
%     MINMAXNORM        : min-max normalization with optional percent clipping
%     MPOWER            : overloaded power for GRIDobj
%     MRDIVIDE          : overloaded right division for GRIDobj
%     MTIMES            : overloaded multiplication for GRIDobj
%     PAD               : add or remove a border of pixels around a GRIDobj
%     POLYGON2GRIDOBJ   : convert polygon to a grid
%     POSTPROCFLATS     : postprocess flat terrain for visualization purpose
%     PRCCLIP           : percentile clipping
%     PROJECT           : transforms a GRIDobj between projected coordinate systems
%     RECLABEL          : labels GRIDobj by rectangular fields
%     RECLASSIFY        : generate univariate class intervals for an instance of GRIDobj
%     REPROJECT2UTM     : Reproject DEM with WGS84 coordinate system to UTM-WGS84 
%     RESAMPLE          : change spatial resolution of a GRIDobj
%     ROUGHNESS         : terrain ruggedness, position and roughness indices of DEMs
%     SHUFFLELABEL      : shufflelabel randomly relabels a label matrix
%     SNAP2STREAM       : snap gauges or pour points to stream raster
%     SUB2COORD         : convert subscripts to x and y coordinates
%     SURF              : surface plot for GRIDobj
%     TANAKACONTOUR     : Relief depiction using Tanaka contours
%     TOPOSHIELDING     : topographic shielding from cosmic rays
%     VALIDATEALIGNMENT : validates whether instances of GRIDobj are spatially aligned
%     ZSCORE            : standardized z-scores for GRIDobj
% 
% FLOWobj methods
% 
%     FLOWOBJ             : create flow direction object
%     FLOWOBJ2GRIDOBJ     : create ESRI ArcGIS flow direction grid from FLOWobj
%     FLOWOBJ2M           : convert instance of FLOWobj to flow direction matrix 
%     FLOWOBJ2CELL        : return cell array of FLOWobjs for individual drainage basins
%     CROP                : crop an instance of FLOWobj
%     DBENTROPY           : entropy of drainage basin delineation
%     DEPENDENCEMAP       : upslope area for specific locations in a DEM
%     DRAINAGEBASINS      : drainage basin delineation/catchments
%     DRAINAGEBASINSTATS  : zonal statistics on drainage basins
%     FIND                : find indices and values of edges in the flow direction graph
%     FLIPDIR             : Flip direction of flow
%     FLOWACC             : flow accumulation (upslope area, contributing area)
%     FLOWCONVERGENCE     : compute flow convergence of a digital elevation model
%     FLOWDISTANCE        : flow distance in upstream and downstream direction
%     FLOWPATHEXTRACT     : extract linear indices of a single flowpath in a DEM
%     FLOWVEC             : velocity vectors from FLOWobj
%     GRADIENT            : gradient along flow direction
%     IMPOSEMIN           : minima imposition (carving) along drainage network
%     IND2COORD           : convert linear index to x and y coordinates
%     INFLUENCEMAP        : downslope area for specific locations in a digital elevation model
%     ISMULTI             : check if FD is multi or single flow direction
%     MAPFROMNAL          : map values from node-attribute list to nearest upstream grid
%     MULTI2SINGLE        : converts multiple to single flow direction
%     MULTI_NORMALIZE     : create flow direction object
%     MULTI_WEIGHTS       : create flow direction object
%     QUANTCARVE          : quantile carving
%     RANDOMIZE           : randomize multiple flow directions
%     SAVEOBJ             : create flow direction object
%     STREAMORDER         : calculates a stream order GRIDobj from FLOWobj
%     STREAMPOI           : stream points of interest
%     UPDATETOPOSORT      : update topological sorting
%     UPSLOPESTATS        : upslope statistics of a variable based on the flow direction matrix
%     VALIDATEALIGNMENT   : validates whether instances of FLOWobj and GRIDobj are spatially aligned
%     VERTDISTANCE2STREAM : vertical distance to streams
% 
% STREAMobj methods
% 
%     STREAMOBJ           : Create stream object (STREAMobj)
%     STREAMOBJ2GRIDOBJ   : convert an instance of STREAMobj to an instance of GRIDobj
%     STREAMOBJ2SWATHOBJ  : Create swath profile (SWATHobj) from stream network
%     STREAMOBJ2XY        : convert instance of STREAMobj to NaN-separated X and Y coordinates
%     STREAMOBJ2CELL      : convert instance of STREAMobj to cell array of stream objects
%     STREAMOBJ2LATLON    : convert instance of STREAMobj to NaN-separated geographic coordinates
%     STREAMOBJ2MAPSTRUCT : convert instance of STREAMobj to mapstruct
%     STREAMOBJ2SHAPE     : Convert STREAMobj to geoshape or mapshape
%     AGGREGATE           : aggregating values of reaches
%     CHIPLOT             : CHI analysis for bedrock river analysis
%     CHITRANSFORM        : Coordinate transformation using the integral approach
%     CLEAN               : Create stream object (STREAMobj)
%     CONNCOMPS           : labels of connected components (individual trees) in a stream network
%     CRS                 : constrained regularized smoothing of the channel length profile
%     CRSAPP              : interactive smoothing of river long profiles
%     CRSLIN              : constrained regularized smoothing of the channel length profile
%     CUMMAXUPSTREAM      : cumulative maximum in upstream direction
%     CUMSUM              : cumulative sum on stream network
%     CUMTRAPZ            : Cumulative trapezoidal numerical integration along a stream network
%     CURVATURE           : curvature or 2nd derivative of a STREAMobj
%     DENSIFY             : Increase number of vertices in stream network using splines
%     DISTANCE            : return node attribute list with distances along the stream network
%     DRAINAGEDENSITY     : drainage density of a stream network
%     EXTRACTCONNCOMPS    : interactive stream network selection
%     GETNAL              : get node-attribute list
%     GETVALUE            : retrieve value from node-attribute list
%     GRADIENT            : along-stream gradient
%     HILLSLOPEAREA       : upslope hillslope area for each stream pixel 
%     IDENTIFYFLATS       : identify flat sections in a river profile
%     IMPOSEMIN           : minima imposition (carving) along stream network
%     INFO                : meta information about STREAMobj
%     INPAINTNANS         : inpaint missing values (nans) in a node attribute list
%     INTERP              : interpolate data on STREAMobj (single river only)
%     INTERSECT           : intersect different instances of STREAMobj 
%     INTERSECTLOCS       : Derive locations where two STREAMobj start to have a common network
%     ISGEOGRAPHIC        : Determines if STREAMobj S has a geographic coordinate system
%     ISNAL               : test whether a vector is a node attribute list of a STREAMobj
%     ISSUBGRAPH          : Create stream object (STREAMobj)
%     KLARGESTCONNCOMPS   : retain k largest connected components in an instance of STREAMobj
%     KNICKPOINTFINDER    : find knickpoints in river profiles
%     KSN                 : normalized steepness index
%     LABELREACH          : create node-attribute list with labelled reaches
%     MAPLATERAL          : map values of regions adjacent to streams to stream network
%     MCHI                : gradient of stream profile in chi space (M_chi)
%     MEANUPSTREAM        : mean (weighted) upstream  values
%     MINCOSTHYDROCON     : minimum cost hydrological conditioning
%     MNOPTIM             : Bayesian optimization of the mn ratio
%     MODIFY              : modify instance of STREAMobj to meet user-defined criteria
%     NAL2NAL             : map one node-attribute list to another 
%     NETDIST             : distance transform on a stream network
%     NETWORKSEGMENT      : Identify river segments and compute segment geometry
%     ORIENTATION         : stream orientation
%     PLOT                : plot instance of STREAMobj
%     PLOT3               : 3d-line plot of a STREAMobj
%     PLOT3D              : 3D plot of a stream network
%     PLOTC               : plot a colored stream network
%     PLOTDZ              : plot upstream distance version elevation of a stream network
%     PLOTDZSHADED        : plot upstream distance version elevation of a stream network
%     PLOTSEGMENTGEOMETRY : Plot segment geometry obtained from the function networksegment
%     PLOTSTREAMORDER     : calculate stream order from STREAMobj
%     PSITRANSFORM        : Parse Inputs
%     QUANTCARVE          : quantile carving
%     RANDLOCS            : Random locations along the stream network
%     REMOVEEDGEEFFECTS   : remove potential edge effects
%     REMOVESHORTSTREAMS  : Remove first order streams with a length less than specified
%     RMEDGE              : Create stream object (STREAMobj)
%     RMNODE              : Create stream object (STREAMobj)
%     SIDEBRANCHING       : side branching classification according to Tokunaga (1978)
%     SINUOSITY           : sinuosity coefficient 
%     SLOPEAREA           : slope-area relation of a stream network
%     SMOOTH              : smoothing of node-attribute lists
%     SNAP2STREAM         : snap locations to nearest stream location
%     SPLIT               : split drainage network at predefined locations
%     STACKEDPLOTDZ       : plot several stacked variables along the stream networks
%     STREAMORDER         : calculate Strahler Stream Order from STREAMobj
%     STREAMPOI           : stream points of interest
%     STREAMPROJ          : project stream elevations based on slope-area scaling
%     SUBGRAPH            : Create stream object (STREAMobj)
%     TRANSFORMCOORDS     : transform coordinates of stream network
%     TRIBDIR             : direction of inflow of tributary
%     TRUNK               : extract trunk stream (longest stream) 
%     UNION               : merge different instances of STREAMobj into a new instance
%     VALIDATEALIGNMENT   : is an instance of STREAMobj is spatially aligned with another object of TopoToolbox
%     WIDENSTREAM         : level elevations adjacent to the stream network
%     WMPLOT              : plot stream network in the webmap browser
%     ZEROBASELEVEL       : set base level to zero
%
% DIVIDEobj methods
%
%     DIVIDEOBJ           : Create divide object (DIVIDEobj)
%     DIVIDEOBJ2MAPSTRUCT : obtain divide properties from GRIDobj
%     ASYMMETRY           : directional asymmetry of divide segments
%     CLEANEDGES          : Delete divides on the edges of DEM
%     COORD2IND           : convert x and y coordinates to linear index
%     DIST2NODE           : network distance to nodes
%     DIVDIST             : DIVDIST   Assign distance to divide segments
%     DIVNET              : Create divide object (DIVIDEobj)
%     DIVORDER            : Assign order to divide segments
%     GETVALUE            : get grid values adjacent to divides
%     IND2COORD           : convert linear index to x and y coordinates
%     JCTANGLE            : angles between divide segments at junctions
%     JCTCON              : compute junction connectivity
%     PLOT                : plot the divide network 
%     PLOTC               : Create colored plot of drainage divide network
%     REMOVESHORTDIVIDES  : Remove short first-order divides
%     SHRINK              : Shrink divide network
%     SORT                : Sort divide segments by network structure.
% 
% SWATHobj methods
% 
%     SWATHOBJ2GRIDOBJ : create a GRIDobj with swath-specific information
%     SWATHOBJ2GDS     : create a geographic data structure from a SWATHobj
%     CONVERT2LATLON   : convert spatial fields in SWATHobj to lat,lon
%     MAPSWATH         : obtain Z values along a SWATHobj from an arbitrary GRIDobj
%     PLOT             : plot instance of SWATHobj
%     PLOTDZ           : create distance-elevation plot of SWATHobj
%     PLOTDZM          : create color-coded distance-elevation plot from SWATHobj and GRIDobj
%     PROFILES         : obtain profiles from a SWATHobj
%     TIDY             : remove overlapping points from SWATHobj
%
% PPS methods
%
%     PPS              : Point patterns on stream networks
%     AGGREGATE        : Aggregate points in PPS to new point pattern
%     AS               : Convert PPS object into various data formats 
%     CLUSTER          : hierarchical clustering of points in PPS
%     CONVHULL         : convex hull around points in PPS
%     DENSITY          : nonparametric density estimation on network
%     ECDF             : Empirical cumulative distribution function of a covariate 
%     EXTENDEDNETWORK  : extend network to account for duplicate points
%     FITLOGLINEAR     : fit loglinear model to point pattern
%     GETMARKS         : extract point marks
%     GFUN             : G-function (nearest inter-point distance distribution)
%     HASDUPLICATES    : checks whether there are duplicate points
%     HISTOGRAM        : histogram of point pattern on stream network
%     INTENSITY        : calculate intensity (density) of points on the stream network
%     NETDIST          : Shortest network distance
%     NPOINTS          : number of points in the point pattern
%     PLOT             : plot instance of PPS
%     PLOTC            : plot a colored stream network
%     PLOTDZ           : plot upstream distance version elevation or covariate of a PPS
%     PLOTEFFECTS      : Plot of slices through fitted generalized linear regression
%     PLOTPOINTS       : plot points of PPS
%     POINTDISTANCES   : inter-point distance
%     POINTS           : extract a list of points from the point pattern
%     QUADRATCOUNT     : Quadrat count and chi2 test
%     RANDOM           : random realisation of a (inhomogeneous) point process
%     REGULARPOINTS    : Generate non-random points on stream network
%     REMOVEDUPLICATES : removes duplicate points
%     REMOVEPOINTS     : Remove points in point pattern
%     RHOHAT           : nonparametric estimation of point pattern dependence on covariate
%     ROC              : receiver-operating characteristics of point pattern
%     SHAPEWRITE       : write point pattern to shapefile
%     SIMULATE         : simulate point pattern using random thinning
%     TLENGTH          : total length of the stream network
%     VORONOI          : nearest neighbor search on a stream network
%     WMPLOT           : plot instance of PPS in webmap browser
%
% TTLEM (TopoToolbox Landscape Evolution Model)
%
%     TTLEM            : TopoToolbox Landscape Evolution Model: TTLEM
%     TTLEMSET         : set options and parameters for TTLEM
%
% Other tools (utitilies)
%
%     compilemexfiles - compile mex-functions that come with TopoToolbox 2
%     coord2ind       - convert xy coordinates to linear index
%     exaggerate      - elevation exaggeration in a 3D surface plot
%     getdistance     - cumulative distance along path defined by vertice coordinates
%     getextent       - get current axis extent
%     interpline      - distribute vertices evenly along line with irregularly spaced vertices
%     ixneighbors     - neighbor indexing in a matrix
%     label2poly      - plot region outlines with polyline
%     landcolor       - colormap for the display of DEMs
%     setextent       - set current axis extent    
%     showmethods     - displays class method names and H1 lines in the command line
%     egm96heights    - read and resample EGM96 geoid heights
%     dpsimplify      - Douglas-Peucker line simplification
%     zonalstats      - Zonal statistics
%
% Documentation
%
%   usersguide_1_intro      - introductory script to the toolbox
%   usersguide_2_mfd        - how to calculate Multiple Flow Direction in
%                             TopoToolbox 2
%   usersguide_3_ksn        - calculate ksn values with TopoToolbox 2
%
%
%   ---------------------------------------------------------------------
%     TopoToolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

