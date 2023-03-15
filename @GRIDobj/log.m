function DEM = log(DEM)

if isinteger(DEM.Z) || islogical(DEM.Z);
    DEM.Z = log(single(DEM.Z));
else
    DEM.Z = log(DEM.Z);
end