function DEM = sqrt(DEM)

if isinteger(DEM.Z) || islogical(DEM.Z);
    DEM.Z = sqrt(single(DEM.Z));
else
    DEM.Z = sqrt(DEM.Z);
end