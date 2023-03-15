function DEM = log2(DEM)

if isinteger(DEM.Z) || islogical(DEM.Z);
    DEM.Z = log2(single(DEM.Z));
else
    DEM.Z = log2(DEM.Z);
end