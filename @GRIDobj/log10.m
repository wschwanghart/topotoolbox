function DEM = log10(DEM)


if isinteger(DEM.Z) || islogical(DEM.Z);
    DEM.Z = log10(single(DEM.Z));
else
    DEM.Z = log10(DEM.Z);
end