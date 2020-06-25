# Erosion Inversion ESPL
Inversion of provenance data and sediment load into spatially varying erosion rates

## Introduction
This repository contains a Matlab script that computes an erosion rate map from detrital data, source area information and sediment load of the river of interest.
The inversion scheme is based on the following forward statement: 
**d** = **G\*edot** 
where **d** is a vector containing the measured tracer distribution at the outlet of the catchment, **G** is the tracer concentration matrix (with a size of n tracers x n pixels) and **edot** is the erosion rate map, with n pixels for which the erosion rate (in mm/y) is known. 
The inversion scheme computes **edot** from **G** and **d**, to deal with the underdeterminedness of the problem, smoothing is imposed in the form of a covariance matrix that takes into account the spatial covariance of the tracer concentration information.

The script can be run in different modes:
1) with geological units as source regions
2) with tributary catchments as source areas
3) with additional subcatchment data (possible in both geological and tributary modes)

## How to use
This script requires Matlab. It was developed in Matlab2019a.
### Input data
The following variables must be updated in the script according to the data the user wants to invert.
1) Source area map
  (geological units or tributaries)
  (geotiff file, projected coordinate system `ProjectedCSTypeGeoKey`)
  (`nanval` (e.g. zero) outside of source areas)
  (source areas have ID's that are equal to the sample ID's in the sample data file)
2) Sample data: file with source and detrital data
  (.csv file, first column: sample ID, other columns: different tracers,
   every row corresponds to tracer information for a different source area)
3) Shapefile with sample locations        --->(source_samples)
   (attribute table: ID corresponds to source area e.g. sample 3 is
   representative for source area 3)
   (contains n_datasamples (detrital) samples for which the data will be inverted,
   default = 1)
4) ID of sample w/ data to be inverted
5) Q_s
(sediment discharge for catchment outlet, m^3/y)

### Inverse parameters
1) prior erosion rate estimation (standard: equal to mean erosion rate derived from Q_s data; constant)
2) standard deviation model covariance
3) smoothing distance
4) uncertainty on data

## Further information
For further information, we refer to the comments in the script as well as the ESPL paper *Inversion of provenance data and sediment load into spatially varying erosion rates* of De Doncker F., Herman F. and Fox M. (in review)
