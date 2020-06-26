# Erosion Inversion ESPL
Inversion of provenance data and sediment load into spatially varying erosion rates

## Introduction
This repository contains 2 Matlab scripts:
1) `Synthetic_tests/Forward_Inverse_ESPL.m` for synthetic tests (forward - inverse tests for a synthetic catchment)
2) `Natural_Example_Marsyandi/INVERSION_Marshyangdi.m` that computes an erosion rate map from detrital data, source area information and sediment load of the river of interest.

The inversion scheme is based on the following forward statement: 
**d** = **G\*edot** 
where **d** is a vector containing the measured tracer distribution at the outlet of the catchment, **G** is the tracer concentration matrix (with a size of n tracers x n pixels) and **edot** is the erosion rate map, with n pixels for which the erosion rate (in mm/y) is known. 
The inversion scheme computes **edot** from **G** and **d**, to deal with the underdeterminedness of the problem, smoothing is imposed in the form of a covariance matrix that takes into account the spatial covariance of the tracer concentration information.

The forward - inverse script can be run in different modes:
1) with different true erosion rate patterns (block, gaussian bump or checkerboard)
2) with or without additional subcatchment data

First, data are computed with the forward model (**d** = **G\*edot**) using a chosen, 'true' erosion rate map and chosen tracer concentrations for the different geological units. These data are then inverted into spatially varying erosion rates with a closed-form, linear least-squares inversion scheme. Posterior uncertainties are evaluated by mapping the spread function, resolution, posterior and reduced variance, and the difference between the true and posterior erosion rate model. See paper *Inversion of provenance data and sediment load into spatially varying erosion rates* (ESPL, in review) for more information and examples of synthetic tests. 

The Marsyangdi script can be run in the following modes:
1) with geological units as source regions
2) with tributary catchments as source areas
3) with additional subcatchment data (possible in both geological and tributary modes)

## How to use
These scripts require Matlab and were developed in Matlab2019a.

### Forward - inversion script
Variables that can be changed:
1) Forward parameters:
   - flags (to set type of erosion rate pattern and for additional data)
   - settings G matrix (number of tracers, number of geological units)
   - number of subcatchment samples for additional data
   - number of analysed ages (needed to construct tracer concentrations)
   
2) Inverse parameters
   - smoothing distance
   - model variance
   - data uncertainty

### Marsyandi script
#### Input data
The following variables must be updated in the script according to the data the user wants to invert.
1) Source area map
  (geological units or tributaries)
  (geotiff file, projected coordinate system `ProjectedCSTypeGeoKey`)
  (`nanval` (e.g. zero) outside of source areas)
  (source areas have ID's that are equal to the sample ID's in the sample data file)
2) Sample data: file with source and detrital data
  (.csv file, first column: sample ID, other columns: different tracers,
   every row corresponds to tracer information for a different source area)
3) Shapefile with sample locations
   (attribute table: ID corresponds to source area e.g. sample 3 is
   representative for source area 3)
   (contains n_datasamples (detrital) samples for which the data will be inverted,
   default = 1)
4) ID of sample w/ data to be inverted
5) Q_s
(sediment discharge for catchment outlet, m^3/y)

#### Inverse parameters
1) prior erosion rate estimation (standard: equal to mean erosion rate derived from Q_s data; constant)
2) standard deviation model covariance
3) smoothing distance
4) uncertainty on data

## Further information
For further information, we refer to the comments in the scripts as well as the ESPL paper *Inversion of provenance data and sediment load into spatially varying erosion rates* of De Doncker F., Herman F. and Fox M. (in review)
