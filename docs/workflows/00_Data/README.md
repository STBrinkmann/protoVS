Data Download
================
Sebastian Brinkmann
(31.3.2022)

## Elevation data

For testing the algorithms we used LiDAR derived Digital Surface Models
(DSM) and Digital Terrain Models (DTM) with 1 m spatial resolution of
the Vancouver metropolitan area (Natural Resources Canada 2019). For
downloading and processing the data, follow the workflow provided in
this folder
([`00_elevation_data.R`](https://github.com/STBrinkmann/protoVS/blob/main/docs/workflows/00_Data/00_elevation_data.R)).

## Land Use and Land Cover

Publicly available LULC data has been acquired by Metro Vancouver
(31.11.2019) at 2m resolution and reclassified to a binary greenness
raster (0=no-green; 1=green). On their website
(<http://www.metrovancouver.org/data>) you can downloaded the data by
clicking on “FGDB,” shown in the screenshot below.

[![Land Cover Classification 2014 - 2m LiDAR (Raster) (Metro Vancouver
31.11.2019).](metroVancouver.PNG)](http://www.metrovancouver.org/data)

## Bibliography

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-MetroVancouver.31.11.2019" class="csl-entry">

Metro Vancouver. 31.11.2019. “Land Cover Classification 2014 - 2m LiDAR
(Raster).” <http://www.metrovancouver.org/data>.

</div>

<div id="ref-NaturalResourcesCanada.2019" class="csl-entry">

Natural Resources Canada. 2019. “High Resolution Digital Elevation Model
(HRDEM) - CanElevation Serie.”
<https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995>.

</div>

</div>
