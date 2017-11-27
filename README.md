# ArcGIS-Protein-Analysis-Toolbar
Python toolbox "Protein" for use in ArcCatalog for the protein concentration gradient
and protein body size-distribution analysis of wheat grain images. As primary inputs
takes a TIFF image file of the stained wheat grain, a shapefile (.shp) of the outline of
the wheat grain drawn in ArcMap, an image classification signature file (.gsg) generated
with training samples defined in ArcMap, and the treatment name/code as a text input.
Inputs are also taken for the number of samples within the signature file that represent
area of interest (default = 10), the number of zones to be drawn (default = 5), an op-
tional input file for the widths of the zones to be drawn (.txt) (file containing negative
distance from the aleurone for each zone boundary to be drawn, with each boundary
measurement on a new line), and a scalebar.txt file (two line file: first line is length
of scalebar in um, second line is length of scalebar in pixels) that should be used to
override the default scaling factor for converting from arbitrary units to micrometers,
which will only produce accurate results on the microscope used by the present authors.
Toolbox produces two outputs as csv files: the result of the protein concentration gra-
dient analysis (treatment zones.csv), and of the protein body size-distribution analysis
(treatment spatial.csv). Additionally, an output of the maximum grain width is pro-
duced for each analysis, but not permanently stored, and is overwritten by subsequent
analyses.
A second "RescalingBatch" toolbox is required to rescale input images to a 1x1 cell
size prior to inputting the images into the protein analysis toolbox. This is required to
prevent inaccuracies in the conversion of measurements from pixels to micrometers that
may arise with different microscopes.
