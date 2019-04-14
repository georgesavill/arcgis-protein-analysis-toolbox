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


Citation:
@article{10.1093/jxb/ery127,
    author = {Savill, George P and Buchner, Peter and Wan, Yongfang and Hawkesford, Malcolm J and Michalski, Adam and Powers, Stephen J and Tosi, Paola},
    title = "{Temperature and nitrogen supply interact to determine protein distribution gradients in the wheat grain endosperm}",
    journal = {Journal of Experimental Botany},
    volume = {69},
    number = {12},
    pages = {3117-3126},
    year = {2018},
    month = {04},
    abstract = "{Gradients exist in the distribution of storage proteins in the wheat (Triticum aestivum) endosperm and determine the milling properties and protein recovery rate of the grain. A novel image analysis technique was developed to quantify both the gradients in protein concentration, and the size distribution of protein bodies within the endosperm of wheat plants grown under two different (20 or 28 °C) post-anthesis temperatures, and supplied with a nutrient solution with either high or low nitrogen content. Under all treatment combinations, protein concentration was greater in the endosperm cells closest to the aleurone layer and decreased towards the centre of the two lobes of the grain, i.e. a negative gradient. This was accompanied by a decrease in size of protein bodies from the outer to the inner endosperm layers in all but one of the treatments. Elevated post-anthesis temperature had the effect of increasing the magnitude of the negative gradients in both protein concentration and protein body size, whilst limiting nitrogen supply decreased the gradients.}",
    issn = {0022-0957},
    doi = {10.1093/jxb/ery127},
    url = {https://doi.org/10.1093/jxb/ery127},
    eprint = {http://oup.prod.sis.lan/jxb/article-pdf/69/12/3117/25088854/ery127.pdf},
}
