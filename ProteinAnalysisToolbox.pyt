import os, sys, string, arcpy
import numpy as np
import glob
from arcpy import env
from arcpy.sa import *
arcpy.env.pyramid = "NONE"  # Improves processing time
arcpy.env.overwriteOutput=True
arcpy.env.cellSize = 1  # Set the output raster cell size

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
        self.label = "Extracting Proteins"
        self.alias = ""
        # List of tool classes associated with this toolbox
        self.tools = [Protein,RescalingBatch]

# Main class for analysis of microscopy images
# Calculates protein concentration gradient and protein body size-distribution
# data
class Protein(object):
    def __init__(self):
        self.label = 'Extract proteins from the grain'
        self.canRunInBackground = False

    def getParameterInfo(self):
        '''Set up the parameters and return the list of parameter objects.'''

        # Input grain raster image
        param0 = arcpy.Parameter()
        param0.name = 'tifFile'
        param0.displayName = '1. Input grain TIF file:'
        param0.parameterType = 'Required'
        param0.direction = 'Input'
        param0.datatype = 'DEFile'
        param0.filter.list = ['tif', 'tiff', 'TIF']	

        # Input image classification signature file		
        param1 = arcpy.Parameter()
        param1.name = 'signature'
        param1.displayName = '2. Input signature file:'
        param1.parameterType = 'Required'
        param1.direction = 'Input'
        param1.datatype = 'DEFile'
        param1.filter.list = ['gsg']		
		
		# Input number of samples representing protein (default = 10)
        param2 = arcpy.Parameter()
        param2.name = 'number_of_samples'
        param2.displayName = '3. Input number of samples:'
        param2.parameterType = 'Required'
        param2.direction = 'Input'
        param2.datatype = 'GPLong'
        param2.value = 10 # Default value
		
		# Input number of zones  to be drawn (default = 5)
        param3 = arcpy.Parameter()
        param3.name = 'number_of_zones'
        param3.displayName = '4. Input number of zones:'
        param3.parameterType = 'Required'
        param3.direction = 'Input'
        param3.datatype = 'GPLong'
        param3.value = 5 # Default value

		# Input text file with zone distances (optional)
        param4 = arcpy.Parameter()
        param4.name = 'textZones'
        param4.displayName = '5. Input distance zones text file:'
        param4.parameterType = 'Optional'
        param4.direction = 'Input'
        param4.datatype = 'DETextfile'
        param4.filter.list = ['txt']	
		
        # Input outline of grain shapefile
        param5 = arcpy.Parameter()
        param5.name = 'Borders'
        param5.displayName = '6. Select shapefile (*.shp) with grain borders:'
        param5.parameterType = 'Required'
        param5.direction = 'Input'
        param5.datatype = 'Shapefile'
        param5.filter.list = ['shp']
		
		# Input text description field
        param6 = arcpy.Parameter()
        param6.name = 'textDescription'
        param6.displayName = '7. Text description column:'
        param6.parameterType = 'Optional'
        param6.direction = 'Input'
        param6.datatype = 'String'
				
        return [param0,param1,param2,param3,param4,param5,param6]

    def isLicensed(self):
        """Prevent the tool from running if the Spatial Analyst extension is not available."""
        if arcpy.CheckExtension('Spatial') == 'Available':
            return True  # The tool can be executed.
        else:
            return False  # The tool can not be executed.
			
    def updateParameters(self, parameters):
        return
		
    def updateMessages(self, parameters):
        return
		
    def execute(self, parameters, messages):
		# Assigning input parameters to variable names
		grainTIF = parameters[0].valueAsText		
		signatureFILE = parameters[1].valueAsText		
		numberSAMPLES = parameters[2].value
		numbZONES = parameters[3].value
		txtFILEdistances = parameters[4].valueAsText		
		shpBORDER = parameters[5].valueAsText
		description = parameters[6].valueAsText
		
		arcpy.env.overwriteOutput=True
		
		# Calculate Maximum Likelihood (ML) of pixels in input raster image using 
		# image classification signature file
		mlRASTER = MLClassify(grainTIF, signatureFILE)
		
		# Con ML
		# Identifying which pixels were identified as protein by ML 
		# classification (according to numberSAMPLES)
		numberSAMPLESstring = str(numberSAMPLES)
		where_clause = "VALUE >= 1 AND VALUE <= " + numberSAMPLESstring
		protCON = Con(mlRASTER, 1, "", where_clause)
		
		# RasterToPolygon
		folderIN = os.path.dirname(grainTIF)
		polySHP = folderIN + "/" + "shape.shp"
		arcpy.RasterToPolygon_conversion(protCON, polySHP, "NO_SIMPLIFY")
		
		# Constants for conversion to micrometers (um).
		constantLINEAR = 0.32059502436522185175686073352141 # for distance
		constantAREA = 0.10278116964773719279201741905662 # for area
		
		# Calculating scaling for distance (constantLINEAR) and area 
		# (constantAREA)
		# scalebar.txt contains length of scalebar in um on first line, and in 
		# pixels on the second line.
		if os.path.isfile(folderIN + "/" + "scalebar.txt"):
			scale = []
			scaleBAR = open(folderIN + "/" + "scalebar.txt", 'r')
			for linear in scaleBAR:
				scale.append(float(linear))
			scaleBAR.close()
			constantLINEAR = (scale[0])/(scale[1]) # ratio between um and pixels
			constantAREA = constantLINEAR ** 2
		constantLINEARstring = str(constantLINEAR)
		constantAREAstring = str(constantAREA)	
		
		# Add fields and calculate protein area in arbitrary units and 
		# micrometers squared
		arcpy.AddField_management(polySHP, "ProtAreaSc", "DOUBLE")
		arcpy.CalculateField_management(polySHP, "ProtAreaSc",'!shape.area!', "PYTHON_9.3")
		arcpy.AddField_management(polySHP, "ProtAreaMi", "DOUBLE")
		arcpy.CalculateField_management(polySHP, "ProtAreaMi",'!shape.area!*' + constantAREAstring, "PYTHON_9.3")
		arcpy.AddField_management(polySHP, "areaProt", "DOUBLE")
		arcpy.CalculateField_management(polySHP, "areaProt",'!shape.area!*' + constantAREAstring, "PYTHON_9.3")
		arcpy.FeatureToPoint_management(polySHP, folderIN + "/" + "centroids.shp")
		
		# Calculating euclidean distances
		featureTOlineTEMP = folderIN + "/" + "ftl.shp"
		arcpy.FeatureToLine_management(shpBORDER, featureTOlineTEMP)
		outEucDistance = EucDistance(featureTOlineTEMP, "", grainTIF)
		outExtractByMask = ExtractByMask(outEucDistance, shpBORDER)
		outZonalTableTEMP = folderIN + "/" + "ozt.dbf"
		ZonalStatisticsAsTable(shpBORDER, "FID", outExtractByMask, outZonalTableTEMP, "DATA", "ALL")
		distRows = arcpy.da.SearchCursor(outZonalTableTEMP, ['RANGE'])
		distRow = distRows.next()
		maximumDIST = distRow[0]  # Maximum width of grain, aleurone to aleurone
		
		dist = []
		if numbZONES == 0:
			# Reading zone distances text file, if present
			txtFILE = open(txtFILEdistances, 'r')
			for line in txtFILE:
				dist.append(float(line)*(-1))
			txtFILE.close()
		else:
			# Automatically calculating zones based on maximumDIST
			dist.append(-0.001)
			zoneWIDTH = maximumDIST/(float(numbZONES))
			for z in range(1,numbZONES): 
				dist.append((-1)*z*zoneWIDTH)	

		zoneSHP = folderIN + "/" + "zones.shp"
		arcpy.MultipleRingBuffer_analysis(shpBORDER, zoneSHP, dist,"", "", "ALL")
		protein_in_ZONES = folderIN + "/" + "proteinZONES.shp"
		arcpy.Intersect_analysis([zoneSHP,polySHP], protein_in_ZONES) 
		arcpy.CalculateField_management (protein_in_ZONES, "areaProt", '!shape.area!*' + constantAREAstring, "PYTHON_9.3") 
		
		# Add fields and calculate zone distance in arbitrary units and 
		# micrometers
		arcpy.AddField_management(zoneSHP, "distZoneSc", "DOUBLE", 18, 10,"","","NULLABLE")  
		arcpy.CalculateField_management (zoneSHP, "distZoneSc", '!distance! * (-1)', "PYTHON_9.3")  
		arcpy.AddField_management(zoneSHP, "distZoneMi", "DOUBLE", 18, 10,"","","NULLABLE")  
		arcpy.CalculateField_management (zoneSHP, "distZoneMi", '!distance!* (-1)*' + constantLINEARstring, "PYTHON_9.3")
		
		# Add fields and calculate zone area in arbitrary units and micrometers
		arcpy.AddField_management(zoneSHP, "zoneAreaSc", "DOUBLE", 18, 10,"","","NULLABLE")
		arcpy.CalculateField_management (zoneSHP, "zoneAreaSc", '!shape.area!', "PYTHON_9.3") 
		arcpy.AddField_management(zoneSHP, "zoneAreaMi", "DOUBLE", 18, 10,"","","NULLABLE")
		arcpy.CalculateField_management (zoneSHP, "zoneAreaMi", '!shape.area!*' + constantAREAstring, "PYTHON_9.3")
		
		# Add fields and calculate protein area in each zone in arbitrary units 
		# and micrometers
		protein_in_zonesdissolved = folderIN + "/" + "proteinsZONES.shp"	
		arcpy.Dissolve_management(protein_in_ZONES, protein_in_zonesdissolved, "distance")
		arcpy.AddField_management(protein_in_zonesdissolved, "protAreaSc", "DOUBLE", 18, 10,"","","NULLABLE")
		arcpy.CalculateField_management (protein_in_zonesdissolved, "protAreaSc", '!shape.area!', "PYTHON_9.3")
		arcpy.AddField_management(protein_in_zonesdissolved, "protAreaMi", "DOUBLE", 18, 10,"","","NULLABLE")
		arcpy.CalculateField_management (protein_in_zonesdissolved, "protAreaMi", '!shape.area!*' + constantAREAstring, "PYTHON_9.3")
		
		# Join zone background and protein areas and distances. Calculate 
		# percentage protein per zone
		joinedRESULT = folderIN + "/" + "result.shp"
		arcpy.SpatialJoin_analysis(zoneSHP, protein_in_zonesdissolved, joinedRESULT,"","","","CONTAINS")
		arcpy.AddField_management(joinedRESULT, "percent", "DOUBLE", 18, 10,"","","NULLABLE")
		arcpy.CalculateField_management (joinedRESULT, "percent", '(!protAreaMi!/!zoneAreaMi!)*100', "PYTHON_9.3")
		arcpy.DeleteField_management(joinedRESULT,["OID", "Join_Count", "TARGET_FID", "distance", "distance_1"])
		
		# Output results
		arcpy.TableToTable_conversion(joinedRESULT, folderIN, description + "_zones.csv")

		# Clip centroids to outline of grain
		arcpy.Clip_analysis(folderIN + "/" + "centroids.shp", shpBORDER, folderIN + "/" + "centroids2.shp")
		
		# Load distance values for each protein body point
		ExtractValuesToPoints(folderIN + "/" + "centroids2.shp", outEucDistance, folderIN + "/" + "proteinCentroids.shp")
		
		# Add fields and calculate protein body distances and sizes, and 
		# description of treatment
		arcpy.AddField_management(folderIN + "/" + "proteinCentroids.shp", "dist_Sc", "FLOAT",18, 10,"","","NULLABLE")
		arcpy.AddField_management(folderIN + "/" + "proteinCentroids.shp", "dist_Mi", "FLOAT",18, 10,"","","NULLABLE")
		arcpy.AddField_management(folderIN + "/" + "proteinCentroids.shp", "treatment", "TEXT","","","","","NULLABLE")
		arcpy.CalculateField_management(folderIN + "/" + "proteinCentroids.shp", "dist_Sc","!RASTERVALU!","PYTHON_9.3")
		arcpy.CalculateField_management(folderIN + "/" + "proteinCentroids.shp", "dist_Mi", "!RASTERVALU!*" + constantLINEARstring, "PYTHON_9.3") 
		arcpy.CalculateField_management(folderIN + "/" + "proteinCentroids.shp", "treatment", "'" + description + "'" , "PYTHON_9.3")
		arcpy.DeleteField_management(folderIN + "/" + "proteinCentroids.shp", ["ID", "GRIDCODE", "ORIG_FID", "RASTERVALU", "areaProt"])	
		
		# Output results
		arcpy.TableToTable_conversion(folderIN + "/" + "proteinCentroids.shp", folderIN, description + "_spatial.csv")

		arcpy.SpatialJoin_analysis(folderIN + "/" + "shape.shp", folderIN + "/" + "proteinCentroids.shp", folderIN + "/" + "shapeTemp.shp")
		arcpy.Clip_analysis(folderIN + "/" + "shapeTemp.shp", shpBORDER, folderIN + "/" + "proteins.shp")
		arcpy.DeleteField_management(folderIN + "/" + "proteins.shp",["TARGET_FID", "ID", "GRIDCODE", "Join_Count", "areaProt_1"])
		
		# Delete temporary files used
		arcpy.Delete_management(folderIN + "/" + "centroids.shp")
		arcpy.Delete_management(folderIN + "/" + "centroids2.shp")
		arcpy.Delete_management(folderIN + "/" + "ftl.shp")
		arcpy.Delete_management(folderIN + "/" + "ozt.dbf")
		arcpy.Delete_management(protein_in_ZONES)
		arcpy.Delete_management(folderIN + "/" + "shapeTemp.shp")
		arcpy.Delete_management(folderIN + "/" + "shape.shp")
		arcpy.Delete_management(joinedRESULT)
		
		arcpy.RefreshCatalog(folderIN)
				
		return

# RescalingBatch class used to rescale a batch of microscopy images to 1x1
# cell sizes
# This is required for input images used in the Protein class to ensure 
# correct scaling from pixels to micrometers
class RescalingBatch(object):
    def __init__(self):
        self.label = 'Batch rescaling image to 1 x 1 pixel size'
        self.canRunInBackground = False

    def getParameterInfo(self):
        '''Set up the parameters and return the list of parameter objects.'''

        # Input grain raster
        param0 = arcpy.Parameter()
        param0.name = 'tifFile'
        param0.displayName = 'Input grain TIF file:'
        param0.parameterType = 'Required'
        param0.direction = 'Input'
        param0.datatype = 'DEFolder'
        param0.filter.list = ['tif', 'tiff', 'TIF']	
		
        return [param0]

    def isLicensed(self):
        """Prevent the tool from running if the Spatial Analyst extension is not available."""
        if arcpy.CheckExtension('Spatial') == 'Available':
            return True  # The tool can be executed.
        else:
            return False  # The tool can not be executed.
			
    def updateParameters(self, parameters):
        return
		
    def updateMessages(self, parameters):
        return
		
    def execute(self, parameters, messages):
		# Calculations
		folder = parameters[0].valueAsText
		arcpy.env.cellSize = 1
		
		for input_file in glob.glob(os.path.join(folder,'*.tif')):
			filename = os.path.basename(input_file)
			input = input_file
			a = Raster(input)
			rescaleRATIO = 1/a.meanCellHeight
			b = folder + "/" + "pxl_1x1_" + filename
			arcpy.Rescale_management(a,b,str(rescaleRATIO),str(rescaleRATIO))
		
		return