import os
import sys
#import pySW4 as sw4
import matplotlib.pyplot as plt
import numpy as np
import json

#if this is false negetive numbers count as below sea level (the oppsite of how they are in SW4 rFiles)
seaLevel = False

#decide if I want to inherit from block or not, may not actually be that useful
#class Block(sw4.prep.rfileIO.Block):
class Block():
    def __init__(self,blockControl):
        self.number = blockControl["BLOCK_NUMBER"]
        self.hh = blockControl["HHb"]
        self.hv = blockControl["HVb"]
        self.z0 = blockControl["Z0"]
        self.nc = blockControl["NCb"]
        self.ni = blockControl["Ni"]
        self.nj = blockControl["Nj"]
        self.nk = blockControl["Nk"]
        
        #-999 is the null data type which is ignored durring interpolation etc.
        if(self.nc != 1):
            self.vp = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
            self.vs = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
            self.p = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
        if(self.nc > 3):
            self.qp = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
            self.qs = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)

        elif(self.nc ==1):
            #topo block
            self.topo = np.full((self.ni,self.nj),-999,dtype=np.float32)
            
        self.x_extent = (0, (self.ni - 1) * self.hh * 1e-3)
        self.y_extent = (0, (self.nj - 1) * self.hh * 1e-3)
        self.z_extent = ((self.z0 + (self.nk - 1) * self.hv) * 1e-3,
                 self.z0 * 1e-3)

        self.xyextent = self.y_extent + self.x_extent
        self.xzextent = self.x_extent + self.z_extent
        self.yzextent = self.y_extent + self.z_extent
        
    def linearlyInterpolate(self):
        #interpolates to fill all of the block data which is missing (i.e for which there is no data)
        pass
        
    def project(self):
        #for use on topography data block 1 only (basically cuts off everything below the surface)
        pass

class Parameterfile():
    def __init__(self,pFname = "rFile.ini"):
        self.pfContents = []            
        #open the Parameterfile and load its contents
        self.loadPF(pFname)        

    def loadPF(self,pFname):
        with open(pFname,'r') as fileObject:
            self.pfContents = [line.strip().encode('string-escape') for line in fileObject if("#" not in line.strip())]
        #attempt to jsonify
        #r is to indicate that this is a raw string (i.e ignore \ escape charecter
        self.pfContents = json.loads(r''.join(self.pfContents))


class Model():
    def __init__(self,pfName = "rFile.ini",fname="untrackedModels/GFM_all_clean"):   
        self.Parameterfile = Parameterfile(pfName)                
        #build the blocks as specified
        self.blocks = [Block(self.Parameterfile.pfContents["BLOCK_CONTROL"][i]) for i in self.Parameterfile.pfContents["BLOCK_CONTROL"].keys() if (i!="HEADER")]
        #load the information from the model descritption
        self.loadModel(fname)
        #find the datums
        self.minX = min(self.ModelFileData[:,0])        
        self.minY = min(self.ModelFileData[:,1])
               
    #clean this up to make it faster at a later point
    def loadModel(self,fname):
        lines = []
        self.header = []
        self.ModelFileData = []
        with open( fname,'r') as fileObject:
            while True:
                #get the header
                line =  fileObject.readline()
                #save the line
                self.header.append(line)
                if("#" not in line):
                    break                    
                    
            for line in fileObject:           
                #assign everything to the nearest point
                #the multiplication by 1000 moves from m/s to km/s
                lines = [i+" " for i in line.split() if i != ""]                
                self.ModelFileData.append((float(lines[0]),float(lines[1]),float(lines[2]),int(lines[3])))
                
        #now convert it to a numpy array for easier use
        self.ModelFileData = np.array(self.ModelFileData,dtype=np.float32)
        #correct the up down orientation
        self.ModelFileData = np.flipud(self.ModelFileData)
        #flip the negetives so that it is in the same standard as the Rfile
        ### VERY IMPORTANT!!! ###
        self.ModelFileData[:,2] = self.ModelFileData[:,2] * -1
        
        
    #fix coordinates (i.e use cartesian coordinates)
    def fixCoordinates(self):
        #find the minimum value for all of the UTM coordinates and assign accordingly
        
        pass
        
    """
    coordsToIndex
    
    Parameters
    ----------
    east: float
        UTM easting coordinate for the bottom left most coordinate in rFile

    north : float
        UTM northing coordinate for the bottom left most coordinate in rFile     

    depth : float
        current depth in model (note that this follows the rFile depth convention, negetiv == above sea level depth
        

    depthDelta : float
        Depth delta (i.e 10 means that every pixel is sepperated by 10M depth)

    horrizontalDelta: float
        horrizontal delta (i.e 10 means that every pixel in xy space is 10M square)

    Ni : int
        #i coordinates for this block
        
    Ni : bool
        #j coordinates for this block
        
    Nk: 
        #k coordinates for this block
        
    """
        
    def coordsToIndex(self,depthDatum,depthDelta,horrizontalDelta,Ni,Nj,Nk,easting,northing,depthVal):        
        #datums which are always the same since each block will always have the same datum  (i.e rfiles describe 
        eastingDatum = self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["LAT0"]
        northingDatum = self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["LON0"]
        #return the index (this should map everything to the nearest pixel? test IT EMPIRICALLY!!
        return int(((easting-eastingDatum)/(Ni*horrizontalDelta)) + ((northing-northingDatum)/(Nj*horrizontalDelta))
        + ((depthVal-depthDatum)/(Nk*depthDelta)))
    
        
    def buildRfile(self,fileObject): 
        blockExtent = []     
        blockIndex = 1  
        rasterDex = 0
        easting = 0
        northing = 0
        #I am just going to do these all in one go to start with and I will deal with memory issues later        
        #construct topography block, assign to 0 for now, then once I have the data assign to the data
        #TODO assign this to topography correctly
        self.blocks[0].topo[:] = 0
        
        #save the rFile header
        
        #save the block headers
        
        #save datums from header and remove all other uneeded crap
        easting = self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]['LAT0'] 
        northing = self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]['LON0']
        
        #save the block data
        for block in self.Parameterfile.pfContents["BLOCK_CONTROL"]:
            #make sure that I am dealing with an actualy BLOCK
            if(block != "TOPO" and block != "HEADER"):
                #subset the input data by the current block (basically get a list of the indexes that correspond to data that is within the current block
                #takes everything within (inclusive) depth range for this block)
                blockExtent = np.where(np.logical_and(np.less_equal(self.ModelFileData[:,2],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"]),np.greater_equal(self.ModelFileData[:,2],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"] 
                - (self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"]*self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"]))))[0]
                
                for i in blockExtent:
                    #compute the correct index for each point in the raster!                
                    #depthDatum,depthDelta,horrizontalDelta,Ni,Nj,Nk,easting,northing,depthVal)
                    #TODO fix this function since it appears to be just returning zero!
                    rasterDex = self.coordsToIndex(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"]
                    ,self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HHb"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Ni"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nj"],
                    self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"],easting,northing,self.ModelFileData[i][2])
                    
                    print(rasterDex)
                    #assign each point to the model raster(s) as specified
                    if(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["NCb"] != 1):
                        self.blocks[blockIndex].vp[rasterDex] = 0
                        self.blocks[blockIndex].vs[rasterDex] = 0
                        self.blocks[blockIndex].p[rasterDex] = 0
                    
                    if(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["NCb"]  > 3):
                        self.blocks[blockIndex].qp[rasterDex] = 0
                        self.blocks[blockIndex].qs[rasterDex] = 0
                
                #linearly interpolate this block
                
                
                #save this block
                
                
                #(optional) return the memory used by this block
                
                #next block 
                blockIndex += 1
                pass
            
        """
        #build the blocks for everything else        
        for block in range(1,len(self.blocks)):
            #build each block and populate with model data
            for j in range()
            pass
        """
        
        #halt
        #free memory held by base model 
        #now linearly interpolate to fill all of the blocks
        #write the rFile
        pass
    
    #get all of the charecteristics
    #TODO add support to compute this from model specified in parameter file, that is WHY THIS IS A FUNCTION INSTEAD OF JUST AND ASSIGNMENT
    #I PLAN TO BUILD ONTO THESE IN ORDER TO MAKE MORE SOPHISTICATED GOELOGIC MODELS!
    def getVP(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit)]["VP"]
    
    def getVS(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit)]["VS"]

    def getQP(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit)]["QP"]
    
    def getQS(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit)]["QS"]
    
    def getP(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit)]["P"]


def main():
    #allocate pixel space
    
    #load the geologic model and map each point to pixel space (at the same time)
    
    #for each block
    
        #interpolate
        
        #save pixel space to the rFile
        
        #reallocate memory 
        
    #done
    
    pass

if __name__=="__main__":
    main()
    pass


