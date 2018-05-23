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
        
        #Do not predclare any of this use a mesh instead!
        """
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
        """
        
        self.x_extent = (0, (self.ni - 1) * self.hh * 1e-3)
        self.y_extent = (0, (self.nj - 1) * self.hh * 1e-3)
        self.z_extent = ((self.z0 + (self.nk - 1) * self.hv) * 1e-3,
                 self.z0 * 1e-3)

        self.xyextent = self.y_extent + self.x_extent
        self.xzextent = self.x_extent + self.z_extent
        self.yzextent = self.y_extent + self.z_extent
        

        
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
                #fix everything according to the datums           
                #also internally I assign 0 to be the bottome of the rfile!b             
                x = int(float(lines[0])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["xDatum"])  
                y = int(float(lines[1])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["yDatum"])
                z = int(float(lines[2]) + self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"] if float(lines[2]) > 0 else abs(float(lines[2]) + self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"]))
                
                self.ModelFileData.append((x,y,z,int(lines[3])))
                
                
                #self.ModelFileData.append((int(float(lines[0])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["xDatum"]),int(float(lines[1])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["yDatum"]),float(lines[2]),int(lines[3]) if int(lines[3])))
                
        #now convert it to a numpy array for easier use
        self.ModelFileData = np.array(self.ModelFileData,dtype=np.float32)
        #correct the up down orientation
        self.ModelFileData = np.flipud(self.ModelFileData)

        
        
        
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
        
    def coordsToIndex(self,depthDatum,depthDelta,horrizontalDelta,Ni,Nj,Nk,xDatum,yDatum,depthVal,xVal,yVal):        
        #datums which are always the same since each block will always have the same datum  (i.e rfiles describe 
        #return the index (this should map everything to the nearest pixel? test IT EMPIRICALLY!!
        return int(((xVal-xDatum)/(Ni*horrizontalDelta)) + ((yVal-yDatum)/(Nj*horrizontalDelta))
        + ((depthVal-abs(depthDatum))/(Nk*depthDelta)))
            
    def buildRfile(self,fileObject): 
        blockExtent = []     
        blockIndex = 1  
        dataArrayDimensions = []
        xCoord = 0
        yCoord = 0
        zCoord = 0
        #I am just going to do these all in one go to start with and I will deal with memory issues later        
        #construct topography block, assign to 0 for now, then once I have the data assign to the data
        #TODO assign this to topography correctly
        self.blocks[0].topo = np.full((self.blocks[0].ni,self.blocks[0].nj),0,dtype=np.float32)
        self.blocks[0].topo[:] = 0        

        #save datums from header and remove all other uneeded crap
   
        #save the block data
        for block in self.Parameterfile.pfContents["BLOCK_CONTROL"]:
            #make sure that I am dealing with an actualy BLOCK
            if(block != "TOPO" and block != "HEADER"):                               
                #subset the input data by the current block (basically get a list of the indexes that correspond to data that is within the current block
                #takes everything within (inclusive) depth range for this block)
                blockExtent = np.where(np.logical_and(np.less_equal(self.ModelFileData[:,2],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"]),np.greater_equal(self.ModelFileData[:,2],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"] 
                - (self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"]*self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"]))))[0]
                
                #compute the data array dimensions (note these may exceed i,j and k)
                dataArrayDimensions = [self.blocks[blockIndex].ni,self.blocks[blockIndex].nj,self.blocks[blockIndex].nk]
                self.blocks[blockIndex].vp = np.nan*np.empty(dataArrayDimensions)
                #now assign to this based on the availibility of points in the underlying model
                for i in range(len(self.ModelFileData[blockExtent][:,0])):    
                    #take the coordinate, scale according to ni and hh and assign
                    xCoord = int(self.ModelFileData[blockExtent][i][0]/(self.blocks[blockIndex].ni*self.blocks[blockIndex].hh))
                    yCoord = int(self.ModelFileData[blockExtent][i][1]/(self.blocks[blockIndex].nj*self.blocks[blockIndex].hh))
                    #this one is a little tricky because I have to fix the datum
                    zCoord = -1
                    #do other things                    
                    #self.blocks[blockIndex].vp[][][]                                        
                    #do this for all components!
                    
                    #done
                                                                                       
                #assign based on the number of components availible                                         
                """
                as shown by https://stackoverflow.com/questions/30764955/python-numpy-create-2d-array-of-values-based-on-coordinates
                x = [0, 0, 1, 1, 2, 2]
                y = [1, 2, 0, 1, 1, 2]
                z = [14, 17, 15, 16, 18, 13]
                z_array = np.nan * np.empty((3,3))
                z_array[y, x] = z
                """           
                
                """
                #actually, downsize the x,y,z indexes based on the datums and deltas in order to make all of this fit correctly
                TODO figure out why this interpolation is failing!!
                ACTUALLY DEFINE DIMENSIONS BASED ON THE WHOLE VOLUME, THEN ASSIGN VALUES IN THAT VOLUME!!! (so each pixel is a meter!)
                THEN IF RESOLUTION IS TO HIGH BLOCK!!
                I think that the problem is that you need to define the dimensions of the mesh correctly, and then compute deltas accordingly!
                I am like 90% sure that will fix this!!
                """                
                #correct deltas to maintain model area (i.e no stretch)
                
                #assign each data point to its respective mesh (if possible) (note assume base resolution of at least 1 meter               
                """
                self.blocks[blockIndex].vp = np.nan*np.empty([self.blocks[blockIndex].ni*self.blocks[blockIndex].hh/( self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"]
                ),self.blocks[blockIndex].nj*self.blocks[blockIndex].hh/( self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"]
                ),self.blocks[blockIndex].nk*self.blocks[blockIndex].hv/( self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"]
                )])
                """
                #divide by deltas and subtract datum to fix these indexes
                #self.blocks[blockIndex].vp[((self.ModelFileData[blockExtent][:,0])/self.blocks[blockIndex].hh).astype(int),((self.ModelFileData[blockExtent][:,1])/self.blocks[blockIndex].hh).astype(int),(self.ModelFileData[blockExtent][:,2]/self.blocks[blockIndex].hv).astype(int)] = [self.getVP(unit) for unit in self.ModelFileData[blockExtent][:,3]]
                                
                """
                if(self.blocks[block].nc != 1):
                    self.vp = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
                    self.blocks[block].vp = np.nan*np.empty((self.blocks[block].ni,self.blocks[block].nj,self.blocks[block].nk))
                    
                    self.vs = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
                    self.p = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
                
                if(self.blocks[block].nc > 3):
                    self.qp = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)
                    self.qs = np.full((self.ni,self.nj,self.nk),-999,dtype=np.float32)


                #assign each data dimension (i.e vp, vs p etc...)
                self.vp = np.nan * np.empty(())
                """
                                
                #linearly interpolate this block
                                                
                #save this block
                                
                #(optional) return the memory used by this block
                
                #next block 
                blockIndex += 1
                pass
            

        #halt
        #free memory held by base model 
        #now linearly interpolate to fill all of the blocks
        #write the rFile
        pass
    
    #get all of the charecteristics
    #TODO add support to compute this from model specified in parameter file, that is WHY THIS IS A FUNCTION INSTEAD OF JUST AND ASSIGNMENT
    #I PLAN TO BUILD ONTO THESE IN ORDER TO MAKE MORE SOPHISTICATED GOELOGIC MODELS!
    def getVP(self,unit):
        #figure out a better way to deal with this SLIT function!
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit).split(".")[0]]["VP"]
    
    def getVS(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit).split(".")[0]]["VS"]

    def getQP(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit).split(".")[0]]["QP"]
    
    def getQS(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit).split(".")[0]]["QS"]
    
    def getP(self,unit):
        return self.Parameterfile.pfContents["GEOL_CONTROL"][str(unit).split(".")[0]]["P"]


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


