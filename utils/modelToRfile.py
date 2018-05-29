import os
import sys
#import pySW4 as sw4
import matplotlib.pyplot as plt
import numpy as np
import json
#from scipy.spatial import KDTree
#from scipy.ndimage import *
#from scipy.interpolate.griddata import *
import scipy
import scipy.ndimage
from scipy.interpolate import griddata

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
        
        #TODO remove these, I do not think that I need them at all!
        self.x_extent = (0, (self.ni - 1) * self.hh * 1e-3)
        self.y_extent = (0, (self.nj - 1) * self.hh * 1e-3)
        self.z_extent = ((self.z0 + (self.nk - 1) * self.hv) * 1e-3,
                 self.z0 * 1e-3)

        self.xyextent = self.y_extent + self.x_extent
        self.xzextent = self.x_extent + self.z_extent
        self.yzextent = self.y_extent + self.z_extent
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
               
    #clean this up to make it faster at a later point
    def loadModel(self,fname):
        lines = []
        self.header = []
        self.ModelFileData = []
        numLines = 0
        dataStart = 0
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
                x = int((float(lines[0])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["xDatum"])/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])
                y = int((float(lines[1])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["yDatum"])/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])
                z = int(abs(float(lines[2]) - self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"]) if float(lines[2]) > 0 else abs(float(lines[2])) + self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"])/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"]
                
                self.ModelFileData.append((x,y,z,int(lines[3])))                                
                #self.ModelFileData.append((int(float(lines[0])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["xDatum"]),int(float(lines[1])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["yDatum"]),float(lines[2]),int(lines[3]) if int(lines[3])))
                
        #now convert it to a numpy array for easier use
        self.ModelFileData = np.array(self.ModelFileData,dtype=np.float32)

        #correct the up down orientation
        #self.ModelFileData = np.flipud(self.ModelFileData)

        #find datums
        self.minX = min(self.ModelFileData[:,0])        
        self.minY = min(self.ModelFileData[:,1])
        #put this into 3d grid space        
        x,y,z=np.mgrid[0:len(np.unique(self.ModelFileData[:,0])),0:len(np.unique(self.ModelFileData[:,1])) ,0:len(np.unique(self.ModelFileData[:,2]))]

        #assign data to "nearest" point on the mesh
        self.inputModel = griddata((self.ModelFileData[:,0].astype(int),self.ModelFileData[:,1].astype(int),self.ModelFileData[:,2]),(self.ModelFileData[:,3]).astype(int),(x,y,z),method='nearest')
        
        #free up memory by removing the origional data
        #TODO uncomment this for the actualy production code!
        #del self.ModelFileData                                
        
    #fix coordinates (i.e use cartesian coordinates)
    def fixCoordinates(self):
        #find the minimum value for all of the UTM coordinates and assign accordingly
        pass
        
    def buildRfile(self,fileObject):   
        i =0
        j=0
        k=0
        x=0
        y=0
        z=0
        vp=[]
        vs=[]
        qp=[]
        qs=[]
        p=[]
        componentMesh = []
        blockIndex = 1  
        baseDepth = 0
        topBlock = 0

        #I am just going to do these all in one go to start with and I will deal with memory issues later        
        #construct topography block, assign to 0 for now, then once I have the data assign to the data
        #TODO assign this to topography correctly
        self.blocks[0].topo = np.full((self.blocks[0].ni,self.blocks[0].nj),0,dtype=np.float32)
        self.blocks[0].topo[:] = 0        
        #save datums from header and remove all other uneeded crap
        
        #save the block header
        
        #now save the block data
        for block in self.Parameterfile.pfContents["BLOCK_CONTROL"]:
            #make sure that I am dealing with an actualy BLOCK
            if(block != "TOPO" and block != "HEADER"):    
                print("Iteriation ",block)                           
                #construct component mesh for points in this space
                meshGrid = np.mgrid[0:self.blocks[blockIndex].ni,0:self.blocks[blockIndex].nj,0:self.blocks[blockIndex].nk]                
                
                #compute the depth datum--using my modified coordinates
                """
                topBlock = (((abs(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"]) - self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"])*self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"]/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])) if (self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"]) > 0 else (abs(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"] 
                *self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"]/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])
                baseDepth = ((abs(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"]) - self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"])/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])) if (self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"]) > 0 else (abs(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"] 
                /self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])
                """
                topBlock = self.computeTop(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"],
                self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"])
                
                baseDepth = self.computeBottom(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"])

                #compute coordinate sets
                x,y,z = np.meshgrid(np.arange(self.inputModel.shape[0]),np.arange(self.inputModel.shape[1]),np.arange(self.inputModel[:,:,baseDepth:topBlock].shape[2]))                           
                #x,y,z = self.inputModel[:,:,topBlock:baseDepth].shape
                #now assign to this based on the availibility of points in the underlying model                                
                #test = scipy.interpolate.griddata((x,y,z,([self.getVP(i) for i in self.inputModel[:,:,topBlock:baseDepth].flatten()]),(meshGrid[0],meshGrid[1],meshGrid[2]),method='nearest'))
                X = meshGrid[0]
                Y = meshGrid[1]
                Z = meshGrid[2]
                #TODO take expand on this so that these functions can read and understand youre itnerpolation scheme
                #vp
                print(baseDepth,topBlock,)
                blockData = np.asarray([self.getVP(i) for i in np.nditer(self.inputModel[:,:,baseDepth:topBlock])])
                #the flattening is because what I essentially have in grid data is a basis and I want EVERY coordinate pair!
                vp = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(meshGrid[0],meshGrid[1],meshGrid[2]),method='nearest')
                #vs
                #blockData = np.asarray([self.getVS(i) for i in np.nditer(self.inputModel[:,:,topBlock:baseDepth])])
                #vs = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(meshGrid[0],meshGrid[1],meshGrid[2]),method='linear')
                
                #p
                #blockData = np.asarray([self.getP(i) for i in np.nditer(self.inputModel[:,:,topBlock:baseDepth])])
                #p = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(meshGrid[0],meshGrid[1],meshGrid[2]),method='linear')
                                
                #qp
                #blockData = np.asarray([self.getQP(i) for i in np.nditer(self.inputModel[:,:,topBlock:baseDepth])])
                #qp = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(meshGrid[0],meshGrid[1],meshGrid[2]),method='linear')
                                
                #qs
                #blockData = np.asarray([self.getQS(i) for i in np.nditer(self.inputModel[:,:,topBlock:baseDepth])])
                #qs = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(meshGrid[0],meshGrid[1],meshGrid[2]),method='linear')
                #now save it all to block data
                
                
                #next block 
                blockIndex += 1

                
        print("finished!")
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
    
    #computes the bottom of the zRange for the current block
    @staticmethod
    def computeBottom(z0,baseResolution,datum):
        if(z0<0):
            return (abs(z0)+datum)/baseResolution        
        else:
            return (z0-datum)/baseResolution
            pass
            
    @staticmethod
    def computeTop(z0,baseResolution,datum,HVb,Nk):
        if(z0<0):
            return (abs(z0)+datum+Nk*HVb)/baseResolution
        else:
            return (z0-datum+Nk*HVb)/baseResolution
            pass





def main():
    # find 10 nearest points
    tree = KDTree(a, leafsize=a.shape[0]+1)
    distances, ndx = tree.query([point], k=10)
    # print 10 nearest points to the chosen one
    print a[ndx]
    pass
