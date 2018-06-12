import os
import sys
#import pySW4 as sw4
import matplotlib.pyplot as plt
import numpy as np
import json
import scipy
import scipy.ndimage
from scipy.interpolate import griddata
#pySW4 dependency
import pySW4Rfile


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
    def __init__(self,pfName = "rFile.ini",fname="untrackedModels/GFM_all_clean",topo="/home/eeckert/git/pySW4Forked/pySW4/utils/untrackedModels/NNSS_VOLUME/Topography.2grd.dat"):   
        self.Parameterfile = Parameterfile(pfName)                
        #build the blocks as specified
        self.blocks = [Block(self.Parameterfile.pfContents["BLOCK_CONTROL"][i]) for i in self.Parameterfile.pfContents["BLOCK_CONTROL"].keys() if (i!="HEADER")]
        #load the information from the model descritption
        self.loadModel(fname)
        self.loadTopo(topo)
        
    def loadTopo(self,topo):
        self.topoData = []
        lines = []
        header = []
        numLines = 0
        dataStart = 0
        if(topo == "NoTopo" or topo == ""):
            print("No topography data given, setting topo to flat surface == sea level")
            self.topoData = ""
            pass
            
        else:
            #load up the topography information
            with open( topo,'r') as fileObject:
                while True:
                    #skip over the header which I do not care about at all
                    line =  fileObject.readline()
                    if("#" not in line):
                        break                    
                       
                for line in fileObject:           
                    #assign everything to the nearest point
                    #the multiplication by 1000 moves from m/s to km/s
                    lines = [i+" " for i in line.split() if i != ""]  
                    #fix everything according to the datums     
                    #TODO make sure that topography always has the same extent!      
                    x = int((float(lines[0])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["xDatum"])/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])
                    y = int((float(lines[1])-self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["yDatum"])/self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"])
                    #also internally I assign 0 to be the bottome of the rfile!b             
                    self.topoData.append((x,y,float(lines[2])))
            #drop everything that is out of extent
            self.topoData = [i for i in self.topoData if (i[0] >= 0 and i[1] >= 0)]
            #fix to np array
            self.topoData = np.array(self.topoData)
            #now stick it to the grid size specified by the user (generating coordinates)
            x,y=np.mgrid[0:len(np.unique(self.ModelFileData[:,0])),0:len(np.unique(self.ModelFileData[:,1]))]
            #assign, I am using linear interpolation because "Averaging" makes more sense to me than just assigning elevations to whatever happend to be closest
            self.topoData = griddata((self.topoData[:,0],self.topoData[:,1]),(self.topoData[:,2]),(x,y),method='linear')
            #free up memory by removing the origional data
            #TODO uncomment this for the actualy production code!
            #del self.ModelFileData
    

            
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
        
        #put this into 3d grid space        
        x,y,z=np.mgrid[0:len(np.unique(self.ModelFileData[:,0])),0:len(np.unique(self.ModelFileData[:,1])) ,0:len(np.unique(self.ModelFileData[:,2]))]

        #assign data to "nearest" point on the mesh
        self.inputModel = griddata((self.ModelFileData[:,0].astype(int),self.ModelFileData[:,1].astype(int),self.ModelFileData[:,2]),(self.ModelFileData[:,3]).astype(int),(x,y,z),method='nearest')
        
                             
        
    #fix coordinates (i.e use cartesian coordinates)
    def fixCoordinates(self):
        #find the minimum value for all of the UTM coordinates and assign accordingly
        pass
        
    def buildRfile(self,fileObject):   

        x,X=0,0
        y,Y=0,0
        z,Z=0,0
        vp=[]
        vs=[]
        qp=[]
        qs=[]
        p=[]
        topo = []
        componentMesh = []
        blockIndex = 1  
        baseDepth = 0
        topBlock = 0
        
        #save the rFile header
        pySW4Rfile.write_hdr(fileObject,1,4,1,self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["AZIMUTH"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["LAT0"],
        self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["LON0"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["PROJECTION_STRING"],len(self.Parameterfile.pfContents["BLOCK_CONTROL"])-1)

        #I am just going to do these all in one go to start with and I will deal with memory issues later        
        #construct topography block, assign to 0 for now, then once I have the data assign to the data
        #TODO assign this to topography correctly
        
        #assign topography header                                                                   1=HHv ,z0 unuesd, 1 component 
        pySW4Rfile.write_block_hdr(fileObject, self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["HHb"], 1, 0, 1, self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["Ni"],
        self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["Nj"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["Nk"])
        
        #assign headers for all other blocks
        for block in self.Parameterfile.pfContents["BLOCK_CONTROL"]:
            #make sure that I am dealing with an actualy BLOCK
            if(block != "TOPO" and block != "HEADER"):
                #write_block_hdr(f, hh, hv, z0, nc, ni, nj, nk):
                pySW4Rfile.write_block_hdr(fileObject, self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HHb"], 
                self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"], self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"],
                self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["NCb"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Ni"],
                self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nj"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"])
        
        #now compute the actual data
        if(type(self.topoData)!=np.ndarray):
            #assign the block contents to the lowest elevation in the input data
            topo = np.full((self.blocks[0].ni,self.blocks[0].nj),0,dtype=np.float32)
            topo[:] = 0
            #save it
            pySW4Rfile.write_topo_block(fileObject, topo)

        else:
            #use the topo loaded form the file
            X,Y = np.meshgrid(np.arange(0, self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["Ni"]),np.arange(0, self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["Nj"])) 
            x,y = np.meshgrid(np.arange(self.inputModel.shape[0]),np.arange(self.inputModel.shape[1]))                           

            vp = griddata((x.ravel(),y.ravel()),self.topoData.ravel(),(X,Y),method='nearest')

            pySW4Rfile.write_topo_block(fileObject,vp)
            print("wrote topogaphy!")
 
        #now save the block data
        for block in self.Parameterfile.pfContents["BLOCK_CONTROL"]:
            #make sure that I am dealing with an actualy BLOCK
            if(block != "TOPO" and block != "HEADER"):    
                print("DEBUG: Iteriation ",block)                           
                #construct component mesh for points in this space
                meshGrid = np.mgrid[0:self.blocks[blockIndex].ni,0:self.blocks[blockIndex].nj,0:self.blocks[blockIndex].nk]                
                
                #compute the depth datum--using my modified coordinates
       
                topBlock = self.computeTop(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"])

                baseDepth = self.computeBottom(topBlock,self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseResolution"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"])

                #compute coordinate sets
                x,y,z = np.meshgrid(np.arange(self.inputModel.shape[0]),np.arange(self.inputModel.shape[1]),np.arange(self.inputModel[:,:,baseDepth:topBlock].shape[2]))                           
                #coordinate that I want
                X,Y,Z = np.meshgrid(np.arange(0, self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Ni"]),np.arange(0, self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nj"]),np.arange(0, self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"])) 
                #now assign to this based on the availibility of points in the underlying model                                
                #test = scipy.interpolate.griddata((x,y,z,([self.getVP(i) for i in self.inputModel[:,:,topBlock:baseDepth].flatten()]),(meshGrid[0],meshGrid[1],meshGrid[2]),method='nearest'))

                #TODO take expand on this so that these functions can read and understand youre itnerpolation scheme
                #vp
                print(baseDepth,topBlock)
                blockData = np.asarray([self.getVP(i) for i in np.nditer(self.inputModel[:,:,baseDepth:topBlock])])
                #the flattening is because what I essentially have in grid data is a basis and I want EVERY coordinate pair!
                
                vp = griddata((x.ravel(),y.ravel(),z.ravel()),blockData,(X,Y,Z),method='nearest')
                #vs
                blockData = np.asarray([self.getVS(i) for i in np.nditer(self.inputModel[:,:,baseDepth:topBlock])])
                vs = griddata((x.ravel(),y.ravel(),z.ravel()),blockData,(X,Y,Z),method='nearest')
                
                #p
                blockData = np.asarray([self.getP(i) for i in np.nditer(self.inputModel[:,:,baseDepth:topBlock])])
                p = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(X,Y,Z),method='nearest')
                                
                #qp
                blockData = np.asarray([self.getQP(i) for i in np.nditer(self.inputModel[:,:,baseDepth:topBlock])])
                qp = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(X,Y,Z),method='nearest')
                                
                #qs
                blockData = np.asarray([self.getQS(i) for i in np.nditer(self.inputModel[:,:,baseDepth:topBlock])])
                qs = griddata((x.flatten(),y.flatten(),z.flatten()),blockData,(X,Y,Z),method='nearest')
                #now save it all to block data                
                pySW4Rfile.write_properties(fileObject,vp.ravel(),5,vs.ravel(),p.ravel(),qp.ravel(),qs.ravel()) 
                #pySW4Rfile.write_properties(fileObject,vp.ravel(),5) 
                #next block 
                blockIndex += 1
                
        #make sure that you got everything
        #if the unit is above the topograophy set it equal to nothing!
        
        fileObject.flush()
        pass
    @staticmethod
    def interpolateUnit(self,x,y,z,X,Y,Z,blockData):
        #interpolate on each x,y palne and save all of the results accordingly
        return griddata(x.ravel(),X.ravel(),blockData)
        #for a first attempt just to a straight linear interpolation and see what you get
        
        #for i in len(Z):
                
        
                
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
    #TODO fix these two!
    @staticmethod
    def computeTop(z0,baseResolution,datum):
        if(z0<0):
            return (abs(z0)+datum)/baseResolution        
        else:
            return abs(z0-datum)/baseResolution
            pass
            
    @staticmethod
    def computeBottom(top,baseResolution,HVb,Nk):
        return top-((HVb*Nk)/baseResolution)


def main():
    outPutFileName = "GFM_SMALL.r"
    
    fileObject = open(outPutFileName,"wab")
    #topo test /home/eeckert/git/pySW4Forked/pySW4/utils/untrackedModels/NNSS_VOLUME/Topography.2grd.dat
    #self = Model("rFile.ini","untrackedModels/GFM_all_clean","")
    #actually do it this way when running for real
    #with open(outPutFile, "wab")
    
    
    pass
