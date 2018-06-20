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
import pySW4.utils
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
            self.topoData = np.asarray(self.topoData)
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
        #underlying coordinates
        self.x,self.y,self.z = [],[],[]
        
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
                z = int(abs(float(lines[2]) - self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"]) if float(lines[2]) > 0 else abs(float(lines[2])) + self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"])
                #fix everything according to the datums           
                #also internally I assign 0 to be the bottome of the rfile!b                             
                self.ModelFileData.append((float(lines[0]),float(lines[1]),z,int(lines[3])))
                
        #convert all of the coordinates to the Carteisian coordinates used by sw4
        self.ModelFileData = np.asarray(self.ModelFileData)
        #note that Shhar use lon lat instead of lat lon
        #I would change my stuff except that I am so behind on this project that I dont want to sacrafice my limited time to that end
        self.ModelFileData[:,1],self.ModelFileData[:,0] = pySW4.utils.simple_lonlat2xy(self.ModelFileData[:,1],self.ModelFileData[:,0],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["LON0"],
        self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["LAT0"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["AZIMUTH"])
        #now sort it so that it sits in memory in the same way that the rfile does
        self.ModelFileData = self.ModelFileData[abs(self.ModelFileData[:,0].argsort(kind='quicksort'))]
        self.ModelFileData = self.ModelFileData[abs(self.ModelFileData[:,1].argsort(kind='quicksort'))]
        self.ModelFileData = self.ModelFileData[abs(self.ModelFileData[:,2].argsort(kind='quicksort'))]
       
        """
        for i in range(len(self.ModelFileData)):
            self.ModelFileData[i][2] = 
        """
        
        """
        I think that this is the problem , you need to decide youre own dimensions and then impose onto that!, this will distory youre model substantially!
        """            
        self.x,self.y,self.z=np.mgrid[0:self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseNi"],0:self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseNj"] ,0:self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["baseNk"]]        
        #self.x,self.y,self.z=np.mgrid[0:len(np.unique(self.ModelFileData[:,0])),0:len(np.unique(self.ModelFileData[:,1])) ,0:len(np.unique(self.ModelFileData[:,2]))]
        #assign data to "nearest" point on the mesh
        self.inputModel = 0
        
    #fix coordinates (i.e use cartesian coordinates)
    def fixCoordinates(self):
        #find the minimum value for all of the UTM coordinates and assign accordingly
        pass
        
    def buildRfile(self,fileObject):   
        vp,Vip=[],[]
        vs,Vis=[],[]
        qp,Qip=[],[]
        qs,Qis=[],[]
        p,Pis=[],[]
        topo = []
        componentMesh = []
        blockIndex = 1  
        baseDepth = 0
        topBlock = 0
        blockRange = 0
        #precompute material properties for all of my points        
        #apparently this actually takes longer than doing the easy way np = np.array([self.getP(i) for i in np.nditer(self.ModelFileData[:,3])])
        #you should still think about this a little and see if there is some clever row op thing that you can do to speed this up (as Shhar does to compute coordinates)
        vp = np.asarray([self.getVP(self.ModelFileData[i][3]) for i in range(len(self.ModelFileData[:,3]))])
        #vp = griddata((x.ravel(),y.ravel(),z.ravel()),blockData,(X.ravel(),Y.ravel(),Z.ravel()),method='nearest')

        #vs
        vs = np.asarray([self.getVS(self.ModelFileData[i][3]) for i in range(len(self.ModelFileData[:,3]))])
        #vs = griddata((x.ravel(),y.ravel(),z.ravel()),blockData,(X,Y,Z),method='nearest')
        
        #p
        p = np.asarray([self.getP(self.ModelFileData[i][3]) for i in range(len(self.ModelFileData[:,3]))])
                        
        #qp
        qp = np.asarray(np.asarray([self.getQP(self.ModelFileData[i][3]) for i in range(len(self.ModelFileData[:,3]))]))
                        
        #qs
        qs = np.asarray(np.asarray([self.getQS(self.ModelFileData[i][3]) for i in range(len(self.ModelFileData[:,3]))]))
        
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
            #use the topo loaded from the file
            X,Y = np.meshgrid(np.arange(0, self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["Ni"]),np.arange(0, self.Parameterfile.pfContents["BLOCK_CONTROL"]["TOPO"]["Nj"])) 
            x,y = np.meshgrid(np.arange(self.inputModel.shape[0]),np.arange(self.inputModel.shape[1]))                           

            vp = griddata((x.ravel(),y.ravel()),self.topoData.ravel(),(X.ravel(),Y.ravel()),method='nearest')

            pySW4Rfile.write_topo_block(fileObject,vp)
            print("wrote topogaphy!")
        
        #now save the block data
        for block in self.Parameterfile.pfContents["BLOCK_CONTROL"]:
            #make sure that I am dealing with an actualy BLOCK
            if(block != "TOPO" and block != "HEADER"):                       
                #compute the depth datum--using my modified coordinates
                topBlock = self.computeTop(self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Z0"],self.Parameterfile.pfContents["BLOCK_CONTROL"]["HEADER"]["maxDepthBelowSeaLvl"])
                baseDepth = self.computeBottom(topBlock,self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"])
                #use these to pull the correct part of the block dimensions, then interpolate and throw this in the rFile
                blockRange = np.where(self.ModelFileData[topBlock - self.ModelFileData[:,2] >= baseDepth])[0]
                coords = self.ModelFileData[blockRange][:,0:3]
                #data = vp[blockRange]
                HHb = self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HHb"]
                HVb = self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["HVb"]
                shape = (self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Ni"],self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nj"],
                self.Parameterfile.pfContents["BLOCK_CONTROL"][block]["Nk"])
                #compute the data

                """
                coords,data, shape,HHb,HVb, method='linear',
                draft=False, corners=None, origin='nw', verbose=False
                """
                Vip = resample3D(coords,vp[blockRange].ravel(),shape,HHb,HVb,"nearest").ravel()
                
                pySW4Rfile.write_properties(fileObject,Vip,5,resample3D(coords,vs[blockRange],shape,HHb,HVb,"nearest").ravel(),
                resample3D(coords,p[blockRange].ravel(),shape,HHb,HVb,"nearest").ravel(),resample3D(coords,qp[blockRange],shape,HHb,HVb,"nearest").ravel(),resample3D(coords,qs[blockRange],shape,HHb,HVb,"nearest").ravel()) 
                
                #pySW4Rfile.write_properties(fileObject,resample3D(coords,data,HHb,HVb,"nearest").ravel(),5,vs.ravel(),p.ravel(),qp.ravel(),qs.ravel()) 
                #pySW4Rfile.write_properties(fileObject,vp.ravel(),5) 
                #next block 
                    
                blockIndex += 1
                
        #make sure that you got everything
        #if the unit is above the topograophy set it equal to nothing!        
        fileObject.flush()
        pass
    
    #TODO implement full support for the direction opperator which instructs the machine how to splice the interpolation (wish I draw picutre in a comment)
    #the basic functionality is that each x y and z holds all of the planes in the volume, spliced in a particular direction, this controls the directionality.
    @staticmethod
    def interpolateUnit(x,y,z,X,Y,Z,blockData,topBlock,baseDepth,direction="z"):
        planeData = []
        dataMin = 0
        dataMax = 0
        
        #to recap the problem here is that blockData is in a different direction than x or y, so you have to transpose it to correct it!
        if(direction=='z' and (direction != 'x' or direction != 'y')):
            for i in range(baseDepth,topBlock):
                #compute the range of block data for each plane that I am dealing with
                dataMin = 0
                dataMax = 0
                print(i)
                #I am going to have to recompute X and Y for each plane if I decide to do this with this process

                test = griddata((self.x[i].ravel(),self.y[i].ravel()),blockData[i].T.ravel(),(X[i],Y[i])) 
                planeData.append()
            
            return 
            pass
        
        #add in the other direction opperations
        
        
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
    def computeTop(z0,datum):
        if(z0<0):
            return (abs(z0)+datum)      
        else:
            return abs(z0-datum)
            pass
            
    @staticmethod
    def computeBottom(top,HVb,Nk):
        return top-((HVb*Nk))
    
    #Tri linear interpolation from  https://stackoverflow.com/questions/6427276/3d-interpolation-of-numpy-arrays-without-scipy

def resample3D(coords,data, shape,HHb,HVb, method='linear',
             draft=False, verbose=True):
    """
    Resample or interpolate an array on to a (3d) grid with new extent and or
    new shape.
    #todo figure out how to do this with openCV or some other (faster) platform
    """
    
    if coords.shape[1] != 3:
        print('Error: data must be a tuple of 3x2d '
              ':class:`numpy.ndarry` (X, Y, Z) or one 2d ``Z`` array')
              

    """
    if X is None and Y is None:
        X, Y = np.meshgrid(np.linspace(w, e, z_nx),
                           np.linspace(s, n, z_ny))
    """
    if draft is False:  # the long and accurate way...
        if verbose:
            message = ('Accurately interpolating data onto a grid of '
                       'shape %d%d%d and '
                       'using X, Y, Z arrays.\n'
                       '                      This may take a while...')
            print(message)
            sys.stdout.flush()

        xi, yi,zi = np.meshgrid(np.linspace(0, shape[0]*HHb,shape[0]),
                             np.linspace(0, shape[1]*HHb,shape[1]),np.linspace(0, shape[2]*HVb,shape[2]))
        di = griddata((coords[:,0].ravel(), coords[:,1].ravel(), coords[:,2].ravel()),data.ravel(), (xi, yi,zi),
                      method=method)
        return di

    elif draft is True:  # the fast and less accurate way...
        nx =0
        ny=0
        try:  # both src and dst are passed
            src, dst = corners
            src, dst = tuple(src), tuple(dst)

        except ValueError:  # only dst corners passed
            src = ((0, 0), (nx, 0), (nx, ny), (0, ny))
            dst = corners

        xc, yc = [p[0] for p in dst], [p[1] for p in dst]
        xc, yc = xy2pixel_coordinates(xc, yc, extent, shape, origin)
        dst = tuple(zip(xc, yc))

        if verbose:
            message = ('Transforming points:\n'
                       '%s\n'
                       'in the data to points:\n'
                       '%s\n'
                       'in the the output grid of shape %dx%d and '
                       'extent %.2f,%.2f,%.2f,%.2f.')
            print(message.format(src, dst, nx, ny, w, e, s, n))
            sys.stdout.flush()

        # Compute the transformation matrix which places
        # the corners Z at the corners points bounding the
        # data in output grid pixel coordinates
        tranformation_matrix = cv2.getPerspectiveTransform(np.float32(src),
                                                           np.float32(dst))

        # Make the transformation
        interpolation = {'nearest' : 0,
                         'linear'  : 1,
                         'cubic'   : 2}
        zi = cv2.warpPerspective(Z, tranformation_matrix,
                                 (nx, ny),
                                 flags=interpolation[method],
                                 borderMode=0,
                                 borderValue=np.nan)

        return zi


"""
try doing it this way
width = 200
height = 200
img_stack_sm = np.zeros((len(img_stack), width, height))

for idx in range(len(img_stack)):
    img = img_stack[idx, :, :]
    img_sm = cv2.resize(img, (width, height), interpolation=cv2.INTER_CUBIC)
    img_stack_sm[idx, :, :] = img_sm


"""



def main():
    outPutFileName = "GFM_SMALL.r"
    
    fileObject = open(outPutFileName,"wab")
    #topo test /home/eeckert/git/pySW4Forked/pySW4/utils/untrackedModels/NNSS_VOLUME/Topography.2grd.dat
    #self = Model("rFile.ini","untrackedModels/GFM_all_clean","")
    self = Model("rFile.ini","untrackedModels/GFM_wgs84.txt","")
    #actually do it this way when running for real
    #with open(outPutFile, "wab")
    self.buildRfile(fileObject)

main()


