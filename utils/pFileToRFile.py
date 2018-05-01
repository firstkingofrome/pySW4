"""
pFile to r file 
Generates an R file (with flat topography) from the given pFile with a constant depth delta
Eric Eckert, Nevada Seismology Lab
4/17/18
"""
import os
import sys
import pySW4 as sw4
import matplotlib.pyplot as plt
import numpy as np


#loads a pfile  for easy conversion to R file standard
class pFile():
    def __init__(self,fname="granite7Q_mix.8.400m.pfile"):   
        self.elementMap = {"lon":0,"lat":1,"depth":2,"CP":3,"CS":4,"p":5,"QP":6,"QS":7}

        dataSectionStart = 0
        #load the input
        self.loadPfile(fname)
        #condense to numpy object
        self.data =np.array(self.pFileData)  
        del self.pFileData #return some memory
        
        #stack it by depth
        self.points =np.dstack((self.data))[0]
        del self.data #return some more memory
        #and sort, actually never mind, a pfile is supposed to be sorted all ready
        
        #and return a little memory
        #del self.data
        pass
    
    def sortData(self):
        #note that I am sorting by actual depth, not the depth coordinates given!
        #this sorts the data file (using cartesian coordinates) so that the data will be arranged by (depth,minx,miny)
        #https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
        #note that this guy was very clever, he got the list indexes each time an then applied them to the data set!
        #I actually dont think that I need a stable sort come to think of it, (there should be no two depths that are identical)
        self.points = self.points[abs(self.points[self.elementMap["lon"]].argsort(kind='quicksort'))] # Y using a mergesort becuase I need the sort to be stable (i.e only one correct solution)
        self.points = self.points[abs(self.points[self.elementMap["lat"]].argsort(kind='quicksort'))] # X
        self.points = self.points[self.points[self.elementMap["depth"]].argsort(kind='quicksort')] # Depth sorted last (so it is the most significant ordering!)
    
    #clean this up to make it faster at a later point
    def loadPfile(self,fname):
        lat = 0.0
        lon = 0.0
        lines = []
        self.header = []
        self.pFileData = []
        with open( fname,'r') as fileObject:
            while True:
                #get the header
                line =  fileObject.readline()
                #save the line
                self.header.append(line)
                if(line[0] == 'T' or '.true.' in line):
                    break                    
            for line in fileObject:           
            #now we are done with the header
                if(len(line.split()) ==3):
                    #then this is a coordinate pair, discard the num points value
                    lines = [i for i in line.split(" ") if i != '']
                    lat = np.float32(lines[0])
                    lon = np.float32(lines[1])   
                                     
                else:
                    #the multiplication by 1000 moves from m/s to km/s
                    self.pFileData.append([])
                    lines = [i+" " for i in line.split() if i != ""]
                    #note that a pfile is in Km/s format instead of the M/s units used by r files
                    self.pFileData[-1].append(lon)
                    self.pFileData[-1].append(lat)
                    self.pFileData[-1].append(np.float32(lines[1])*1E3) #depth
                    self.pFileData[-1].append(np.float32(lines[2])*1E3) #Cp
                    self.pFileData[-1].append(np.float32(lines[3])*1E3) #Cs
                    self.pFileData[-1].append(np.float32(lines[4])*1E3) #P
                    self.pFileData[-1].append(np.float32(lines[5]))     #Qp
                    self.pFileData[-1].append(np.float32(lines[6]))     #Qs
       
#compute number of meters for degree delta (extremly crappy approximation)
def distanceDelta(deltaDegrees):
    #           Earth Circumference
    return round(40075.0517E3*(deltaDegrees /360.0),-1)

def depthDelta(depths,nDepths):
    return round(((((max(depths)-min(depths))))/nDepths),-1)
    


def write_properties(f,vp, nc, vs=None, rho=None, qp=None, qs=None):
    """
    Write material properties at a point `i`, `j` in block `b`.

    This is a convinient function to use while looping over `i`, `j` in
    a specific block `b` for writing out material properties at each
    index `k`. At the very least `vp` should be provided. If only `vp`
    is provided, the other properties are calculated using the
    :mod:`~..material_model` module.

    Parameters
    ----------
    f : file
        Open file handle in ``'wa'`` mode for appending data to the end
        of a file in construction ot in ``'r+b'`` mode for overwriting
        existing data.

        When overwriting existing data, the user must take care to place
        the cursor in the right place in the file.

    vp : array-like
        P wave velocity at indicies of `k` in m/s.

    nc : int
        Number of components to write out. Either 3 (`rho`, `vp`, and
        `vs` if ``attenuation=0``) or 5 (also `qp` and `qs` if
        ``attenuation=1``).

    vs : array-like, optional
        S wave velocity at indicies of `k` in m/s. If not given, `vs` is
        calculated from :func:`~..material_model.get_vs`.

    rho : array-like, optional
        Density at indicies of `k` in kg/m^3. If not given, `rho` is
        calculated from :func:`~..material_model.get_rho`.

    qp : array-like, optional
        P quality factor at indicies of `k`. If not given, `qp` is
        calculated from :func:`~..material_model.get_qp`.

    qs : array-like, optional
        S quality factor at indicies of `k`. If not given, `qs` is
        calculated from :func:`~..material_model.get_qs`.
    """

    k_array = np.empty((vp.size, nc), np.float32)
    vs = vs
    rho = rho 
    qs = qs 
    qp = qp 

    k_array[:, 0] = rho 
    k_array[:, 1] = vp 
    k_array[:, 2] = vs 
    k_array[:, 3] = qp
    k_array[:, 4] = qs
    k_array.tofile(f)

#TODO figure out a clever way to dynamically generate resolutions (based on the maximums availbile)
def main():

    topography = ""
    inPutFile = "granite7Q_mix.8.400m.pfile"
    outPutFile = "ericRFile.r"
            #   max hieght    min hieght, unique x,unique y, grid spacing horrizontal (meters), gridSpace vertical meters)
    projectionString = ('+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs')
    az = 0.0
    DistanceDelta=0
    DepthDelta = 0
    #I am using attentuation
    attenuation = 1
    components = 5
    #assign lat and lon 0 form p file
    lat0 =0.0
    lon0 = 0.0
    #load the file and sort the data
    data = pFile(inPutFile)
    DistanceDelta = distanceDelta(float(data.header[1]))
    DepthDelta = depthDelta(data.points[2],int(data.header[4].split()[0]))
    lat0= min(data.points[data.elementMap["lat"]])
    lon0 = min(data.points[data.elementMap["lon"]])
    #prepare the block strucutre, since this is a pfile, I assume that they only want two blocks
    #           top of blocl        bottom of block (z0) i,j,k, HH,HV
    BLOCKS = [[min(data.points[2])*1E3,max(data.points[2])*1E3,len(np.unique(data.points[1])),len(np.unique(data.points[0])),len(np.unique(data.points[2])),DistanceDelta,DepthDelta],
    [min(data.points[2])*1E3,max(data.points[2])*1E3,len(np.unique(data.points[1])),len(np.unique(data.points[0])),len(np.unique(data.points[2])),DistanceDelta,DepthDelta]]   
    BLOCK_MAPPING = {"BOTTOM":0,"TOP":1,"i":2,"j":3,"k":4,"HH":5,"HV":6}   
    pointRanges = []
    depthPoints = []
    horizontalPoints = []
    #fill sw4 structures
        #open file object
    with open(outPutFile, "wab") as rFile:
            #write file header header, rFile = open("DAGRfile.r", "wab")
            sw4.prep.rfileIO.write_hdr(rFile, 1, 4, attenuation, az,lon0, lat0, projectionString.encode(),len(BLOCKS))
            #write_hdr(rFile, 1, 4, attenuation, az,lon0, lat0, projectionString,len(BLOCKS)+1)
            #sw4.prep.rfileIO.write_hdr(rFile, 1, 4, attenuation, az,lon0, lat0, projectionString,len(BLOCKS)+1)
            #write block headers                
            #first write the topography header
            #write_block_hdr(f, hh, hv, z0, nc, ni, nj, nk)            
            sw4.prep.rfileIO.write_block_hdr(rFile,BLOCKS[0][5], BLOCKS[0][6], BLOCKS[0][0], 1, BLOCKS[0][2], BLOCKS[0][3], 1.0)
            #now write the block headers
            for i in range(1,len(BLOCKS)):
                #compute the best possible resolution for each block (based on the depth range)
                #see if the i,j,k # points make sense or not                
                sw4.prep.rfileIO.write_block_hdr(rFile,BLOCKS[i][5], BLOCKS[i][6], BLOCKS[i][0], components, BLOCKS[i][2], BLOCKS[i][3], BLOCKS[i][4])                                                
            #topo data block
            if(topography == ""):
                #assign to zero (flat)
                sw4.prep.rfileIO.write_topo_block(rFile,np.zeros((BLOCKS[0][2],BLOCKS[0][3]),dtype=np.float32))
            else:
                #write from the description in the topograophy file
                pass
                                                                                       
            for i in range(len(data.points[1])):          
                
                #sw4.prep.rfileIO.write_properties(rFile, data.points[data.elementMap["CP"]][i], components)
                #note that you dont have to do this in steps of 1 as you have (dispense with the for loop and do it in bulk later)
                write_properties(rFile, data.points[data.elementMap["CP"]][i], components,data.points[data.elementMap["CS"]][i],
                data.points[data.elementMap["p"]][i], data.points[data.elementMap["QP"]][i],
                data.points[data.elementMap["QS"]][i])                                       

    #done
    pass
    
if __name__=="__main__":
    #main()
    pass
