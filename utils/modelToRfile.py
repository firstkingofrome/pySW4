import os
import sys
import pySW4 as sw4
import matplotlib.pyplot as plt
import numpy as np
import json


class parameterFile():
    def __init__(self,pFname = "rFile.ini"):
        self.pfContents = []            
        #open the parameterFile and load its contents
        self.loadPF(pFname)

    
    def loadPF(self,pFname):
        with open(pFname,'r') as fileObject:
            self.pfContents = [line.strip().encode('string-escape') for line in fileObject if("#" not in line.strip())]
        #attempt to jsonify
        #r is to indicate that this is a raw string (i.e ignore \ escape charecter
        self.pfContents = json.loads(r''.join(self.pfContents))


class geologicModel():
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

def main():
    #allocate pixel space
    
    #load the geologic model and map each point to pixel space (at the same time)
    
    #interpolate
    
    #
    
    #save pixel space to the rFile
    pass

if __name__=="__main__":
    main()
    pass


