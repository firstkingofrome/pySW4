"""
process all sw4 images in the specified directory
Generates an R file (with flat topography) from the given pFile
Eric Eckert, Nevada Seismology Lab
4/17/18
"""


import os
import sys

import pySW4 as sw4
import obspy
import matplotlib.pyplot as plt
import numpy as np

from argparse import ArgumentParser 


#incomplete
def genImage(f,vmax='3',cmap='jet'):
    image = sw4.read_image(f)
    pass

#build movie
#include input file from sw4 run to populate with metadata
def buildMovie(imageDir="/home/eeckert/git/rFileArben/results/myRfileHigherTimestep",inputFile = None, outPut="movie.mp4",fps=1,movieType = 'velocity'):
    #set the path
    os.chdir(imageDir)
    #get all of the sw4 images in the directory
    files = os.listdir('.')
    #only consume things that are actually images
    #files = [i for i in files if(".sw4img" in i)]
    #look for the input file
    if([i for i in files if(".in" in  i)] != []):
        inputFile = [i for i in files if(".in" in  i)][0]
    #digest
    sw4.plotting.image.create_image_plots(inputFile, imageDir, movieType, fps, None, True)


def main():


    return 0
    
if __name__=="__main__":
    buildMovie()
    #main()
    pass

