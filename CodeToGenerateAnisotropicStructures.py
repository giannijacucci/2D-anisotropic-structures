from matplotlib import pyplot as plt
from matplotlib import interactive
from matplotlib.pyplot import figure, show
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry
from shapely.geometry.point import Point
from shapely import affinity
from matplotlib.patches import Polygon
import numpy as np
import pylab as pl
import time
from termcolor import colored, cprint
from scipy.stats import norm
from datetime import datetime
import os

def main():
    bbox=(-2.5,2.5) #sample size (µm)
    step=7 #number of, statistical identical, random structures to generate
    alphaMean=0 #particles' orientation
    alphaSigma=0 #deviation from the average orientation
    
    if index<(indexMin+step):
        c=7 #aspect ratio
        r=0.1 #radius of the corresponding (same area) isotropic particle (in µm)
        rsigma=r*30/100 #polydispersity
        ff=0.60 #filling fraction

    ellipses=[]
    currentFF=0
    particlesToAdd=[]

    while currentFF<ff:
        rdd=normDist(r,rsigma)
        particlesToAdd.append(rdd)
        currentFF=currentFF+(rdd*rdd*c*3.1427)/(((bbox[1]-bbox[0])**2))

    
    particlesToAdd=np.array(particlesToAdd)
    particlesToAdd=np.flip(np.sort(particlesToAdd),0) #reodered to start adding from the largest particles
    aOccupied=0

    k=0
    iteration=0
    resampled=0 #checks for eventual size-rejection
    aOccupied=0
    while k<len(particlesToAdd):
        while(iteration < 20000):
            iteration += 1
            rd=particlesToAdd[k]
            alpha = normDist(alphaMean, alphaSigma)
            center = boxForcing(bbox[0], bbox[1], rd, c, alpha)
            x = center[0][0]
            y = center[0][1]
            ellips = create_ellipse((x, y), rd, (1, c), alpha)
            overlap = pl.array([ellips.intersects(ellipses[i][4]) for i in range(len(ellipses))])

            if all(overlap[i] == False for i in range(len(ellipses))):
                one_row = []
                one_row.extend([x, y, rd, alpha, ellips])
                ellipses.append(one_row)
                aOccupied +=rd*rd*c*3.1427
                cprint(aOccupied/((bbox[1]-bbox[0])**2), 'green')
                k+=1
                iteration=0
                break

            if iteration == 20000:
                cprint('rejected', 'yellow')
                resampled+=1
                break
        else:
            iteration=0
            rd = normDistRejected(r,rd, rsigma) #samples a smaller particle


    radDist=[]
    for i in range(len(ellipses)):
        radDist.append(ellipses[i][2]*2*c)
   
    (mu, sigma) = norm.fit(radDist)  
    
    if 0: #visualize the size distribution
        cprint("Plotting legth distribution","cyan")
        
        (mu, sigma) = norm.fit(radDist)
        xmin, xmax = (mu-5*sigma,mu+5*sigma)
        x = pl.linspace(xmin, xmax, 1000)
        pl.hist(radDist,bins=15,normed=True)
        p = norm.pdf(x, mu, sigma)
        pl.plot(x, p, 'k', linewidth=2)
        pl.xlabel('Length (µm)')
        pl.ylabel('Probability')
        pl.title(r'$\mathrm{Radii\ distribution\:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
        pl.savefig("RadiusDist"+str(index)+".pdf")
        pl.close('all')
        pl.xlim([0,0.4])    
    
    f = open("2D"+str(index)+".txt",'w')
    f.write("x-pos\ty-pos\ta\tb\talpha\n")
    for elm in ellipses:
        f.write("{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n".format(elm[0],elm[1],elm[2],elm[2]*c,elm[3]))
    f.close()
    
    if 1: #visualise the final structure
        vertices=pl.array([ellipses[i][4].exterior.coords.xy for i in range(len(ellipses))])
        patches=pl.array([Polygon(vertices[i].T,color = 'black', alpha = 1) for i in range(len(ellipses))])
        fig,ax = plt.subplots()
        ax.set_xlim([bbox[0],bbox[1]])

        ax.set_ylim([bbox[0],bbox[1]])
        ax.set_aspect('equal')
    
        for i in range(len(patches)):
             ax.add_patch(patches[i])

        pl.title("Distribution:("+str(mu)[0:5]+", "+str(sigma)[0:5]+"); Filling Fraction:"+str(ff))
        pl.savefig("Ensemble"+str(index)+".pdf")
        pl.close('all')

def normDist(r,rsigma): #defines an always positive normal distribution (for very broad size distributions)
	x=pl.normal(r,rsigma)
	while x<0:
		x=pl.normal(r,rsigma)
	return x

def normDistRejected(r,rd,rsigma): #looks for smaller particles in case of rejection
	x=pl.normal(r,rsigma)
	while x>rd:
		x=pl.normal(r,rsigma)
	return x

def create_ellipse(center,radius,lengths,angle=0): #creates an ellipse using Shapely
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    
    circ = Point(center).buffer(radius)
    ell = affinity.scale(circ, int(lengths[0]), int(lengths[1]))
    ellr = affinity.rotate(ell, angle)
    return ellr

def boxForcing(bMin,bMax,radius,c,angle): #confines the particles inside the sample
    
    x=pl.uniform(bMin,bMax)
    y=pl.uniform(bMin,bMax)
    
    if x<(bMin+c*radius):
        x=x+(radius)*(1+c*pl.absolute(pl.sin(pl.radians(angle))))
        if x>(bMax-c*radius):
            x=bMin+c*radius
    if x>(bMax-c*radius):
        x=x-(radius)*(1+c*pl.absolute(pl.sin (pl.radians(angle))))

        if x<(bMin+c*radius) :
            x=bMax-c*radius

    if y>(bMax-c*radius): 
        y=y-(radius)*(1+c*pl.absolute(pl.cos(pl.radians(angle))))
        if y<(bMin+c*radius) :
            y=bMax-c*radius

    if y<(bMin+c*radius):
        y=y+(radius)*(1+c*pl.absolute(pl.cos(pl.radians(angle))))
        if y>(bMax-c*radius):
            y=bMin+c*radius

    center=[]
    center.append([x,y])

    return center
    
indexMin=1
indext=1
indexMax=2
for i in range(indext,indexMax):
		index=i
		main()
