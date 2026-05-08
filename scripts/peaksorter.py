# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 16:49:35 2024

@author: lwg513
"""
from tkinter import filedialog
import sys
#from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.ticker as mtick
import matplotlib

class peaksort():
    
    def roundup(x,r):
        return int(math.ceil(x/r))*r
    
    def mass_jitter(refmass,refjitter,ppm,mass):
        constant=(refjitter-ppm*refmass*1e-6)/refmass**2
        PDMT=(constant*mass**2)+(ppm*mass*1e-6)
        return PDMT
    
    def TwoDMS_plot(X,Y,Z):
        revY=Y.iloc[::-1]
        revX=X.iloc[::-1]#reversed as want the most intense plotted last
        revZ=Z.iloc[::-1]
        fig, ax=plt.subplots(dpi=75)
        cm=matplotlib.colormaps.get_cmap('RdPu')
        sc = ax.scatter(revX, revY, c=np.log(revZ), s=0.5, cmap=cm)
        ax.set_ylim([min(Y), max(Y)])
        ax.set_xlim([min(X), max(X)])
        ax.set_xlabel("Fragment m/z", size =20,style="italic")
        ax.set_ylabel("Precursor m/z", size=20, style="italic")
        plt.colorbar(sc, ax=ax, label="log10 intensity", spacing='proportional') #log10 to make lower intensiites more visible
        fig.set_figheight(12)
        fig.set_figwidth(14)
        return fig
    
    def TwoDMS_plot_neutral_loss(X,Y,Z):
        neutrallossy=Y-X
        neutrallossx=Y
#        for i in range (0,len(X)):
#            tmpNEUTRALLOSSX=Y[i]-X[i]
#            neutrallossx.append(tmpNEUTRALLOSSX)
#        neutrallossx=np.array(neutrallossx)
        revneuy=neutrallossy.iloc[::-1]
        revneux=neutrallossx.iloc[::-1]
        revZ=Z.iloc[::-1]
        fig, ax=plt.subplots(dpi=75)
        cm = plt.cm.get_cmap('RdPu')
        sc = ax.scatter(revneux, revneuy, c=np.log(revZ), s=0.5, cmap=cm)
        ax.set_ylim([0, max(neutrallossy)])
        ax.set_xlim([min(neutrallossx), max(neutrallossx)])
        ax.set_xlabel("Precursor m/z", size =20, style="italic")
        ax.set_ylabel("Neutral loss Da", size=20)
        plt.colorbar(sc, ax=ax, label="log10 intensity", spacing='proportional')
        fig.set_figheight(12)
        fig.set_figwidth(14)
        return fig
        
    def msms_peaklister(data,index):
        if not data:
            print("must run peak sorter first")
            exit()
        peaklist=data[index]
        peakstmp=[]
        xvalues=[]
        yvalues=[]
        zvalues=[]
        for i in range(1,len(peaklist)):
            peakstmp.append(peaklist[i])
        for i in range (0,len(peakstmp)):
            xvalues.append(peakstmp[i][0])
            yvalues.append(peakstmp[i][1])
            zvalues.append(peakstmp[i][2])
        return xvalues, yvalues, zvalues
    
    def msms_plotter(xvalues, yvalues, zvalues,header, lowmass, highmass, roundto=1e8):
#        ymax=peaksort.roundup(max(zvalues),roundto)
#        plt.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
#        plt.xlabel("m/z", size =15,style="italic")
#        plt.ylabel("Intensity", size=15)
#        plt.title(str(header), size=20)
#        plt.xticks(np.arange(peaksort.roundup(lowmass,xstep), highmass, step=xstep), fontsize=15) 
#        plt.yticks(fontsize = 15) 
#        plt.ylim([0, ymax])
#        plt.figure(figsize=(150,60))
        header=round(header,6)
        y_max=peaksort.roundup(max(zvalues),roundto)
        x_max=peaksort.roundup(max(xvalues),50)
        fig, ax=plt.subplots(dpi=75)
        ax.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
        ax.set_xlabel("m/z", size =20,style="italic")
        ax.set_ylabel("Intensity", size=20)
        ax.set_title(str(header), size=20)
#        ax.set_xticks(np.arange(lowmass, xmax, step=50), fontsize=15) 
#        ax.set_yticks(fontsize = 15) 
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(8)
        fig.set_figwidth(12)
        return fig
    
    def NL_plotter(xvalues, yvalues, zvalues,header, lowmass, highmass, roundto=1e8):
        header=round(header,6)
        y_max=peaksort.roundup(max(zvalues),roundto)
        x_max=peaksort.roundup(max(xvalues),50)
        fig, ax=plt.subplots(dpi=75)
        ax.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
        ax.set_xlabel("m/z", size =20,style="italic")
        ax.set_ylabel("Intensity", size=20)
        ax.set_title("Neutral loss of "+str(header), size=20) 
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(8)
        fig.set_figwidth(12)
        return fig
    
    def VL_plotter(xvalues, yvalues, zvalues,header, lowmass, highmass, roundto=1e8):
        header=round(header,6)
        y_max=peaksort.roundup(max(zvalues),roundto)
        x_max=peaksort.roundup(max(xvalues),50)
        fig, ax=plt.subplots(dpi=75)
        ax.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
        ax.set_xlabel("m/z", size =20,style="italic")
        ax.set_ylabel("Intensity", size=20)
        ax.set_title("Vertical line at "+str(header), size=20) 
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(8)
        fig.set_figwidth(12)
        return fig
    
    def isotope_msms_peaklister(data,indexes,jitter):
        peaklist=[]
        xvalues=[]
        yvalues=[]
        zvalues=[]
        for index in indexes:
            peaklist=data[index]
            tmppeaklist=peaklist[1:]
            for i in range (0,len(tmppeaklist)):
                xvalues.append(tmppeaklist[i][0])
                yvalues.append(tmppeaklist[i][1])
                zvalues.append(tmppeaklist[i][2])
        for i in xvalues:
            for j in xvalues:
                if i==j:
                    continue
                elif abs(i-j)<= jitter[xvalues.index(i)]:
                    del yvalues[xvalues.index(j)]
                    zvalues[xvalues.index(i)]=zvalues[xvalues.index(i)]+zvalues[xvalues.index(j)]
                    del zvalues[xvalues.index(j)]
                    del xvalues[xvalues.index(j)]
                else:
                    continue
        return xvalues, yvalues, zvalues
                
    def isotope_msms_plotter(xvalues, yvalues, zvalues, lowmass, highmass, indexes, roundto=1e8):
        y_max=peaksort.roundup(max(zvalues),roundto)
        x_max=peaksort.roundup(max(xvalues),50)
        fig, ax=plt.subplots(dpi=75)

        ax.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
        ax.set_xlabel("m/z", size =20,style="italic")
        ax.set_ylabel("Intensity", size=20)
#        ax.set_xticks(np.arange(lowmass, xmax, step=50), fontsize=15) 
#        ax.set_yticks(fontsize = 15) 
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(8)
        fig.set_figwidth(12)
        return fig
    
    
    def identify_ridges(xvalues, yvalues, zvalues, header, precursorsx, refmass, refjitter, ppm=0):
        ridgesx=[]
        ridgesz=[]
        NOPX=[]
        NOPY=[]
        NOPZ=[]
        jitter=[]
        for i in range(0,len(yvalues)):
            if yvalues[i]!=header:
                NOPX.append(xvalues[i])
                NOPY.append(yvalues[i])
                NOPZ.append(zvalues[i])
        for i in range(0,len(NOPX)):
            for j in range(0,len(precursorsx)):
                tmpjitter=peaksort.mass_jitter(refmass,refjitter,ppm,precursorsx[j])
                jitter.append(tmpjitter)
                if NOPX[i]<(precursorsx[j]+jitter[j]) and NOPX[i]>(precursorsx[j]-jitter[j]):
                    ridgesx.append(NOPX[i])
                    ridgesz.append(NOPZ[i])
        return ridgesx,ridgesz
        
        
    
    def highlight_msms_plotter(xvalues, yvalues, zvalues, header, lowmass, highmass, ridges, roundto=1e8):
#        ymax=peaksort.roundup(max(zvalues),roundto)
#        plt.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
#        plt.xlabel("m/z", size =15,style="italic")
#        plt.ylabel("Intensity", size=15)
#        plt.title(str(header), size=20)
#        plt.xticks(np.arange(peaksort.roundup(lowmass,xstep), highmass, step=xstep), fontsize=15) 
#        plt.yticks(fontsize = 15) 
#        plt.ylim([0, ymax])
#        plt.figure(figsize=(150,60))
        y_max=peaksort.roundup(max(zvalues),roundto)
        x_max=peaksort.roundup(max(xvalues),50)
        fig, ax=plt.subplots(dpi=75)
        ax.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
        ax.vlines(x = ridges[0], ymin = 0, ymax = ridges[1], colors="#FF10F0")
        ax.set_xlabel("m/z", size =20,style="italic")
        ax.set_ylabel("Intensity", size=20)
        ax.set_title(str(header), size=20)
#        ax.set_xticks(np.arange(lowmass, xmax, step=50), fontsize=15) 
#        ax.set_yticks(fontsize = 15) 
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(8)
        fig.set_figwidth(12)
        return fig

    
    def identify_precursors(fragmentmass, precursormass, intensity, refmass, reffwhm, ppm, refxy):
        precursors=[]  
        precursorintensity=[] 
        precursorjitter=[]
        xyresolution=[]
        for i in range (0,len(precursormass)):
            tempxy=peaksort.mass_jitter(refmass, refxy, ppm, precursormass[i])
            xyresolution.append(tempxy)
            if fragmentmass[i] < (precursormass[i]+tempxy) and fragmentmass[i] > (precursormass[i]-tempxy):
                precursors.append(precursormass[i])
                precursorintensity.append(intensity[i])
        for i in range (0,len(precursors)):
            tempjitter=peaksort.mass_jitter(refmass,reffwhm,ppm,precursors[i])
            precursorjitter.append(tempjitter)
        return precursors,precursorintensity,precursorjitter
    
    def precursor_plotter(precursors,precursorintensity, lowmass, highmass, roundto=1e8):
        y_max=peaksort.roundup(max(precursorintensity),roundto)
        x_max=peaksort.roundup(max(precursors),50)
        fig, ax=plt.subplots(dpi=75)
        ax.vlines(x = precursors, ymin = 0, ymax = precursorintensity, colors="black")
        ax.set_xlabel("m/z", size =20,style="italic")
        ax.set_ylabel("Intensity", size=20)
        ax.set_title("Precursors", size=20)
#        ax.set_xticks(np.arange(lowmass, xmax, step=50), fontsize=15) 
#        ax.set_yticks(fontsize = 15) 
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(8)
        fig.set_figwidth(12)
        return fig

    def mini_precursor_plotter(index,precursors,precursorintensity, lowmass, highmass, roundto=1e8):
        y_max=peaksort.roundup(max(precursorintensity),roundto)
        x_max=peaksort.roundup(max(precursors),50)
        fig, ax=plt.subplots(dpi=75)
        ax.vlines(x = precursors, ymin = 0, ymax = precursorintensity, colors="black")
        ax.vlines(x = precursors[index], ymin = 0, ymax = precursorintensity[index], colors="#FF10F0")
        ax.set_xlabel("m/z", size =10,style="italic")
        ax.set_ylabel("Intensity", size=10)
        ax.set_title("Precursors", size=10)
#        ax.set_xticks(np.arange(lowmass, xmax, step=50), fontsize=15) 
#        ax.set_yticks(fontsize = 15) 
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(3)
        fig.set_figwidth(5.25)
        return fig
        
        
    def peak_sort(precursors,fragmentmass,precursormass,intensity, PDMT):
        assigned=[]
        data=[]
        for i in range (0,len(precursors)): #cycles through the avilable precursor values
            tmp=[[precursors[i]]] #sets the precursor as the 'header'
            for j in range (0,len(fragmentmass)): #cycles through fragments (x values)
                if precursormass[j] >= (precursors[i]-PDMT[i]) and precursormass[j] <= (precursors[i]+PDMT[i]): 
                    #checks if in list of avilable precursor values,
                    assigned.append(fragmentmass[j])
                    tmp.append([fragmentmass[j],precursormass[j],intensity[j]])# if it is then adds to the tmp file 
                    #to contains the peak list for each precursor
            data.append(tmp)#adds to a large overall data set with i rows containing the x and z values
        unassigned=[]
        for i in range (0,len(fragmentmass)):
            if fragmentmass[i] not in assigned:
                unassigned.append([fragmentmass[i],precursormass[i],intensity[i]])
        return data, unassigned
    
    def mgf_text(data,inputname):
        mgftext=""
        for i in range (0,len(data)):
            mgftext=mgftext+"BEGIN IONS"
            mgftext=mgftext+"\nTITLE="+str(inputname)+" "+str(data[i][0]).strip("[]")+"scan=1 msLevel=2"
            mgftext=mgftext+"\nCHARGE=1+"
            mgftext=mgftext+"\nPEPMASS="+str(data[i][0]).strip("[]")
            for j in range (1,len(data[i])):
                mgftext=mgftext+"\n"+str(data[i][j][0])+" "+str(data[i][j][2])
            mgftext=mgftext+"\n\nEND IONS \n\n"
        return mgftext
    
    def sirius_text(data,compound_name,charge,indexes,xygap,xvalues, yvalues, zvalues,frag_type):
        ms1_data=[]
        peaklist=[]
        tmppeaklist=[]
        for i in indexes:
            peaklist=data[i]
            tmppeaklist=peaklist[1:]
            for j in range (0,len(tmppeaklist)):
                if abs(tmppeaklist[j][0] - tmppeaklist[j][1]) < xygap:
                    tmpms1=[tmppeaklist[j][0],tmppeaklist[j][2]]
                    ms1_data.append(tmpms1)
        siriustext=""
        siriustext=siriustext+">compound "+str(compound_name)
        siriustext=siriustext+"\n>ionization [M+?]"+str(charge)
        siriustext=siriustext+"\n>parentmass="+str(ms1_data[0][-1]).strip("[]")
        siriustext=siriustext+"\n\n>ms1"
        for i in range (0,len(ms1_data)):
            siriustext=siriustext+"\n"+str(ms1_data[i][0])+" "+str(ms1_data[i][1])
        siriustext=siriustext+"\n\n>ms2"
        for j in range (0,len(xvalues)):
            siriustext=siriustext+"\n"+str(xvalues[j])+" "+str(zvalues[j])
        return siriustext
    
    def experiment_txt_file(data,lowmass,highmass):
        exptext=""
        exptext=exptext+"Experiment\tMS Type\tMin m/z\tMax m/z"
        for i in range (0,len(data)):
            exptext=exptext+"\n"+str(i)+"\tMSMS\t"+str(lowmass)+"\t"+str(highmass)
        return exptext

        