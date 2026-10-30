# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 10:03:32 2025

@author: lwg513
"""
import numpy as np
import time
import csv
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk,StringVar,OptionMenu
from tkinter import *
from tkinter import filedialog,PhotoImage
from tkinter.font import Font
import pandas as pd
#from PIL import ImageTk, Image
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib
pythonicon = "resources\\PINK_logo_transparent.gif"
import h5py
import lmfit
from scipy.optimize import curve_fit
from scipy.signal import decimate

"""Peak INterpretation via picKing (P.I.N.K) pick for the identification of peak apexes from 2DMS surface plots
    Copyright (C) 2026  Anna Cordiner, Tony Larson, Jackie Mosely, University of York

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>."""

class TwoDpeakpickergui():
    def __init__(self):
        self.root=tk.Tk()
        self.root.title("P.I.N.K pick")
        self.root.geometry("1600x850")
        photo_image = tk.PhotoImage(file=pythonicon)
        #print(type(photo_image))
        self.root.wm_iconphoto(False, photo_image)
        
        style=ttk.Style(self.root)
        style.configure("TFrame", background="#E69", foreground="black")
        
        self.datadec=[100]
        nb = ttk.Notebook(self.root)
        self.homeframe = ttk.Frame(nb)
        nb.add(self.homeframe, text = "Home")
        
        self.inputvalues = ttk.Frame(self.homeframe, height=750, width=250)
        self.inputvalues.grid(row=1, column=0, padx=5,pady=5)
        self.plotframe = ttk.Frame(self.homeframe, height=750, width=250)
        self.plotframe.grid(row=1, column=1, padx=5,pady=5)
        self.pickvalues = ttk.Frame(self.homeframe, height=750, width=250)
        self.pickvalues.grid(row=1, column=2, padx=5,pady=5)
        
        nb.pack(expand=1, fill="both")
        self.nameprint=tk.Label(self.homeframe,text="")
        self.nameprint.grid(row=0,column=0,padx=5,pady=5)
        self.importbutton=tk.Button(self.inputvalues,text="Import surface plot",command=self.open_surface_plot)
        self.importbutton.grid(row=1, column=0, padx=50,pady=50)
        self.importoptionspre=["Touchstone","SPIKE"]
        self.importtypepre = StringVar(self.inputvalues)
        self.importtypepre.set(self.importoptionspre[0])
        self.drop=OptionMenu(self.inputvalues, self.importtypepre, *self.importoptionspre)
        self.drop.grid(row=0, column=0, padx=50,pady=50)
        self.plotcontourbutton=tk.Button(self.inputvalues,text="plot contour plot\n50x decimated",command=self.plot_contour)
        self.plotcontourbutton.grid(row=2, column=0, padx=50,pady=50)
        self.clearbutton=tk.Button(self.inputvalues,text="clear data",command=self.clearimport)
        self.clearbutton.grid(row=6,column=0,padx=50,pady=50)
        self.shortlistbutton=tk.Button(self.inputvalues,text="find local maxima",command=self.find_shortlist)
        self.shortlistbutton.grid(row=3,column=0,padx=50,pady=50)
        self.listbutton=tk.Button(self.inputvalues,text="find peak centre",command=self.find_peakcentres)
        self.listbutton.grid(row=4,column=0,padx=50,pady=50)
        
        self.con= tk.IntVar()
        self.checkcontour=tk.Checkbutton(self.inputvalues, text='Do you want a contour plot?\n50x decimated',variable=self.con, onvalue=1, offvalue=0)
        self.checkcontour.grid(row=5, column=0, padx=5,pady=5)
        
        self.highmassxlabel=tk.Label(self.pickvalues,text="high x limit for peak picking")
        self.highmassxlabel.grid(row=0, column=5, padx=5,pady=5)
        self.highmassxentry=tk.Entry(self.pickvalues)
        self.highmassxentry.insert(0, "1500") 
        self.highmassxentry.grid(row=1, column=5, padx=5,pady=5)
        self.lowmassxlabel=tk.Label(self.pickvalues,text="low x limit for peak picking")
        self.lowmassxlabel.grid(row=0, column=6, padx=5,pady=5)
        self.lowmassxentry=tk.Entry(self.pickvalues)
        self.lowmassxentry.insert(0, "100") 
        self.lowmassxentry.grid(row=1, column=6, padx=5,pady=5)
        self.highmassylabel=tk.Label(self.pickvalues,text="high y limit for peak picking")
        self.highmassylabel.grid(row=2, column=5, padx=5,pady=5)
        self.highmassyentry=tk.Entry(self.pickvalues)
        self.highmassyentry.insert(0, "1500") 
        self.highmassyentry.grid(row=3, column=5, padx=5,pady=5)
        self.lowmassylabel=tk.Label(self.pickvalues,text="low y limit for peak picking")
        self.lowmassylabel.grid(row=2, column=6, padx=5,pady=5)
        self.lowmassyentry=tk.Entry(self.pickvalues)
        self.lowmassyentry.insert(0, "100") 
        self.lowmassyentry.grid(row=3, column=6, padx=5,pady=5)
        
        self.gridlabel=tk.Label(self.pickvalues,text="set peak picking grid size")
        self.gridlabel.grid(row=6, column=5, padx=5,pady=5)
        self.gridoptions=["5x5","7x7","9x9","11x11"]
        self.gridtype = StringVar(self.pickvalues)
        self.gridtype.set(self.gridoptions[1])
        self.griddrop=OptionMenu(self.pickvalues, self.gridtype, *self.gridoptions)
        self.griddrop.grid(row=6, column=6, padx=5,pady=5)
        
        self.thresholdlabel=tk.Label(self.pickvalues, text="threhold for peak picking")
        self.thresholdlabel.grid(row=4, column=5, padx=25,pady=25, columnspan=2)
        self.thresholdentry=tk.Entry(self.pickvalues)
        self.thresholdentry.insert(0, "5E6")
        self.thresholdentry.grid(row=5, column=5,padx=25,pady=25,columnspan=2)
        
        self.ML1label=tk.Label(self.pickvalues, text="ML1\nenter before importing SPIKE data\nif calibrated through T2D")
        self.ML1label.grid(row=9, column=5, padx=5,pady=5)
        self.ML1entry=tk.Entry(self.pickvalues)
        self.ML1entry.insert(0, "0")
        self.ML1entry.grid(row=9, column=6,padx=5,pady=5)
        
        self.ML2label=tk.Label(self.pickvalues, text="ML2")
        self.ML2label.grid(row=10, column=5, padx=5,pady=5)
        self.ML2entry=tk.Entry(self.pickvalues)
        self.ML2entry.insert(0, "0")
        self.ML2entry.grid(row=10, column=6,padx=5,pady=5)
        
        self.ML3label=tk.Label(self.pickvalues, text="ML3")
        self.ML3label.grid(row=11, column=5, padx=5,pady=5)
        self.ML3entry=tk.Entry(self.pickvalues)
        self.ML3entry.insert(0, "0")
        self.ML3entry.grid(row=11, column=6,padx=5,pady=5)
        
        self.lowxinput=tk.Label(self.pickvalues,text="Low x mass is m/z 0")
        self.lowxinput.grid(row=12, column=5,padx=5,pady=5)
        
        self.lowyinput=tk.Label(self.pickvalues,text="Low y mass is m/z 0")
        self.lowyinput.grid(row=12, column=6,padx=5,pady=5)
        
        self.highinput=tk.Label(self.pickvalues,text="high mass is m/z 0")
        self.highinput.grid(row=13, column=5,padx=5,pady=5)
        
        self.sizeinput=tk.Label(self.pickvalues,text="The number of scans is 0")
        self.sizeinput.grid(row=13, column=6,padx=5,pady=5)
        
        self.savebutton=tk.Button(self.pickvalues, text="save", command=self.save_peaklist)
        self.savebutton.grid(row=7,column=5, columnspan=2, padx=5,pady=25)
        
        self.modellabel=tk.Label(self.pickvalues, text="Which model do you want to save?")
        self.modellabel.grid(row=8,column=5, padx=5,pady=25)
        self.modeloptions=["Local max","2D Gauss","3D Gauss"]
        self.modeltype = StringVar(self.pickvalues)
        self.modeltype.set(self.modeloptions[1])
        self.modeldrop=OptionMenu(self.pickvalues, self.modeltype, *self.modeloptions)
        self.modeldrop.grid(row=8, column=6, padx=5,pady=25)
        
        
        self.gausplt= tk.IntVar()
        self.gaussplots= tk.Checkbutton(self.inputvalues, text='Do want to make peak picking plots?',variable=self.gausplt, onvalue=1, offvalue=0)
        self.gaussplots.grid(row=8, column=0, padx=5,pady=5)
        
        self.contour=plt.figure(dpi=75)
        self.TwoDcanvas = FigureCanvasTkAgg(self.contour, master=self.plotframe)  
        self.TwoDcanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,columnspan=3,rowspan=6)
        self.TwoDtoolbar = NavigationToolbar2Tk(self.TwoDcanvas,self.plotframe, pack_toolbar=False)
        self.TwoDtoolbar.update()
        self.TwoDtoolbar.grid(row=7,column=1,padx=0,pady=0,columnspan=3)
        
        self.root.mainloop()
    
    def dummy(self):
        print("you have run a dummy")
        
    def open_surface_plot(self):
        self.nameprint.destroy()
        t0=time.time()
        try:
            del self.data
            del self.fragmentaxis
            del self.precursoraxis
            self.datadec=[100]
            del self.fragmentaxisdec
        except:
            print(" ")
        if self.importtypepre.get() =="Touchstone": #all code to import TS data orginated from Touchstone, Dr Callan Littlejohn
            self.dirname=tk.filedialog.askdirectory()
            self.data=np.abs(np.load(str(self.dirname+"//spectrum.npy"),allow_pickle=True))
            with open(str(self.dirname+"//axesinfo.ax"),"r") as f:
                csv_reader=csv.reader(f, delimiter=",")
                axes=[]
                for row in csv_reader:
                    axes.append(row)
                for i in range(len(axes)):
                    for j in range(len(axes[i])):
                        if axes[i][j]=="":
                            axes[i][j]=0
                        axes[i][j]=float(axes[i][j])
            self.fragmentaxis=axes[0]
            self.precursoraxis=axes[1]
            #print(self.precursoraxis[0:10],self.precursoraxis[-10:-1])
            #print(self.fragmentaxis[0:10],self.fragmentaxis[-10:-1])
            lowx=(self.fragmentaxis[-2])
            lowy=(self.precursoraxis[-2])
            high=(self.precursoraxis[1])
            self.highmass=high
            self.lowxinput=tk.Label(self.pickvalues,text="Low x mass is m/z "+str(lowx))
            self.lowxinput.grid(row=12, column=5,padx=5,pady=5)
            self.lowyinput=tk.Label(self.pickvalues,text="Low y mass is m/z "+str(lowy))
            self.lowyinput.grid(row=12, column=6,padx=5,pady=5)
            self.highinput=tk.Label(self.pickvalues,text="High mass is m/z "+str(high))
            self.highinput.grid(row=13, column=5,padx=5,pady=5)
            self.size_2D=len(self.precursoraxis)
            self.sizeinput=tk.Label(self.pickvalues,text="The number of scans is "+str(self.size_2D))
            self.sizeinput.grid(row=13, column=6,padx=5,pady=5)
            self.root.update_idletasks()
            self.nameprint=tk.Label(self.homeframe,text=self.dirname)
            self.nameprint.grid(row=0,column=0,padx=5,pady=5)
            f.close()
        
        elif self.importtypepre.get() =="SPIKE":
            self.filename=filedialog.askopenfilename(defaultextension=".msh5", filetypes=[("SPIKE Files", "*.msh5")])
            self.methodfile=filedialog.askopenfilename(defaultextension=".method", filetypes=[("Bruker method files", "*.method")])
            hdf = h5py.File(self.filename, 'r')
            i=0
            self.axis_param1=hdf['resol1/axes'][0]
            self.axis_param2=hdf['resol1/axes'][1]
            ML1=float(hdf['resol1/axes'][0][1])
            ML2=float(hdf['resol1/axes'][0][2])
            ML3=float(hdf['resol1/axes'][0][3])
            self.size_2D=int(hdf['resol1/axes'][0][11])
            TD=float(hdf['resol1/axes'][1][11])
            EXC_Freq_High=float(hdf['resol1/axes'][0][4])
            self.data=np.array(hdf['resol1/data'][i:i+self.size_2D],dtype="float32")
            #print(ML1, ML2, ML3, self.size_2D)
            #print("got to here")
            with open((str(self.methodfile)),"r") as f:
                lines=f.readlines()
                linecounter=0
                for line in lines:
                    linecounter=linecounter+1
                    if line.find('<param name="IN_26"><value>') != -1:
                        line=line.strip('<param name="IN_26">')
                        line=line.strip('</value></param>\n')
                        IN_26=float(line)
                    if line.find('<param name="EXC_Freq_Low"><value>') != -1:
                        line=line.strip('<param name="EXC_Freq_Low">')
                        line=line.strip('</value></param>\n')
                        EXC_Freq_Low=float(line)
                    if line.find('<param name="FR_low"><value>') != -1:
                        line=line.strip('<param name="FR_low">')
                        line=line.strip('</value></param>\n')
                        FR_low=float(line)
                    if line.find('<param title="Data Processing Size (SI)"  value=') != -1:
                        line=line.strip('\t<param title="Data Processing Size (SI)"  value="')
                        line=line.strip('" />\n')
                        filledatapoints=float(line)
                    if line.find('<param name="SW_h_Broadband"><value>') != -1:
                        line=line.strip('\t<param name="SW_h_Broadband"><value>')
                        line=line.strip('"</value></param>\n')
                        sw=float(line)
            nozerofills=(filledatapoints/TD)-1
            #print(EXC_Freq_High,FR_low,EXC_Freq_Low,IN_26)
            #print(self.axis_param1)
            #print(self.axis_param2)
            lowylist=in_26_mz(9.4,IN_26,EXC_Freq_Low)
            lowy=lowylist[0]
            lowx=vc_mz(9.4,EXC_Freq_High,EXC_Freq_Low)
            high=vc_mz(9.4,FR_low,EXC_Freq_Low)
            self.highmass=high
            #print(lowy)
            if float(self.ML1entry.get())==0:
                self.ML1entry.insert(0, str(round(ML1,5)))
                self.ML2entry.insert(0, str(round(ML2,5)))
                self.ML3entry.insert(0, str(round(ML3,5)))
            else:
                ML1=float(self.ML1entry.get())
                ML2=float(self.ML2entry.get())
                ML3=float(self.ML3entry.get())
            self.lowxinput=tk.Label(self.pickvalues,text="Low x mass is m/z "+str(round(lowx,5)))
            self.lowxinput.grid(row=12, column=5,padx=5,pady=5)
            self.lowyinput=tk.Label(self.pickvalues,text="Low y mass is m/z "+str(round(lowy,5)))
            self.lowyinput.grid(row=12, column=6,padx=5,pady=5)
            self.highinput=tk.Label(self.pickvalues,text="High mass is m/z "+str(round(high,5)))
            self.highinput.grid(row=13, column=5,padx=5,pady=5)
            self.sizeinput=tk.Label(self.pickvalues,text="The number of scans is "+str(self.size_2D))
            self.sizeinput.grid(row=13, column=6,padx=5,pady=5)
            self.root.update_idletasks()
            self.fragmentaxisfreq=np.linspace(0, sw, int(TD))#constructing axis here, FR_Low is unlikely to be the actual high axis limit, but not sure what is
            #self.zerofilledaxis=np.zeros(filledatapoints)
            #self.zerofilledaxis[:int(TD)]=self.fragmentaxisfreq
            #self.fragmentaxisfreq=self.zerofilledaxis
            #print(nozerofills)
            #self.fragmentaxisfreq=self.fragmentaxisfreq+self.zerofill
            #print(len(self.fragmentaxisfreq))
            self.precursoraxisfreq=np.linspace(FR_low,lowylist[1],num=int(self.size_2D))
            #print(self.fragmentaxisfreq)
            self.fragmentaxis = []
            self.fragmentaxisnotcal = []
            for i in self.fragmentaxisfreq:
                if i==0:
                    self.fragmentaxis.append(0)
                    self.fragmentaxisnotcal.append(0)
                else:
                    self.fragmentaxis.append((ML1/i)+(ML3/(i**2))+ML2)
                    self.fragmentaxisnotcal.append(vc_mz(9.4,i,EXC_Freq_Low))
            #print(self.fragmentaxis[0:25],self.fragmentaxis[-25:-1])
            #self.fragmentaxisnotcal=vc_mz(9.4,self.fragmentaxisfreq,EXC_Freq_Low) # i think i need to add calibration to the mz converter? only for the x axis though
            #self.fragmentaxis=vc_mz_cal(9.4,self.fragmentaxisfreq,EXC_Freq_Low,ML1,ML2,ML3)
            self.precursoraxis=vc_mz(9.4,self.precursoraxisfreq,EXC_Freq_Low)
            #print(self.fragmentaxiscal)
            axes=[self.precursoraxis,self.fragmentaxis]
            #plt.plot(self.fragmentaxisfreq)
            #plt.show()
            self.nameprint=tk.Label(self.homeframe,text=self.filename)
            self.nameprint.grid(row=0,column=0,padx=5,pady=5)
            f.close()
        elif self.importtypepre.get() =="T2D":
            print("not written yet, go bug anna about it if you want this feature")
            exit()

        else:
            #self.infile=filedialog.askopenfilename(defaultextension=".csv", filetypes=[("csv", "*.csv")])
            print("not written yet, go bug anna about it if you want this feature")
            exit()
        print(np.shape(self.data))
        print(len(axes[1]))
        print(len(axes[0]))
        t1=time.time()
        print(t1-t0)
        #plt.plot(axes[1])
        #plt.show()

        
    def plot_contour(self):
        t0=time.time()
        for item in self.TwoDcanvas.get_tk_widget().find_all():
            self.TwoDcanvas.get_tk_widget().delete(item)
        self.TwoDtoolbar.destroy()
        if type(self.datadec[0])==int:
            self.datadec=decimate(self.data,5,axis=1)
            self.fragmentaxisdec=decimate(self.fragmentaxis,5)
        print(np.shape(self.datadec))
        self.contour=contour_plot(self.fragmentaxisdec,self.precursoraxis,self.datadec,self.highmass)
        self.TwoDcanvas = FigureCanvasTkAgg(self.contour, master=self.plotframe)  
        self.TwoDcanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,columnspan=3,rowspan=6)
        self.TwoDtoolbar = NavigationToolbar2Tk(self.TwoDcanvas,self.plotframe, pack_toolbar=False)
        self.TwoDtoolbar.update()
        self.TwoDtoolbar.grid(row=7,column=1,padx=0,pady=0,columnspan=3)
        t1=time.time()
        print(t1-t0)
        
    def clearimport(self):
        del self.data
        del self.fragmentaxis
        del self.precursoraxis
        del self.datadec
        del self.fragmentaxisdec
        
    def find_shortlist(self):
        t0=time.time()
        for item in self.TwoDcanvas.get_tk_widget().find_all():
            self.TwoDcanvas.get_tk_widget().delete(item)
        self.TwoDtoolbar.destroy()
        #print(self.gridtype.get(),type(self.gridtype.get()))
        self.x_localmax=[]
        self.y_localmax=[]
        self.z_localmax=[]
        self.threshold=float(self.thresholdentry.get())
        self.lowxmass=float(self.lowmassxentry.get())
        self.highxmass=float(self.highmassxentry.get())
        self.lowymass=float(self.lowmassyentry.get())
        self.highymass=float(self.highmassyentry.get())
        self.masslimits=[self.lowxmass,self.highxmass,self.lowymass,self.highymass]
        self.localmaxindexes=find_local_maxima(self.data,self.fragmentaxis,self.precursoraxis,self.threshold,self.masslimits)
        #print(self.localmaxindexes)
        for i in range(0,len(self.localmaxindexes)):
            self.x_localmax.append(self.fragmentaxis[self.localmaxindexes[i,1]])
            self.y_localmax.append(self.precursoraxis[self.localmaxindexes[i,0]])
            self.z_localmax.append(self.data[self.localmaxindexes[i,0],self.localmaxindexes[i,1]])
        print(np.shape(self.x_localmax),np.shape(self.y_localmax),np.shape(self.z_localmax))
        self.shortlist=np.array([[self.x_localmax],[self.y_localmax],[self.z_localmax]])
        #print(self.shortlist)
        #self.combinedplot=scatter_and_contour_plot(self.fragmentaxis,self.precursoraxis,self.data,self.x_localmax,self.y_localmax,self.z_localmax)
        self.combinedplot=scatter_plot(self.fragmentaxis,self.precursoraxis,self.data,self.x_localmax,self.y_localmax,self.z_localmax)
        self.TwoDcanvas = FigureCanvasTkAgg(self.combinedplot, master=self.plotframe)  
        self.TwoDcanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,columnspan=3,rowspan=6)
        self.TwoDtoolbar = NavigationToolbar2Tk(self.TwoDcanvas,self.plotframe, pack_toolbar=False)
        self.TwoDtoolbar.update()
        self.TwoDtoolbar.grid(row=7,column=1,padx=0,pady=0,columnspan=3)
        t1=time.time()
        print(t1-t0)
        
    def find_peakcentres(self):
        t0=time.time()
        for item in self.TwoDcanvas.get_tk_widget().find_all():
            self.TwoDcanvas.get_tk_widget().delete(item)
        self.TwoDtoolbar.destroy()
        self.gridsize=str(self.gridtype.get())
        #print(self.gridsize)
        self.x_peakcentre=[]
        self.y_peakcentre=[]
        self.z_peakcentre=[]
        self.xcentres=[]
        self.ycentres=[]
        self.zcentres=[]
        self.slicex=[]
        self.slicey=[]
        self.slicez=[]
        self.allowedindexs=[]
        not_completed_counter=0
        completed_counter=0
        self.plotsyn=self.gausplt.get()
        #print(self.localmaxindexes[10])
        #print(self.localmaxindexes[10][0],self.localmaxindexes[10][1])
        #print(len(self.fragmentaxis))
        for i in range (0,len(self.localmaxindexes)):
            if self.gridsize == "5x5":
                if self.localmaxindexes[i][1]<2:
                    not_completed_counter+=1
                    #print(self.localmaxindexes[i][1])
                elif ((len(self.fragmentaxis))-self.localmaxindexes[i][1])<3:
                    not_completed_counter+=1
                    #print(self.localmaxindexes[i][1])
                elif self.localmaxindexes[i][0]<2:
                    y_indexes=[(self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2)]
                    x_indexes=[(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                          (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2)]
                    completed_counter+=1
                elif ((len(self.precursoraxis))-self.localmaxindexes[i][0])<3:
                     y_indexes=[(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                           (self.localmaxindexes[i][0])]
                     x_indexes=[(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                           (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2)]
                     completed_counter+=1
                else:
                    y_indexes=[(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                               (self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2)]
                    x_indexes=[(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                              (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2)]
                    completed_counter+=1
            elif self.gridsize== "7x7":
                if self.localmaxindexes[i][1]<3:
                    not_completed_counter+=1
                    #print(self.localmaxindexes[i][1])
                elif ((len(self.fragmentaxis))-self.localmaxindexes[i][1])<4:
                    not_completed_counter+=1
                    #print(self.localmaxindexes[i][1])
                elif self.localmaxindexes[i][0]<3:
                    y_indexes=[(self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2),(self.localmaxindexes[i][0]+3)]
                    x_indexes=[(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                          (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3)]
                    completed_counter+=1
                elif ((len(self.precursoraxis))-self.localmaxindexes[i][0])<4:
                     y_indexes=[(self.localmaxindexes[i][0]-3),(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                           (self.localmaxindexes[i][0])]
                     x_indexes=[(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                           (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3)]
                     completed_counter+=1
                else:
                    y_indexes=[(self.localmaxindexes[i][0]-3),(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                               (self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2),(self.localmaxindexes[i][0]+3)]
                    x_indexes=[(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                              (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3)]
                    completed_counter+=1
            elif self.gridsize == "9x9":
                if self.localmaxindexes[i][1]<4:
                    not_completed_counter+=1
                    #print(self.localmaxindexes[i][1])
                elif ((len(self.fragmentaxis))-self.localmaxindexes[i][1])<5:
                    not_completed_counter+=1
                    #print(self.localmaxindexes[i][1])
                elif self.localmaxindexes[i][0]<4:
                    y_indexes=[(self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2),(self.localmaxindexes[i][0]+3),(self.localmaxindexes[i][0]+4)]
                    x_indexes=[(self.localmaxindexes[i][1]-4),(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                          (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3),(self.localmaxindexes[i][1]+4)]
                    completed_counter+=1
                elif ((len(self.precursoraxis))-self.localmaxindexes[i][0])<5:
                     y_indexes=[(self.localmaxindexes[i][0]-4),(self.localmaxindexes[i][0]-3),(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                           (self.localmaxindexes[i][0])]
                     x_indexes=[(self.localmaxindexes[i][1]-4),(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                           (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3),(self.localmaxindexes[i][1]+4)]
                     completed_counter+=1
                else:
                    y_indexes=[(self.localmaxindexes[i][0]-4),(self.localmaxindexes[i][0]-3),(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                               (self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2),(self.localmaxindexes[i][0]+3),(self.localmaxindexes[i][0]+4)]
                    x_indexes=[(self.localmaxindexes[i][1]-4),(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                              (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3),(self.localmaxindexes[i][1]+4)]
                    completed_counter+=1
            elif self.gridsize == "11x11":
                if self.localmaxindexes[i][1]<5:
                    not_completed_counter+=1
                    print(self.localmaxindexes[i][1])
                elif ((len(self.fragmentaxis))-self.localmaxindexes[i][1])<6:
                    not_completed_counter+=1
                    print(self.localmaxindexes[i][1])
                elif self.localmaxindexes[i][0]<5:
                    y_indexes=[(self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2),(self.localmaxindexes[i][0]+3),(self.localmaxindexes[i][0]+4),(self.localmaxindexes[i][0]+5)]
                    x_indexes=[(self.localmaxindexes[i][1]-5),(self.localmaxindexes[i][1]-4),(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                          (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3),(self.localmaxindexes[i][1]+4),(self.localmaxindexes[i][1]+5)]
                    completed_counter+=1
                elif ((len(self.precursoraxis))-self.localmaxindexes[i][0])<6:
                     y_indexes=[(self.localmaxindexes[i][0]-5),(self.localmaxindexes[i][0]-4),(self.localmaxindexes[i][0]-3),(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                           (self.localmaxindexes[i][0])]
                     x_indexes=[(self.localmaxindexes[i][1]-5),(self.localmaxindexes[i][1]-4),(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                           (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3),(self.localmaxindexes[i][1]+4),(self.localmaxindexes[i][1]+5)]
                     completed_counter+=1
                else:
                    y_indexes=[(self.localmaxindexes[i][0]-5),(self.localmaxindexes[i][0]-4),(self.localmaxindexes[i][0]-3),(self.localmaxindexes[i][0]-2),(self.localmaxindexes[i][0]-1),
                               (self.localmaxindexes[i][0]),(self.localmaxindexes[i][0]+1),(self.localmaxindexes[i][0]+2),(self.localmaxindexes[i][0]+3),(self.localmaxindexes[i][0]+4),(self.localmaxindexes[i][0]+5)]
                    x_indexes=[(self.localmaxindexes[i][1]-5),(self.localmaxindexes[i][1]-4),(self.localmaxindexes[i][1]-3),(self.localmaxindexes[i][1]-2),(self.localmaxindexes[i][1]-1),
                              (self.localmaxindexes[i][1]),(self.localmaxindexes[i][1]+1),(self.localmaxindexes[i][1]+2),(self.localmaxindexes[i][1]+3),(self.localmaxindexes[i][1]+4),(self.localmaxindexes[i][1]+5)]
                    completed_counter+=1
            else:
                print ("invalid grid size 1")
                return
            if self.plotsyn==1:
                self.peak_centre=TwoD_Gauss(x_indexes,y_indexes,self.data,self.fragmentaxis,self.precursoraxis,i,self.gridtype.get(),len(self.precursoraxis))
            else:
                self.peak_centre=TwoD_Gauss_noplots(x_indexes,y_indexes,self.data,self.fragmentaxis,self.precursoraxis,i,self.gridtype.get(),len(self.precursoraxis))
            if len (self.peak_centre)!=1:
                self.xcentres.append(self.peak_centre[1])
                self.ycentres.append(self.peak_centre[2])
                self.zcentres.append(self.peak_centre[3])
                self.slicex.append(self.peak_centre[4])
                self.slicey.append(self.peak_centre[5])
                self.slicez.append(self.peak_centre[6])
            self.allowedindexs.append(self.peak_centre[0])
        #print(self.allowedindexs)
        #print(y_indexes)
        print(not_completed_counter,completed_counter)
        self.on=self.con.get()
        #print(self.on)
        if self.on==1:
            try:
                self.datadec=decimate(self.data,10,axis=1)
                self.fragmentaxisdec=decimate(self.fragmentaxis,10)
                self.combinedplot=scatter_and_contour_plot_allfit(self.fragmentaxisdec,self.precursoraxis,self.datadec,self.x_localmax,self.y_localmax,self.z_localmax,self.slicex,self.slicey,self.slicez,self.xcentres,self.ycentres,self.zcentres,self.highmass)
            except Exception as error:
                print(type(error).__name__,'-',error)
                self.combinedplot=scatter_plot_localand2fit(self.fragmentaxis,self.precursoraxis,self.data,self.x_localmax,self.y_localmax,self.z_localmax,self.slicex,self.slicey,self.slicez,self.xcentres,self.ycentres,self.zcentres)
        else:
            self.combinedplot=scatter_plot_localand2fit(self.fragmentaxis,self.precursoraxis,self.data,self.x_localmax,self.y_localmax,self.z_localmax,self.slicex,self.slicey,self.slicez,self.xcentres,self.ycentres,self.zcentres)
        self.TwoDcanvas = FigureCanvasTkAgg(self.combinedplot, master=self.plotframe)  
        self.TwoDcanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,columnspan=3,rowspan=6)
        self.TwoDtoolbar = NavigationToolbar2Tk(self.TwoDcanvas,self.plotframe, pack_toolbar=False)
        self.TwoDtoolbar.update()
        self.TwoDtoolbar.grid(row=7,column=1,padx=0,pady=0,columnspan=3)
        t1=time.time()
        print(t1-t0)
        print(len(self.slicex))
        
    def save_peaklist(self):
        #print(self.zcentres)
        self.results=[]
        if self.modeltype.get() =="Local max":
            for i in range (0,len(self.x_localmax)):
                self.results.append([self.x_localmax[i],self.y_localmax[i],self.z_localmax[i]])
        elif self.modeltype.get() =="2D Gauss":
            for i in range (0,len(self.slicex)):
                self.results.append([self.slicex[i],self.slicey[i],self.slicez[i]])
        elif self.modeltype.get() =="3D Gauss":
            for i in range (0,len(self.xcentres)):
                self.results.append([self.xcentres[i],self.ycentres[i],self.zcentres[i]])
        self.csvfile=filedialog.asksaveasfile(defaultextension=".csv",filetypes=[("CSV Files", "*.csv")])
        pd.DataFrame(self.results).to_csv(self.csvfile, header=['X','Y','Z'])
        self.csvfile.close()
        print("saved")
        
def contour_plot(x,y,array,highmass):
    indexlist=[]
    red_x=[]
    red_array=array#[:,indexlist]
    levels=[5E6,1E7,5E7,1E8,5E8,1E9,5E9]
    #red_array_log=np.log(red_array)
    fig,ax=plt.subplots(dpi=75)
    X, Y = np.meshgrid(x, y)
    cm = matplotlib.colormaps["RdPu"]
    sc =plt.contour(X, Y,red_array,levels=levels,cmap=cm)
    ax.set_ylim([min(y), max(y)+25])
    ax.set_xlim([min(x), highmass])
    ax.set_xlabel("Fragment m/z", size =20,style="italic")
    ax.set_ylabel("Precursor m/z", size=20, style="italic")
    plt.colorbar(sc, ax=ax, label="intensity", spacing='uniform')
    fig.set_figheight(10)
    fig.set_figwidth(10)
    return fig

def scatter_and_contour_plot(x,y,array,x_lm,y_lm,z_lm,highmass):
    indexlist=[]
    red_x=[]
    levels=[1E7,5E7,1E8,5E8,1E9]
    cm = matplotlib.colormaps["RdPu"]
    red_array=array#[:,indexlist]
    #red_array_log=np.log(red_array)
    fig,ax=plt.subplots(dpi=75)
    X, Y = np.meshgrid(x, y)
    #cm = matplotlib.colormaps["RdPu"]
    sc =plt.contour(X, Y,red_array,levels=levels)
    ax.scatter(x_lm, y_lm, c=np.log(z_lm), s=20, cmap=cm)
    ax.set_ylim([min(y), max(y)+25])
    ax.set_xlim([min(x), highmass])
    ax.set_xlabel("Fragment m/z", size =20,style="italic")
    ax.set_ylabel("Precursor m/z", size=20, style="italic")
    #plt.colorbar(sc, ax=ax, label="intensity", spacing='proportional')
    fig.set_figheight(10)
    fig.set_figwidth(10)
    return fig

def scatter_plot(x,y,array,x_lm,y_lm,z_lm):
    indexlist=[]
    red_x=[]
    red_array=array#[:,indexlist]
    #red_array_log=np.log(red_array)
    fig,ax=plt.subplots(dpi=75)
    cm = matplotlib.colormaps["RdPu"]
    ax.scatter(x_lm, y_lm, c=np.log(z_lm), s=20, cmap=cm)
    ax.set_ylim([min(y), max(y_lm)+25])
    ax.set_xlim([min(x), max(x_lm)+25])
    ax.set_xlabel("Fragment m/z", size =20,style="italic")
    ax.set_ylabel("Precursor m/z", size=20, style="italic")
    #plt.colorbar(sc, ax=ax, label="intensity", spacing='proportional')
    fig.set_figheight(10)
    fig.set_figwidth(10)
    return fig

def scatter_plot_localandfit(x,y,array,x_lm,y_lm,z_lm,x_f,y_f,z_f):
    #indexlist=[]
    #red_x=[]
    #red_array=array#[:,indexlist]
    #red_array_log=np.log(red_array)
    fig,ax=plt.subplots(dpi=75)
    cm = matplotlib.colormaps["RdPu"]
    ax.scatter(x_lm, y_lm, c=np.log(z_lm), s=10, cmap=cm)
    ax.scatter(x_f, y_f, c=np.log(z_f), s=10, cmap=cm, marker="x")
    ax.set_ylim([min(y), max(y_lm)+25])
    ax.set_xlim([min(x), max(x_lm)+25])
    ax.set_xlabel("Fragment m/z", size =20,style="italic")
    ax.set_ylabel("Precursor m/z", size=20, style="italic")
    #plt.colorbar(sc, ax=ax, label="intensity", spacing='proportional')
    fig.set_figheight(10)
    fig.set_figwidth(10)
    return fig

def scatter_plot_localand2fit(x,y,array,x_lm,y_lm,z_lm,x_s,y_s,z_s,x_f,y_f,z_f):
    #indexlist=[]
    #red_x=[]
    #red_array=array#[:,indexlist]
    #red_array_log=np.log(red_array)
    fig,ax=plt.subplots(dpi=75)
    cm = matplotlib.colormaps["RdPu"]
    cm2 = matplotlib.colormaps["OrRd"]
    ax.scatter(x_lm, y_lm, c=np.log(z_lm), s=10)
    ax.scatter(x_s, y_s, c=np.log(z_s), s=10, cmap=cm, marker="P")
    ax.scatter(x_f, y_f, c=np.log(z_f), s=10, cmap=cm, marker="D")
    ax.set_ylim([min(y), max(y_lm)+25])
    ax.set_xlim([min(x), max(x_lm)+25])
    ax.set_xlabel("Fragment m/z", size =20,style="italic")
    ax.set_ylabel("Precursor m/z", size=20, style="italic")
    #plt.colorbar(sc, ax=ax, label="intensity", spacing='proportional')
    fig.set_figheight(10)
    fig.set_figwidth(10)
    return fig

def scatter_and_contour_plot_allfit(x,y,array,x_lm,y_lm,z_lm,x_s,y_s,z_s,x_f,y_f,z_f,highmass):
    indexlist=[]
    red_x=[]
    cm = matplotlib.colormaps["RdPu"]
    levels=[5E6,1E7,5E7,1E8,5E8,1E9]
    red_array=array#[:,indexlist]
    #red_array_log=np.log(red_array)
    fig,ax=plt.subplots(dpi=75)
    X, Y = np.meshgrid(x, y)
    #cm = matplotlib.colormaps["RdPu"]
    sc =plt.contour(X, Y,red_array,levels=levels)
    ax.scatter(x_lm, y_lm, c=np.log(z_lm), s=10, cmap=cm)
    ax.scatter(x_s, y_s, c=np.log(z_s), s=10, cmap=cm, marker="P")
    ax.scatter(x_f, y_f, c=np.log(z_f), s=10, cmap=cm, marker="X")
    ax.set_ylim([min(y), max(y)+25])
    ax.set_xlim([min(x), highmass])
    ax.set_xlabel("Fragment m/z", size =20,style="italic")
    ax.set_ylabel("Precursor m/z", size=20, style="italic")
    #plt.colorbar(sc, ax=ax, label="intensity", spacing='proportional')
    fig.set_figheight(10)
    fig.set_figwidth(10)
    return fig

def find_local_maxima(data,fragmentaxis,precursoraxis,threshold,masslimits):
    local_maxima_index=[]
    for i in range(0,len(data)):
        if i==0:
            for j in range(0,len(data[i])):
                if data[i,j]>threshold and precursoraxis[i]>masslimits[2] and precursoraxis[i]<masslimits[3] and fragmentaxis[j]>masslimits[0] and fragmentaxis[j]<masslimits[1]:
                    if j==0:
                        if data[i,j]> data[i+1,j] and data[i,j]> data[i,j+1] and data[i,j]> data[i+2,j] and data[i,j]> data[i,j+2]:
                            local_maxima_index.append([i,j])
                    elif j==(len(data[i])-1) or j==(len(data[i])-2):
                        if data[i,j]> data[i+1,j] and data[i,j]> data[i,j-1] and data[i,j]> data[i+2,j] and data[i,j]> data[i,j-2]:
                            local_maxima_index.append([i,j]) 
                    else:
                        if data[i,j]> data[i+1,j] and data[i,j]> data[i,j+1] and data[i,j]> data[i,j-1] and data[i,j]> data[i+2,j] and data[i,j]> data[i,j+2] and data[i,j]> data[i,j-2]:
                            local_maxima_index.append([i,j]) 
        elif i==(len(data)-1) or i==(len(data)-2):
            for j in range(0,len(data[i])):
                if data[i,j]>threshold and precursoraxis[i]>masslimits[2] and precursoraxis[i]<masslimits[3] and fragmentaxis[j]>masslimits[0] and fragmentaxis[j]<masslimits[1]:
                    if j==0:
                        if data[i,j]> data[i,j+1] and data[i,j]> data[i-1,j] and data[i,j]> data[i,j+2] and data[i,j]> data[i-2,j]:
                            local_maxima_index.append([i,j])
                    elif j==(len(data[i])-1) or j==(len(data[i])-2):
                        if data[i,j]> data[i-1,j] and data[i,j]> data[i,j-1] and data[i,j]> data[i-2,j] and data[i,j]> data[i,j-2]:
                            local_maxima_index.append([i,j])
                    else:
                        if data[i,j]> data[i,j+1] and data[i,j]> data[i-1,j] and data[i,j]> data[i,j-1] and data[i,j]> data[i,j+2] and data[i,j]> data[i-2,j] and data[i,j]> data[i,j-2]:
                            local_maxima_index.append([i,j])
        else:
            for j in range(0,len(data[i])):
                if data[i,j]>threshold and precursoraxis[i]>masslimits[2] and precursoraxis[i]<masslimits[3] and fragmentaxis[j]>masslimits[0] and fragmentaxis[j]<masslimits[1]:
                    if j==0:
                        if data[i,j]> data[i+1,j] and data[i,j]> data[i,j+1] and data[i,j]> data[i-1,j] and data[i,j]> data[i+2,j] and data[i,j]> data[i,j+2] and data[i,j]> data[i-2,j]:
                            local_maxima_index.append([i,j])
                    elif j==(len(data[i])-1):
                        if data[i,j]> data[i+1,j] and data[i,j]> data[i-1,j] and data[i,j]> data[i,j-1] and data[i,j]> data[i+2,j] and data[i,j]> data[i-2,j] and data[i,j]> data[i,j-2]:
                            local_maxima_index.append([i,j])
                    else:
                        if data[i,j]> data[i+1,j] and data[i,j]> data[i,j+1] and data[i,j]> data[i-1,j] and data[i,j]> data[i,j-1] and data[i,j]> data[i+2,j] and data[i,j]> data[i,j+2] and data[i,j]> data[i-2,j] and data[i,j]> data[i,j-2]:
                            local_maxima_index.append([i,j])
    local_maxima_index_ar=np.array(local_maxima_index)
    print(len(local_maxima_index_ar),np.shape(local_maxima_index_ar),local_maxima_index_ar[0])
    return local_maxima_index_ar

def cal_sigma(B0,ind,datapoints,mz):
    sigma=(2.3548*mz**2)/(8.58e6*B0*ind*datapoints)
    return sigma

def oneD_gauss_eq(x,height,mean,sigma):
    f=height *np.exp(-((x-abs(mean))**2)/(2*sigma**2))
    return f

def Gaussian_2d(xy,A,mu_x,mu_y,sigma_x,sigma_y):
    x,y=xy
    results=np.zeros((len(x),len(y)))
    for i in range (0,len(x)):
        for j in range(0,len(y)):
            G= A*np.exp(-(x[i]-mu_x)**2/(2*sigma_x**2) - (y[j]-mu_y)**2/(2*sigma_y**2) )
            results[i,j]=G
    return np.ravel(results)
    
def TwoD_Gauss(x_indexlist,y_indexlist,data,fragmentaxis,precursoraxis,localmaxindex,gridsize,scans):
    xlist=[]
    ylist=[]
    allowed_indexs=[]
    #print(precursoraxis[512])
    for i in x_indexlist:
        #print(x_indexlist)
        xlist.append(fragmentaxis[i])
    for j in y_indexlist:
        #print(j)
        ylist.append(precursoraxis[j])
    #print(xlist)
    #print(ylist)
    int_grid=np.empty(shape=(len(ylist), len(xlist)), dtype='object')
    for i in range(0,len(x_indexlist)):
        for j in range(0,len(y_indexlist)):
            int_grid[j][i]=data[y_indexlist[j]][x_indexlist[i]]
            allowed_indexs.append([y_indexlist[j],x_indexlist[i]])
    #print(int_grid)
    #print(allowed_indexs)
    int_grid = np.array(int_grid, dtype=float)
    print(np.shape(int_grid))
    xlistext=np.arange(start=(xlist[-1]-0.005),stop=(xlist[0]+0.005),step=0.00001)
    ylistext=np.arange(start=(ylist[-1]-1),stop=(ylist[0]+1),step=0.1)
    if gridsize == "5x5":
        centreindex=2
        try:
            centerx=fragmentaxis[x_indexlist[2]]
            centery=precursoraxis[y_indexlist[2]]
            sigmaxguess=0.1*(centerx**2/1000**2)
            sigmayguess=10*(cal_sigma(21,5e-6,scans,centery))
            #print(sigmaxguess)
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[2,:],p0 = (int_grid[2,2],fragmentaxis[x_indexlist[2]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            #print(poptfx[2])
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,2],p0 = (int_grid[2,2],precursoraxis[y_indexlist[2]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            #print(poptfy)
            if poptfy[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    elif gridsize =="7x7":
        centreindex=3
        centerx=fragmentaxis[x_indexlist[3]]
        centery=precursoraxis[y_indexlist[3]]
        sigmaxguess=0.02*(centerx**2/1000**2)
        sigmayguess=10*cal_sigma(21,5e-6,scans,centery)
        #print(sigmaxguess)
        try:
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[3,:],p0 = (int_grid[3,3],fragmentaxis[x_indexlist[3]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            #print(poptfx[2])
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,3],p0 = (int_grid[3,3],precursoraxis[y_indexlist[3]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            if poptfy[2]>(sigmayguess*1.1):
               allowed_indexs=[np.NaN]
               return allowed_indexs
 
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    elif gridsize =="9x9":
        centreindex=4
        centerx=fragmentaxis[x_indexlist[4]]
        centery=precursoraxis[y_indexlist[4]]
        sigmaxguess=0.1*(centerx**2/1000**2)
        sigmayguess=10*cal_sigma(21,5e-6,scans,centery)
        #print(sigmaxguess)
        try:
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[4,:],p0 = (int_grid[4,4],fragmentaxis[x_indexlist[4]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            #print(poptfx[2])
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,4],p0 = (int_grid[4,4],precursoraxis[y_indexlist[4]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            #print(poptfy)
            if poptfy[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    elif gridsize =="11x11":
        centreindex=5
        centerx=fragmentaxis[x_indexlist[5]]
        centery=precursoraxis[y_indexlist[5]]
        sigmaxguess=0.1*(centerx**2/1000**2)
        sigmayguess=10*cal_sigma(21,5e-6,scans,centery)
        try:
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[5,:],p0 = (int_grid[5,5],fragmentaxis[x_indexlist[5]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            #print(poptfx[2])
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,5],p0 = (int_grid[5,5],precursoraxis[y_indexlist[5]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            #print(poptfy)
            if poptfy[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    else:
        print ("invalid grid size 2")
        return
    initial_guess=(poptfx[0]+1e6,poptfx[1]+0.0001,poptfy[1]+0.05,abs(2*poptfx[2]),abs(2*poptfy[2])) #initials guess contains intensity, x value, y value,y sigma, x sigma
    try:
        popt2d, pcov2d = curve_fit(Gaussian_2d, (xlist, ylist), np.ravel(int_grid), p0=initial_guess)
    except:
        allowed_indexs=[np.NaN]
        return allowed_indexs
    print(popt2d)
    
    #try:
    initial=lmfit.Parameters()
    initial.add("amplitude",value=poptfx[0]+1e6)
    initial.add("centerx",value=poptfx[1]+0.001)
    initial.add("centery",value=poptfy[1]+0.5)
    initial.add("sigmax",value=abs(2*poptfx[2]))
    initial.add("sigmay",value=abs(2*poptfy[2]))
    model = lmfit.models.Gaussian2dModel()
    result = model.fit(int_grid, x=xlist, y=ylist, params=initial)
    #lmfit.report_fit(result)
    #print(result.params.valuesdict())
    #print(result.params['centerx'].value)popt2d[1],popt2d[2],c=popt2d[0]
    MODELCENTREX=popt2d[1]
    MODELCENTREY=popt2d[2]
    MODELCENTREZ=popt2d[0]
    
    #popt3d, pcov = curve_fit(twoD_gauss_eq, (xlist, ylist), int_grid,  p0 = (poptfx[0], poptfx[1], poptfy[1], abs(poptfx[2]), abs(poptfy[2])))
    #print(popt3d)
    if result.params['sigmax'].value>0 or result.params['sigmay'].value<0:
        #fig, ax = plt.subplots(subplot_kw={"projection": "3d"},dpi=75)
        #print(ylist)
        #ax.plot_surface(xlist, ylist, int_grid)
        #ax.scatter(fragmentaxis[x_indexlist[3]],precursoraxis[y_indexlist[3]],1.2E6,color="red", s=20)
        #print(int_grid[3,3])
        #ax.scatter(data[yindexlist[3]],fragmentaxis[xindexlist[3]],data[yindexlist[3],xindexlist[3]])
        #plt.show()
        fig, ax = plt.subplots(dpi=75)
        ax.contour(xlist, ylist, int_grid)
        ax.autoscale(False)
        if gridsize =="5x5":
            ax.scatter(fragmentaxis[x_indexlist[2]],precursoraxis[y_indexlist[2]],c=int_grid[2,2], s=20)
        elif gridsize =="7x7":
            ax.scatter(fragmentaxis[x_indexlist[3]],precursoraxis[y_indexlist[3]],c=int_grid[3,3], s=20)
        elif gridsize =="9x9":
            ax.scatter(fragmentaxis[x_indexlist[4]],precursoraxis[y_indexlist[4]],c=int_grid[4,4], s=20)
        elif gridsize =="11x11":
            ax.scatter(fragmentaxis[x_indexlist[5]],precursoraxis[y_indexlist[5]],c=int_grid[5,5], s=20)
        else:
            print ("invalid grid size 3")
            return
        ax.scatter(result.params['centerx'].value,result.params['centery'].value,c=result.params['amplitude'].value, marker="x")
        ax.scatter(popt2d[1],popt2d[2],c=popt2d[0], marker="D")
        ax.scatter(poptfx[1],poptfy[1],c=poptfy[0], marker="+")
        plt.show()
        
        #print(result.params.valuesdict())
        
        #fig, ax = plt.subplots()
        #ax.scatter(xlist,int_grid[:,3])
        #ax.plot(xlist,int_grid[:,3])
        #plt.show()
        #fig, ax = plt.subplots()
        #ax.scatter(ylist,int_grid[3,:])
        #ax.plot(ylist,int_grid[3,:])
        #plt.show()
        #
        #print(precursoraxis[y_indexlist[3]],fragmentaxis[x_indexlist[3]])
        
        
        
        #popt, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[:,3],p0 = (int_grid[3,3],fragmentaxis[x_indexlist[3]],10))
        #print(popt)
        #fit_h,fit_m,fit_s=parameters
        #fit_fx = oneD_gauss_eq(xlistext, *popt)
        if gridsize =="5x5":
            plt.plot(xlist,int_grid[2,:], 'o', label='Data')
        elif gridsize =="7x7":
            plt.plot(xlist,int_grid[3,:], 'o', label='Data')
        elif gridsize =="9x9":
            plt.plot(xlist,int_grid[4,:], 'o', label='Data')
        elif gridsize =="11x11":
            plt.plot(xlist,int_grid[5,:], 'o', label='Data')
        else:
            print ("invalid grid size 4")
            return
        plt.plot(xlistext, fit_fx, '-', label=('Fit\nsigma '+str(abs(poptfx[2]))+"\ncenter "+str(abs(poptfx[1]))))
        plt.legend()
        plt.show()
        
        #print(xlistext,xlist[-1])
        #popt, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[3,:],p0 = (int_grid[3,3],precursoraxis[y_indexlist[3]],0.1))
        #print(popt)
        #fit_h,fit_m,fit_s=parameters
        #fit_fy = oneD_gauss_eq(ylistext, *popt)
        
        if gridsize =="5x5":
            plt.plot(ylist,int_grid[:,2], 'o', label='Data')
        elif gridsize =="7x7":
            plt.plot(ylist,int_grid[:,3], 'o', label='Data')
        elif gridsize =="9x9":
            plt.plot(ylist,int_grid[:,4], 'o', label='Data')
        elif gridsize =="11x11":
            plt.plot(ylist,int_grid[:,5], 'o', label='Data')
        else:
            print ("invalid grid size 5")
            return
        plt.plot(ylistext, fit_fy, '-', label=('Fit\nsigma '+str(abs(poptfy[2]))+"\ncenter "+str(abs(poptfy[1]))))
        plt.legend()
        plt.show()
        

        #except:
            #plt.plot(xlist,int_grid[:,3], 'o', label='Data')
           # plt.plot(xlistext, fit_fx, '-', label='Fit')
            #plt.legend()
           # plt.show()
            
          #  plt.plot(ylist,int_grid[3,:], 'o', label='Data')
          #  plt.plot(ylistext, fit_fy, '-', label='Fit')
          #  plt.legend()
          #  plt.show()
            
           # allowed_indexs=[np.NaN]
           # return allowed_indexs
        
        
    return allowed_indexs,MODELCENTREX,MODELCENTREY,MODELCENTREZ,abs(poptfx[1]),(poptfy[1]),(poptfx[0])

def TwoD_Gauss_noplots(x_indexlist,y_indexlist,data,fragmentaxis,precursoraxis,localmaxindex,gridsize,scans):
    xlist=[]
    ylist=[]
    allowed_indexs=[]
    #print(precursoraxis[512])
    for i in x_indexlist:
        #print(x_indexlist)
        xlist.append(fragmentaxis[i])
    for j in y_indexlist:
        #print(j)
        ylist.append(precursoraxis[j])
    #print(xlist)
    #print(ylist)
    int_grid=np.empty(shape=(len(ylist), len(xlist)), dtype='object')
    for i in range(0,len(x_indexlist)):
        for j in range(0,len(y_indexlist)):
            int_grid[j][i]=data[y_indexlist[j]][x_indexlist[i]]
            allowed_indexs.append([y_indexlist[j],x_indexlist[i]])
    #print(int_grid)
    #print(allowed_indexs)
    int_grid = np.array(int_grid, dtype=float)
    #print(np.shape(int_grid))
    xlistext=np.arange(start=(xlist[-1]-0.005),stop=(xlist[0]+0.005),step=0.00001)
    ylistext=np.arange(start=(ylist[-1]-1),stop=(ylist[0]+1),step=0.1)
    if gridsize == "5x5":
        centreindex=2
        try:
            centerx=fragmentaxis[x_indexlist[2]]
            centery=precursoraxis[y_indexlist[2]]
            sigmaxguess=0.1*(centerx**2/1000**2)
            sigmayguess=10*cal_sigma(21,5e-6,scans,centery)
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[2,:],p0 = (int_grid[2,2],fragmentaxis[x_indexlist[2]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            #print(poptfx[2])
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,2],p0 = (int_grid[2,2],precursoraxis[y_indexlist[2]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            #print(poptfy)
            if poptfy[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    elif gridsize =="7x7":
        centreindex=3
        centerx=fragmentaxis[x_indexlist[3]]
        centery=precursoraxis[y_indexlist[3]]
        sigmaxguess=0.1*(centerx**2/1000**2)
        sigmayguess=10*(cal_sigma(21,5e-6,scans,centery))
        try:
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[3,:],p0 = (int_grid[3,3],fragmentaxis[x_indexlist[3]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            #print(poptfx[2])
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,3],p0 = (int_grid[3,3],precursoraxis[y_indexlist[3]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            if poptfy[2]>(sigmayguess*1.1):
               allowed_indexs=[np.NaN]
               return allowed_indexs
 
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    elif gridsize =="9x9":
        centreindex=4
        centerx=fragmentaxis[x_indexlist[4]]
        centery=precursoraxis[y_indexlist[4]]
        sigmaxguess=0.1*(centerx**2/1000**2)
        sigmayguess=2*cal_sigma(21,5e-6,scans,centery)
        try:
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[4,:],p0 = (int_grid[4,4],fragmentaxis[x_indexlist[4]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            #print(poptfx[2])
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,4],p0 = (int_grid[4,4],precursoraxis[y_indexlist[4]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            #print(poptfy)
            if poptfy[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    elif gridsize =="11x11":
        centreindex=5
        centerx=fragmentaxis[x_indexlist[5]]
        centery=precursoraxis[y_indexlist[5]]
        sigmaxguess=0.1*(centerx**2/1000**2)
        sigmayguess=10*cal_sigma(7,1e-5,scans,centery)
        try:
            poptfx, pcov = curve_fit(oneD_gauss_eq, xlist, int_grid[5,:],p0 = (int_grid[5,5],fragmentaxis[x_indexlist[5]],sigmaxguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fx = oneD_gauss_eq(xlistext, *poptfx)
            
            #print(poptfx[2])
            if poptfx[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
            
            #print(xlistext,xlist[-1])
            poptfy, pcov = curve_fit(oneD_gauss_eq, ylist, int_grid[:,5],p0 = (int_grid[5,5],precursoraxis[y_indexlist[5]],sigmayguess))
            #print(popt)
            #fit_h,fit_m,fit_s=parameters
            fit_fy = oneD_gauss_eq(ylistext, *poptfy)
            #print(poptfy)
            if poptfy[2]>(sigmayguess*1.1):
                allowed_indexs=[np.NaN]
                return allowed_indexs
        except:
            allowed_indexs=[np.NaN]
            return allowed_indexs
    else:
        print ("invalid grid size 2")
        return
    initial_guess=(poptfx[0]+1e6,poptfx[1]+0.0001,poptfy[1]+0.05,abs(2*poptfx[2]),abs(2*poptfy[2])) #initials guess contains intensity, x value, y value,y sigma, x sigma
    try:
        popt2d, pcov2d = curve_fit(Gaussian_2d, (xlist, ylist), np.ravel(int_grid), p0=initial_guess)
    except:
        allowed_indexs=[np.NaN]
        return allowed_indexs
    #print(popt2d)
    MODELCENTREX=popt2d[1]
    MODELCENTREY=popt2d[2]
    MODELCENTREZ=popt2d[0]
    #print(popt2d)
        
    return allowed_indexs,MODELCENTREX,MODELCENTREY,MODELCENTREZ,poptfx[1],poptfy[1],poptfx[0]

def in_26_mz(B0,in_26,freq_corr):
    sf=1/in_26
    Vc=0.5*sf+freq_corr
    #print(Vc)
    mzc=(1.535611E7*B0)/Vc
    #print(mzc)
    return mzc,Vc

def vc_mz(B0,Vc,freq_corr):
    mzc=(1.535611E7*B0)/(Vc) #   f =  ml3/mz^2 + ml1/mz + ml2
    #print(mzc)
    return mzc

def vc_mz_cal(B0,Vc,freq_corr,ML1,ML2,ML3):
    mzc=((ML1/Vc)+(ML3/Vc**2)+ML2)
    #print(mzc)
    return mzc
        
if __name__=="__main__":
    r=TwoDpeakpickergui()