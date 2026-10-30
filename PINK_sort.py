# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:00:56 2024

@author: lwg513
"""

"""Peak INterpretation via picKing (P.I.N.K) sort for the extraction of MS2 peak lists from 2DMS peak lists
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

import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk,StringVar,OptionMenu
from tkinter import *
from tkinter import filedialog,PhotoImage
from scripts.peaksorter import peaksort
import pandas as pd
from PIL import ImageTk, Image
from scripts.exportmzml import exportasmzml
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
pythonicon = "resources\\PINK_logo_transparent.gif"

        
class TwoDpeaksortergui():
    def __init__(self):
        self.root=tk.Tk()
        self.root.title("P.I.N.K sort")
        self.root.geometry("1600x850")
        photo_image = tk.PhotoImage(file=pythonicon)
        #print(type(photo_image))
        self.root.wm_iconphoto(False, photo_image)
        
        #print(screen_height,screen_width)
        
        style=ttk.Style(self.root)
        style.configure("TFrame", background="#E69", foreground="black")
        
        nb = ttk.Notebook(self.root)
        self.precursorframe = ttk.Frame(nb)
        nb.add(self.precursorframe, text = "Precursors")
        self.msmsframe = ttk.Frame(nb)
        nb.add(self.msmsframe, text = "MSMS spectrum")
        self.contourframe = ttk.Frame(nb)
        nb.add(self.contourframe, text = "2DMS spectrum")
        self.neutralcontourframe = ttk.Frame(nb)
        nb.add(self.neutralcontourframe, text = "2DMS Neutral loss spectrum")
        self.neutralverticalframe = ttk.Frame(nb)
        nb.add(self.neutralverticalframe, text = "Neutral loss and common fragment spectrum")
        self.comparisonframe = ttk.Frame(nb)
        nb.add(self.comparisonframe, text = "Comparing MSMS spectra")
        self.isotopeframe = ttk.Frame(nb)
        nb.add(self.isotopeframe, text = "Combining spectra")
        self.exportframe = ttk.Frame(nb)
        nb.add(self.exportframe, text = "exporting as...")
        
        self.precursorvalues = ttk.Frame(self.precursorframe, height=750, width=250)
        self.precursorvalues.grid(row=0, column=0, padx=5,pady=5 )
        self.precursorcentre = ttk.Frame(self.precursorframe, height=750, width=100)
        self.precursorcentre.grid(row=0, column=1, padx=5,pady=5 )
        self.precursormasslist = ttk.Frame(self.precursorframe, height=750, width=250)
        self.precursormasslist.grid(row=0, column=2, padx=5,pady=5 )
        self.precursorplots=ttk.Frame(self.precursorcentre, height=500, width=100)
        self.precursorplots.grid(row=1, column=0, padx=5,pady=5 )
        self.precursorplotting=ttk.Frame(self.precursorcentre, height=100, width=100)
        self.precursorplotting.grid(row=0, column=0, padx=5,pady=5 )
        self.precursorsorting=ttk.Frame(self.precursorcentre, height=100, width=100)
        self.precursorsorting.grid(row=2, column=0, padx=5,pady=5 )
        nb.pack(expand=1, fill="both")

        self.interprecursor = ttk.Frame(self.msmsframe)
        self.interprecursor.grid(row=0, column=0, padx=5,pady=5 )
        self.msmscentre = ttk.Frame(self.msmsframe)
        self.msmscentre.grid(row=0, column=1, padx=5,pady=5)
        self.msmsplotssaving = ttk.Frame(self.msmscentre)
        self.msmsplotssaving.grid(row=0, column=0, padx=5,pady=0 )
        self.msmsplots = ttk.Frame(self.msmscentre)
        self.msmsplots.grid(row=1, column=0, padx=5,pady=0 )
        self.msmsfancy = ttk.Frame(self.msmscentre)
        self.msmsfancy.grid(row=2, column=0, padx=5,pady=0 )
        self.msmsmasslist = ttk.Frame(self.msmsframe)
        self.msmsmasslist.grid(row=0, column=2, padx=5,pady=5 )
        
        self.comparisonplots = ttk.Frame(self.comparisonframe)
        self.comparisonplots.grid(row=0, column=0, padx=5,pady=5 )
        self.comparisonplotsparam = ttk.Frame(self.comparisonframe)
        self.comparisonplotsparam.grid(row=0, column=1, padx=5,pady=0 )
        
        self.controlisotopeframe = ttk.Frame(self.isotopeframe)
        self.controlisotopeframe.grid(row=0, column=0, padx=5,pady=5 )
        self.isotopeplotframe = ttk.Frame(self.isotopeframe)
        self.isotopeplotframe.grid(row=0, column=1, padx=5,pady=5)
        self.isotopelistframe = ttk.Frame(self.isotopeframe)
        self.isotopelistframe.grid(row=0, column=2, padx=5,pady=0 )
        self.singleexportframe=ttk.Frame(self.isotopeframe)
        self.singleexportframe.grid(row=0, column=3, padx=5,pady=0 )
        
        self.NLVLselecter = ttk.Frame(self.neutralverticalframe)
        self.NLVLselecter.grid(row=0, column=0, padx=5,pady=5 )
        self.NLVLplotter = ttk.Frame(self.neutralverticalframe)
        self.NLVLplotter.grid(row=0, column=1, padx=5,pady=5 )
        self.NLVLmasslist = ttk.Frame(self.neutralverticalframe)
        self.NLVLmasslist.grid(row=0, column=2, padx=5,pady=5 )
        
        self.TS = tk.IntVar()
        self.TSbox=tk.Checkbutton(self.precursorvalues, text='Does this data need cleaning?',variable=self.TS, onvalue=1, offvalue=0)
        self.TSbox.grid(row=1, column=0, padx=5,pady=5 )
        self.highylabel=tk.Label(self.precursorvalues,text="Cut off for y axis, only needed for cleaning data")
        self.highylabel.grid(row=2, column=0, padx=5,pady=5 )
        self.highyentry=tk.Entry(self.precursorvalues)
        self.highyentry.insert(0, "5000") 
        self.highyentry.grid(row=3, column=0, padx=5,pady=5 )
        self.highxlabel=tk.Label(self.precursorvalues,text="Cut off for x axis, only needed for cleaning data")
        self.highxlabel.grid(row=4, column=0, padx=5,pady=5 )
        self.highxentry=tk.Entry(self.precursorvalues)
        self.highxentry.insert(0, "5000") 
        self.highxentry.grid(row=5, column=0, padx=5,pady=5 )
        
        self.nameprint=tk.Label(self.precursorvalues,text="")
        self.nameprint.grid(row=6,column=0,padx=5,pady=5)
        self.openfile=tk.Button(self.precursorvalues,text="Pick your .csv",command=self.open3Dpeaklist)
        self.openfile.grid(row=0, column=0, padx=5,pady=5 )
        
        self.MFlabel=tk.Label(self.precursorsorting,text="mass error/ppm")
        self.MFlabel.grid(row=0, column=0, padx=5,pady=5 )
        self.MFentry=tk.Entry(self.precursorsorting)
        self.MFentry.insert(0, "1") 
        self.MFentry.grid(row=1, column=0, padx=5,pady=5 )
        self.refmasslabel=tk.Label(self.precursorsorting,text="refernce m/z")
        self.refmasslabel.grid(row=0, column=1, padx=5,pady=5 )
        self.refmassentry=tk.Entry(self.precursorsorting)
        self.refmassentry.insert(0, "363") 
        self.refmassentry.grid(row=1, column=1, padx=5,pady=5 )
        self.refjitterlabel=tk.Label(self.precursorsorting,text="reference jitter, horizontal")
        self.refjitterlabel.grid(row=0, column=2, padx=5,pady=5 )
        self.refjitterentry=tk.Entry(self.precursorsorting)
        self.refjitterentry.insert(0, "0.012") 
        self.refjitterentry.grid(row=1, column=2, padx=5,pady=5)

        self.xyreslabel=tk.Label(self.precursorsorting,text="reference jitter, x=y")
        self.xyreslabel.grid(row=0, column=3, padx=5,pady=5 )
        self.xyresentry=tk.Entry(self.precursorsorting)
        self.xyresentry.insert(0, "0.04") 
        self.xyresentry.grid(row=1, column=3, padx=5,pady=5)

        self.yaxisroundlabel=tk.Label(self.precursorplotting,text="y axis rounding")
        self.yaxisroundlabel.grid(row=0, column=0, padx=5,pady=5 )
        self.yaxisroundentry=tk.Entry(self.precursorplotting)
        self.yaxisroundentry.insert(0, "1e8") 
        self.yaxisroundentry.grid(row=1, column=0, padx=5,pady=5 )

        self.lowmasslabel=tk.Label(self.precursorplotting,text="lowmass in x")
        self.lowmasslabel.grid(row=0, column=1, padx=5,pady=5 )
        self.lowmassentry=tk.Entry(self.precursorplotting)
        self.lowmassentry.insert(0, "100") 
        self.lowmassentry.grid(row=1, column=1, padx=5,pady=5 )

        self.highmasslabel=tk.Label(self.precursorplotting,text="highmass")
        self.highmasslabel.grid(row=0, column=2, padx=5,pady=5 )
        self.highmassentry=tk.Entry(self.precursorplotting)
        self.highmassentry.insert(0, "1500") 
        self.highmassentry.grid(row=1, column=2, padx=5,pady=5 )

        self.peaks=[[],[],[]]
        self.prepeaks=[[],[]]
        self.preplot=plt.figure(figsize=(12,8),dpi=75)
        self.minipreplot=plt.figure(figsize=(3,2),dpi=75)
        self.msmsplot=plt.figure(figsize=(12,8),dpi=75)
        self.isotopeplot=plt.figure(figsize=(12,8),dpi=75)
        self.NLVLplot=plt.figure(figsize=(12,8),dpi=75)
        self.contourplot=plt.figure(dpi=75)
        self.neutralcontourplot=plt.figure(dpi=75)
        self.assignnonassign=[[],[]]

        self.startsortbutton=tk.Button(self.precursorsorting,text="sort",command=self.startsort,height= 3, width=10)
        self.startsortbutton.grid(row=0, column=6, padx=5,pady=5,rowspan=2 )
        self.indictorcanvas = tk.Canvas(self.precursorsorting,bg="#E69",height=50,width=50,highlightthickness=0)
        self.indictorcanvas.grid(row=0, column=7,padx=5,pady=5,rowspan=2)
        self.indicator = self.indictorcanvas.create_oval(20, 10, 10, 20)
        self.indictorcanvas.itemconfig(self.indicator, fill="green")
        
        self.replotbutton=tk.Button(self.precursorplotting,text="replot with new parameters",command=self.new_preplot_param)
        self.replotbutton.grid(row=0, column=3, padx=5,pady=5 )
        
        self.savepreplotbutton=tk.Button(self.precursorplotting,text="save your figure",command=self.save_prefigure)
        self.savepreplotbutton.grid(row=0, column=5, padx=5,pady=5 )
        self.saveoptionspre=["SVG","PDF","JPG","PNG"]
        self.savetypepre = StringVar(self.precursorplotting)
        self.savetypepre.set(self.saveoptionspre[0])
        self.drop=OptionMenu(self.precursorplotting, self.savetypepre, *self.saveoptionspre)
        self.drop.grid(row=0, column=4, padx=5,pady=5)
        self.saveprelistbutton=tk.Button(self.precursorplotting,text="save your precursor peaklist",command=self.save_prelist)
        self.saveprelistbutton.grid(row=0, column=6, padx=5,pady=5 )
        
        self.precanvas = FigureCanvasTkAgg(self.preplot, master=self.precursorplots)  
        self.precanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0)
        
        self.miniprecanvas = FigureCanvasTkAgg(self.minipreplot, master=self.interprecursor)  
        self.miniprecanvas.get_tk_widget().grid(row=4,column=0,padx=0,pady=0)
        
        self.msmscanvas = FigureCanvasTkAgg(self.msmsplot, master=self.msmsplots)  
        self.msmscanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0)
        
        self.contourcanvas = FigureCanvasTkAgg(self.contourplot, master=self.contourframe)
        self.contourcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0)
        self.contourtoolbar = NavigationToolbar2Tk(self.contourcanvas,self.contourframe, pack_toolbar=False)
        self.contourtoolbar.update()
        self.contourtoolbar.grid(row=0,column=1,padx=0,pady=0)
        
        self.neutralcontourcanvas = FigureCanvasTkAgg(self.neutralcontourplot, master=self.neutralcontourframe)
        self.neutralcontourcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0)
        self.neutralcontourtoolbar = NavigationToolbar2Tk(self.neutralcontourcanvas,self.neutralcontourframe, pack_toolbar=False)
        self.neutralcontourtoolbar.update()
        self.neutralcontourtoolbar.grid(row=0,column=1,padx=0,pady=0)
        
        self.isotopecanvas = FigureCanvasTkAgg(self.isotopeplot, master=self.isotopeplotframe)  
        self.isotopecanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0,columnspan=3)
        
        self.NLVLcanvas = FigureCanvasTkAgg(self.NLVLplot, master=self.NLVLplotter)  
        self.NLVLcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0,columnspan=3)
        
        self.NVbg=plt.figure(figsize=(3,2))
        self.NLcanvas = FigureCanvasTkAgg(self.NVbg,master=self.NLVLselecter)  
        self.NLcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=10,rowspan=3)
        
        self.VLbg=plt.figure(figsize=(3,2))
        self.VLcanvas = FigureCanvasTkAgg(self.VLbg,master=self.NLVLselecter)  
        self.VLcanvas.get_tk_widget().grid(row=3,column=0,padx=0,pady=10,rowspan=3)
        
        self.NLlabel1=tk.Label(self.NLVLselecter,text="Neutral loss\nFind all peaks:")
        self.NLlabel1.grid(row=0,column=0,padx=5,pady=5)
        self.NLEntry=tk.Entry(self.NLVLselecter)
        self.NLEntry.insert(0, "0")
        self.NLEntry.grid(row=1,column=0,padx=5,pady=5)
        self.NLlabel2=tk.Label(self.NLVLselecter,text="away from the precursor")
        self.NLlabel2.grid(row=2,column=0,padx=5,pady=5)
        
        self.VLlabel1=tk.Label(self.NLVLselecter,text="Common Fragments\nFind all precursors that have:")
        self.VLlabel1.grid(row=3,column=0,padx=5,pady=5)
        self.VLEntry=tk.Entry(self.NLVLselecter)
        self.VLEntry.insert(0, "0")
        self.VLEntry.grid(row=4,column=0,padx=5,pady=5)
        self.VLlabel2=tk.Label(self.NLVLselecter,text="as a fragment")
        self.VLlabel2.grid(row=5,column=0,padx=5,pady=5)
        
        self.NLVLoptions=["NL","VL"]
        self.NLVLtype = StringVar(self.NLVLselecter)
        self.NLVLtype.set(self.NLVLoptions[0])
        self.dropNLVL=OptionMenu(self.NLVLselecter, self.NLVLtype, *self.NLVLoptions)
        self.dropNLVL.grid(row=6, column=0, padx=5,pady=5)
        
        self.plotNLVLbutton=tk.Button(self.NLVLselecter,text="PLOT",command=self.NLVLsortandplot)
        self.plotNLVLbutton.grid(row=7, column=0, padx=5,pady=5)
        
        self.yaxisroundNLVLlabel=tk.Label(self.NLVLplotter,text="y axis rounding")
        self.yaxisroundNLVLlabel.grid(row=1, column=0, padx=5,pady=5 )
        self.yaxisroundNLVLentry=tk.Entry(self.NLVLplotter)
        self.yaxisroundNLVLentry.insert(0, "1e8") 
        self.yaxisroundNLVLentry.grid(row=2, column=0, padx=5,pady=5 )

        self.lowmassNLVLlabel=tk.Label(self.NLVLplotter,text="lowmass in x")
        self.lowmassNLVLlabel.grid(row=1, column=1, padx=5,pady=5 )
        self.lowmassNLVLentry=tk.Entry(self.NLVLplotter)
        self.lowmassNLVLentry.insert(0, "100") 
        self.lowmassNLVLentry.grid(row=2, column=1, padx=5,pady=5 )

        self.highmassNLVLlabel=tk.Label(self.NLVLplotter,text="highmass")
        self.highmassNLVLlabel.grid(row=1, column=2, padx=5,pady=5 )
        self.highmassNLVLentry=tk.Entry(self.NLVLplotter)
        self.highmassNLVLentry.insert(0, "1500") 
        self.highmassNLVLentry.grid(row=2, column=2, padx=5,pady=5 )
        
        self.saveNLVLplotbutton=tk.Button(self.NLVLselecter,text="save your figure",command=self.save_NLVLfigure)
        self.saveNLVLplotbutton.grid(row=8, column=0, padx=5,pady=5 )
        self.saveoptionsNLVL=["SVG","PDF","JPG","PNG"]
        self.savetypeNLVL = StringVar(self.NLVLselecter)
        self.savetypeNLVL.set(self.saveoptionsNLVL[0])
        self.dropNLVL=OptionMenu(self.NLVLselecter, self.savetypeNLVL, *self.saveoptionsNLVL)
        self.dropNLVL.grid(row=9, column=0, padx=5,pady=5)
        self.saveNLVLlistbutton=tk.Button(self.NLVLselecter,text="save your peaklist",command=self.save_NLVLlist)
        self.saveNLVLlistbutton.grid(row=10, column=0, padx=5,pady=5)
        
        self.lowmassiso=tk.Entry(self.isotopeplotframe)
        self.lowmassiso.insert(0, "100") 
        self.lowmassiso.grid(row=1, column=0, padx=5, pady=5)
        self.highmassiso=tk.Entry(self.isotopeplotframe)
        self.highmassiso.insert(0, "1500") 
        self.highmassiso.grid(row=1, column=1, padx=5, pady=5)
        self.roundtoiso=tk.Entry(self.isotopeplotframe)
        self.roundtoiso.insert(0, "1e8") 
        self.roundtoiso.grid(row=1, column=2, padx=5, pady=5)
        self.lowmassisolabel=tk.Label(self.isotopeplotframe,text="low mass")
        self.lowmassisolabel.grid(row=2, column=0, padx=5, pady=5)
        self.highmassisolabel=tk.Label(self.isotopeplotframe,text="high mass")
        self.highmassisolabel.grid(row=2, column=1, padx=5, pady=5)
        self.roundupisolabel=tk.Label(self.isotopeplotframe,text="round to")
        self.roundupisolabel.grid(row=2, column=2, padx=5, pady=5)
        
        self.savemsmsplotbutton=tk.Button(self.msmsplotssaving,text="save your figure",command=self.save_msmsfigure)
        self.savemsmsplotbutton.grid(row=0, column=1, padx=5,pady=5 )
        self.saveoptionsmsms=["SVG","PDF","JPG","PNG"]
        self.savetypemsms = StringVar(self.msmsplotssaving)
        self.savetypemsms.set(self.saveoptionsmsms[0])
        self.drop=OptionMenu(self.msmsplotssaving, self.savetypemsms, *self.saveoptionsmsms)
        self.drop.grid(row=0, column=0, padx=5,pady=5)
        self.savemsmslistbutton=tk.Button(self.msmsplotssaving,text="save your tandem peaklist",command=self.save_msmslist)
        self.savemsmslistbutton.grid(row=0, column=2, padx=5,pady=5 )
        self.addtocomparebutton=tk.Button(self.msmsfancy,text="add to comparison",command=self.add_to_compare)
        self.addtocomparebutton.grid(row=0, column=0, padx=5,pady=5 )
        self.highlightbutton=tk.Button(self.msmsfancy,text="highlight potential ridges",command=self.highlight_ridges)
        self.highlightbutton.grid(row=0, column=1, padx=5,pady=5 )
        self.removebutton=tk.Button(self.msmsfancy,text="remove potential ridges",command=self.remove_ridges)
        self.removebutton.grid(row=0, column=2, padx=5,pady=5 )
        self.lowmassmsms=tk.Entry(self.msmsfancy)
        self.lowmassmsms.insert(0, "100") 
        self.lowmassmsms.grid(row=1, column=0, padx=5, pady=5)
        self.highmassmsms=tk.Entry(self.msmsfancy)
        self.highmassmsms.insert(0, "1500") 
        self.highmassmsms.grid(row=1, column=1, padx=5, pady=5)
        self.roundtomsms=tk.Entry(self.msmsfancy)
        self.roundtomsms.insert(0, "1e8") 
        self.roundtomsms.grid(row=1, column=2, padx=5, pady=5)
        self.lowmassmsmslabel=tk.Label(self.msmsfancy,text="low mass")
        self.lowmassmsmslabel.grid(row=2, column=0, padx=5, pady=5)
        self.highmassmsmslabel=tk.Label(self.msmsfancy,text="high mass")
        self.highmassmsmslabel.grid(row=2, column=1, padx=5, pady=5)
        self.roundupmsmslabel=tk.Label(self.msmsfancy,text="round to")
        self.roundupmsmslabel.grid(row=2, column=2, padx=5, pady=5)
        
        self.twovalueslabel=tk.Label(self.comparisonplots,text="You need at least 2 selected precursors to plot them here")
        self.twovalueslabel.grid(row=0, column=0, padx=5,pady=5 )
        self.compareclearbutton=tk.Button(self.comparisonplotsparam,text="clear your comparison",command=self.clearcomparison)
        self.compareclearbutton.grid(row=0, column=0, padx=5, pady=5)
        self.highmasscomparelabel=tk.Label(self.comparisonplotsparam,text="high mass")
        self.highmasscomparelabel.grid(row=3, column=0, padx=5,pady=5 )
        self.highmasscompareentry=tk.Entry(self.comparisonplotsparam)
        self.highmasscompareentry.insert(0, "1500") 
        self.highmasscompareentry.grid(row=4, column=0, padx=5,pady=5 )
        self.lowmasscomparelabel=tk.Label(self.comparisonplotsparam,text="low mass")
        self.lowmasscomparelabel.grid(row=1, column=0, padx=5,pady=5 )
        self.lowmasscompareentry=tk.Entry(self.comparisonplotsparam)
        self.lowmasscompareentry.insert(0, "100") 
        self.lowmasscompareentry.grid(row=2, column=0, padx=5,pady=5)
        self.ymaxcomparelabel=tk.Label(self.comparisonplotsparam,text="common y max value")
        self.ymaxcomparelabel.grid(row=5, column=0, padx=5,pady=5 )
        self.ymaxcompareentry=tk.Entry(self.comparisonplotsparam)
        self.ymaxcompareentry.insert(0, "1e10") 
        self.ymaxcompareentry.grid(row=6, column=0, padx=5,pady=5)
        self.replotcomparebutton=tk.Button(self.comparisonplotsparam,text="replot with new parameters",command=self.replot_compare)
        self.replotcomparebutton.grid(row=7, column=0, padx=5, pady=5)
        self.savecomparebutton=tk.Button(self.comparisonplotsparam,text="save comparison figure",command=self.save_comparefigure)
        self.savecomparebutton.grid(row=8, column=0, padx=5, pady=5)
        self.saveoptionscompare=["SVG","PDF","JPG","PNG"]
        self.savetypecompare = StringVar(self.comparisonplotsparam)
        self.savetypecompare.set(self.saveoptionscompare[0])
        self.drop=OptionMenu(self.comparisonplotsparam, self.savetypecompare, *self.saveoptionscompare)
        self.drop.grid(row=9, column=0, padx=5,pady=5)
        
        self.plotisotopesbutton=tk.Button(self.controlisotopeframe,text="plot the selected isotopolugues",command=self.plot_isotopes)
        self.plotisotopesbutton.grid(row=2, column=0, padx=5, pady=5)
        self.clearisotopesbutton=tk.Button(self.controlisotopeframe,text="clear the list of isotopolugues",command=self.clearisotope)
        self.clearisotopesbutton.grid(row=3, column=0, padx=5, pady=5)
        self.saveisofigurebutton=tk.Button(self.controlisotopeframe,text="save the combined isotope figure",command=self.save_isotopefigure)
        self.saveisofigurebutton.grid(row=4, column=0, padx=5, pady=5)
        self.saveisocsvbutton=tk.Button(self.controlisotopeframe,text="save the combined isotope peak list",command=self.save_isotopepeaklist)
        self.saveisocsvbutton.grid(row=5, column=0, padx=5, pady=5)
        self.saveoptionsiso=["SVG","PDF","JPG","PNG"]
        self.savetypeiso = StringVar(self.controlisotopeframe)
        self.savetypeiso.set(self.saveoptionsiso[0])
        self.drop=OptionMenu(self.controlisotopeframe, self.savetypeiso, *self.saveoptionsiso)
        self.drop.grid(row=6, column=0, padx=5,pady=5)
        
        self.exportoptionssingle=["Sirius","mzML"]
        self.singleexporttype = StringVar(self.singleexportframe)
        self.singleexporttype.set(self.exportoptionssingle[0])
        self.drop=OptionMenu(self.singleexportframe, self.singleexporttype, *self.exportoptionssingle)
        self.drop.grid(row=1, column=0, padx=5,pady=5)
        self.compoundlabel=tk.Label(self.singleexportframe,text="Compound Name\n or Name for file")
        self.compoundlabel.grid(row=2, column=0, padx=5,pady=5)
        self.compound=tk.Entry(self.singleexportframe)
        self.compound.grid(row=3, column=0, padx=5,pady=5) 
        self.chargelabel=tk.Label(self.singleexportframe,text="Charge of ion")
        self.chargelabel.grid(row=4, column=0, padx=5,pady=5)
        self.chargeinput=tk.Entry(self.singleexportframe)
        self.chargeinput.insert(0, "1+")
        self.chargeinput.grid(row=5, column=0, padx=5,pady=5)
        self.fraglabel=tk.Label(self.singleexportframe,text="Type of fragmentation")
        self.fraglabel.grid(row=6, column=0, padx=5,pady=5)
        self.frag_input=tk.Entry(self.singleexportframe)
        self.frag_input.insert(0, "IRMPD")
        self.frag_input.grid(row=7, column=0, padx=5,pady=5)
        self.exportsinglebutton=tk.Button(self.singleexportframe,text="Export",command=self.single_export)
        self.exportsinglebutton.grid(row=8, column=0, padx=5,pady=5)
        self.mzmlnote=tk.Label(self.singleexportframe,text="If exporting as an mzML \nplease enter full parameters on the \nexporting as... tab")
        self.mzmlnote.grid(row=9, column=0, padx=5, pady=5)
        
        self.exportbutton=tk.Button(self.exportframe,text="export as...",command=self.export_as)
        self.exportbutton.grid(row=2, column=0, padx=5,pady=5 )
        self.exportoptions=["MASCOT","mzML"]
        self.exporttype = StringVar(self.exportframe)
        self.exporttype.set(self.exportoptions[0])
        self.drop=OptionMenu(self.exportframe, self.exporttype, *self.exportoptions)
        self.drop.grid(row=0, column=0, padx=5,pady=5)
        self.exp_text=tk.Button(self.exportframe,text="exp file for LipidDial",command=self.exp_file_LipidDial)
        self.exp_text.grid(row=1, column=0, padx=5,pady=5)
        
        self.export_analyser_type_label=tk.Label(self.exportframe,text="type of analyser")
        self.export_source_type_label=tk.Label(self.exportframe,text="type of source")
        self.export_Polarity_label=tk.Label(self.exportframe,text="positive or negative mode")
        self.export_fragmentation_technique_label=tk.Label(self.exportframe,text="fragmentation technique")
        self.export_fragmentation_power_label=tk.Label(self.exportframe,text="laser power %")
        self.export_Fragmentation_duration_label=tk.Label(self.exportframe,text="fragmenttaion duration/s")
        self.export_MS_level_label=tk.Label(self.exportframe,text="MS2 or MS3")
        self.export_ExD_bias_label=tk.Label(self.exportframe,text="ExD bias voltage")
        self.export_ExD_lens_label=tk.Label(self.exportframe,text="ExD lens voltage")
        self.sourceoptions=["ESI","nESI","APCI","MALDI"]
        self.sourcetype = StringVar(self.exportframe)
        self.sourcetype.set(self.sourceoptions[0])
        self.export_source_type=OptionMenu(self.exportframe, self.sourcetype, *self.sourceoptions)
        self.analyseroptions=["FT_ICR"]
        self.analysertype = StringVar(self.exportframe)
        self.analysertype.set(self.analyseroptions[0])
        self.export_analyser_type=OptionMenu(self.exportframe, self.analysertype, *self.analyseroptions)
        self.polarityoptions=["positive scan","negative scan"]
        self.polaritytype = StringVar(self.exportframe)
        self.polaritytype.set(self.polarityoptions[0])
        self.export_Polarity=OptionMenu(self.exportframe, self.polaritytype, *self.polarityoptions)
        self.fragmentation_techniqueoptions=["IRMPD","UVPD","ExD"]
        self.fragmentation_techniquetype = StringVar(self.exportframe)
        self.fragmentation_techniquetype.set(self.fragmentation_techniqueoptions[0])
        self.export_fragmentation_technique=OptionMenu(self.exportframe, self.fragmentation_techniquetype, *self.fragmentation_techniqueoptions)
        self.export_fragmentation_power=tk.Entry(self.exportframe)
        self.export_fragmentation_power.insert(0, "50") 
        self.export_Fragmentation_duration=tk.Entry(self.exportframe)
        self.export_Fragmentation_duration.insert(0, "0.1")
        self.MS_leveloptions=["2","3"]
        self.MS_leveltype = StringVar(self.exportframe)
        self.MS_leveltype.set(self.MS_leveloptions[0])
        self.export_MS_level=OptionMenu(self.exportframe, self.MS_leveltype, *self.MS_leveloptions)
        self.export_ExD_bias=tk.Entry(self.exportframe)
        self.export_ExD_bias.insert(0, "10")
        self.export_ExD_lens=tk.Entry(self.exportframe)
        self.export_ExD_lens.insert(0, "10")
        self.export_analyser_type.grid(row=1, column=1, padx=5,pady=5)
        self.export_source_type.grid(row=3, column=1, padx=5,pady=5)
        self.export_Polarity.grid(row=5, column=1, padx=5,pady=5)
        self.export_fragmentation_technique.grid(row=7, column=1, padx=5,pady=5)
        self.export_fragmentation_power.grid(row=9, column=1, padx=5,pady=5)
        self.export_Fragmentation_duration.grid(row=15, column=1, padx=5,pady=5)
        self.export_MS_level.grid(row=17, column=1, padx=5,pady=5)
        self.export_ExD_bias.grid(row=11, column=1, padx=5,pady=5)
        self.export_ExD_lens.grid(row=13, column=1, padx=5,pady=5)
        self.export_analyser_type_label.grid(row=0, column=1, padx=5,pady=5)
        self.export_source_type_label.grid(row=2, column=1, padx=5,pady=5)
        self.export_Polarity_label.grid(row=4, column=1, padx=5,pady=5)
        self.export_fragmentation_technique_label.grid(row=6, column=1, padx=5,pady=5)
        self.export_fragmentation_power_label.grid(row=8, column=1, padx=5,pady=5)
        self.export_Fragmentation_duration_label.grid(row=14, column=1, padx=5,pady=5)
        self.export_MS_level_label.grid(row=16, column=1, padx=5,pady=5)
        self.export_ExD_bias_label.grid(row=10, column=1, padx=5,pady=5)
        self.export_ExD_lens_label.grid(row=12, column=1, padx=5,pady=5)
        
        self.premasslistlabel=tk.Label(self.precursormasslist,text="precusor peaklist")
        self.premasslistlabel.grid(row=0, column=0,padx=5,pady=5)
        self.premasslist=ttk.Treeview(self.precursormasslist,height=30)
        vsb = ttk.Scrollbar(self.precursormasslist, orient="vertical", command=self.premasslist.yview)
        vsb.grid(row=1, column =0, sticky=N+S+W)
        self.premasslist.configure(yscrollcommand=vsb.set)
        self.premasslist["columns"]=("m/z","I")
        self.premasslist.column("#0",width=0)
        self.premasslist.column("m/z",width=90)
        self.premasslist.column("I",width=100)
        self.premasslist.heading("m/z",text="Precursor m/z")
        self.premasslist.heading("I",text="Intensity")
        self.premasslist.grid(row=1, column=0, padx=5,pady=5)

        self.selectpremasslabel=tk.Label(self.interprecursor,text="list of precusrors\ndoubleclick to select")
        self.selectpremasslabel.grid(row=0, column=0,padx=5,pady=5)
        self.sortbymzbutton=tk.Button(self.interprecursor,text="sort by m/z",command=self.sort_by_mz)
        self.sortbymzbutton.grid(row=1, column=0,padx=5,pady=5)
        self.sortbyIbutton=tk.Button(self.interprecursor,text="sort by intensity",command=self.sort_by_intensity)
        self.sortbyIbutton.grid(row=2, column=0,padx=5,pady=5)
        self.selectpremass=ttk.Treeview(self.interprecursor,height=20)
        vsb2 = ttk.Scrollbar(self.interprecursor, orient="vertical", command=self.selectpremass.yview)
        vsb2.grid(row=3, column =0, sticky=N+S+E)
        self.selectpremass.configure(yscrollcommand=vsb2.set)
        self.selectpremass["columns"]=("m/z","I")
        self.selectpremass.column("#0",width=0)
        self.selectpremass.column("m/z",width=100)
        self.selectpremass.column("I",width=100)
        self.selectpremass.heading("m/z",text="Precursor m/z")
        self.selectpremass.heading("I",text="Intensity")
        self.selectpremass.grid(row=3, column=0, padx=5,pady=5) 
        self.selectpremass.bind("<Double-1>", self.selected_precusor)
        
        self.selectisotopelabel=tk.Label(self.controlisotopeframe,text="list of precusrors\ndoubleclick to select isotopolugues")
        self.selectisotopelabel.grid(row=0, column=0,padx=5,pady=5)
        self.selectisotope=ttk.Treeview(self.controlisotopeframe,height=20)
        vsb4 = ttk.Scrollbar(self.controlisotopeframe, orient="vertical", command=self.selectisotope.yview)
        vsb4.grid(row=1, column =0, sticky=N+S+E)
        self.selectisotope.configure(yscrollcommand=vsb4.set)
        self.selectisotope["columns"]=("m/z","I")
        self.selectisotope.column("#0",width=0)
        self.selectisotope.column("m/z",width=100)
        self.selectisotope.column("I",width=100)
        self.selectisotope.heading("m/z",text="Precursor m/z")
        self.selectisotope.heading("I",text="Intensity")
        self.selectisotope.grid(row=1, column=0, padx=5,pady=5) 
        self.selectisotope.bind("<Double-1>", self.selected_isotope)
        
        self.isotopelabel=tk.Label(self.isotopelistframe,text="list of isotopes\ndoubleclick to delete isotopolugues")
        self.isotopelabel.grid(row=0, column=0,padx=5,pady=5)
        self.isotopepeaks=ttk.Treeview(self.isotopelistframe,height=10)
        vsb4 = ttk.Scrollbar(self.isotopelistframe, orient="vertical", command=self.isotopepeaks.yview)
        vsb4.grid(row=1, column =0, sticky=N+S+E)
        self.isotopepeaks.configure(yscrollcommand=vsb4.set)
        self.isotopepeaks["columns"]=("m/z")
        self.isotopepeaks.column("#0",width=0)
        self.isotopepeaks.column("m/z",width=100)
        self.isotopepeaks.heading("m/z",text="peak m/z")
        self.isotopepeaks.grid(row=1, column=0, padx=5,pady=5) 
        self.isotopepeaks.bind("<Double-1>", self.delete_isotope)
        
        self.isotopepeaklist=ttk.Treeview(self.isotopelistframe,height=15)
        vsb4 = ttk.Scrollbar(self.isotopelistframe, orient="vertical", command=self.isotopepeaks.yview)
        vsb4.grid(row=2, column =0, sticky=N+S+E)
        self.isotopepeaklist.configure(yscrollcommand=vsb4.set)
        self.isotopepeaklist["columns"]=("m/z","I","P")
        self.isotopepeaklist.column("#0",width=0)
        self.isotopepeaklist.column("m/z",width=100)
        self.isotopepeaklist.column("I",width=100)
        self.isotopepeaklist.column("P",width=100)
        self.isotopepeaklist.heading("m/z",text="peak m/z")
        self.isotopepeaklist.heading("I",text="peak intensity")
        self.isotopepeaklist.heading("P",text="Precursor")
        self.isotopepeaklist.grid(row=2, column=0, padx=5,pady=5) 
        
        self.msmspeakslabel=tk.Label(self.msmsmasslist,text="tandem peaklist")
        self.msmspeakslabel.grid(row=0, column=0,padx=5,pady=5)
        self.msmspeaks=ttk.Treeview(self.msmsmasslist,height=30)
        vsb3 = ttk.Scrollbar(self.msmsmasslist, orient="vertical", command=self.msmspeaks.yview)
        vsb3.grid(row=1, column =1, sticky=N+S+E)
        self.msmspeaks.configure(yscrollcommand=vsb3.set)
        self.msmspeaks["columns"]=("m/z","I","pre")
        self.msmspeaks.column("#0",width=0)
        self.msmspeaks.column("m/z",width=80)
        self.msmspeaks.column("I",width=70)
        self.msmspeaks.column("pre",width=100)
        self.msmspeaks.heading("m/z",text="fragment m/z")
        self.msmspeaks.heading("I",text="intensity")
        self.msmspeaks.heading("pre",text="precursor m/z")
        self.msmspeaks.grid(row=1, column=0, padx=5,pady=5) 
        
        self.NLVLpeakslabel=tk.Label(self.NLVLmasslist,text="Neutral loss peaklist")
        self.NLVLpeakslabel.grid(row=0, column=0,padx=5,pady=5)
        self.NLVLpeaks=ttk.Treeview(self.NLVLmasslist,height=30)
        vsb3 = ttk.Scrollbar(self.NLVLmasslist, orient="vertical", command=self.NLVLpeaks.yview)
        vsb3.grid(row=1, column =1, sticky=N+S+E)
        self.NLVLpeaks.configure(yscrollcommand=vsb3.set)
        self.NLVLpeaks["columns"]=("m/z","pre","I")
        self.NLVLpeaks.column("#0",width=0)
        self.NLVLpeaks.column("m/z",width=80)
        self.NLVLpeaks.column("I",width=70)
        self.NLVLpeaks.column("pre",width=100)
        self.NLVLpeaks.heading("m/z",text="fragment m/z")
        self.NLVLpeaks.heading("I",text="intensity")
        self.NLVLpeaks.heading("pre",text="precursor m/z")
        self.NLVLpeaks.grid(row=1, column=0, padx=5,pady=5)
        
        self.root.mainloop()
        
    def open3Dpeaklist(self):
        self.infile=filedialog.askopenfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])
        self.highymass=float(self.highyentry.get())
        self.highxmass=float(self.highxentry.get())
        self.indictorcanvas.itemconfig(self.indicator, fill="red")
        with open(str(self.infile),"r") as f:
            self.peaks=[[],[],[]]
            reader=pd.read_csv(f)
            reader.head()
            self.fragmentmass=reader["X"]
            self.precursormass=reader["Y"]
            self.intensity=reader["Z"]
        f.close()
        print("imported "+str(len(self.fragmentmass))+" peaks")
        self.TSon=self.TS.get()
        if self.TSon==1:
            tmpindexvalues=[]
            print("Any x values above "+str(self.highxmass)+" or below 0 will be removed and any y values above "+str(self.highymass)+" or below 0 will be removed.")
            for i in range(0,len(self.fragmentmass)):
                if self.fragmentmass[i] >=self.highxmass or self.fragmentmass[i]<=0 or self.precursormass[i] >=self.highymass or self.precursormass[i]<=0:
                    tmpindexvalues.append(i)
            self.fragmentmass=self.fragmentmass.drop(index=tmpindexvalues)
            self.precursormass=self.precursormass.drop(index=tmpindexvalues)
            self.intensity=self.intensity.drop(index=tmpindexvalues)
            self.fragmentmass = self.fragmentmass.reset_index(drop=True)
            self.precursormass = self.precursormass.reset_index(drop=True)
            self.intensity = self.intensity.reset_index(drop=True)
            print("imported "+str(len(self.fragmentmass))+" peaks, with "+str(len(tmpindexvalues))+" removed")
        self.peaks[0]=self.fragmentmass
        self.peaks[1]=self.precursormass
        self.peaks[2]=self.intensity
        self.contourplot=peaksort.TwoDMS_plot(self.peaks[0],self.peaks[1],self.peaks[2])
        self.contourcanvas = FigureCanvasTkAgg(self.contourplot, master=self.contourframe)
        self.contourcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0)
        self.contourtoolbar = NavigationToolbar2Tk(self.contourcanvas,self.contourframe, pack_toolbar=False)
        self.contourtoolbar.update()
        self.contourtoolbar.grid(row=0,column=1,padx=0,pady=0)
        self.neutralcontourplot=peaksort.TwoDMS_plot_neutral_loss(self.peaks[0],self.peaks[1],self.peaks[2])
        self.neutralcontourcanvas = FigureCanvasTkAgg(self.neutralcontourplot, master=self.neutralcontourframe)
        self.neutralcontourcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0)
        self.neutralcontourtoolbar = NavigationToolbar2Tk(self.neutralcontourcanvas,self.neutralcontourframe, pack_toolbar=False)
        self.neutralcontourtoolbar.update()
        self.neutralcontourtoolbar.grid(row=0,column=1,padx=0,pady=0)
        plt.close()
        self.indictorcanvas.itemconfig(self.indicator, fill="green")
        self.nameprint=tk.Label(self.precursorvalues,text=self.infile)
        self.nameprint.grid(row=6,column=0,padx=5,pady=5)
        return self.fragmentmass, self.precursormass, self.intensity
    
    def dummy(self):
        print("you've run a dummy")
        
    def startsort(self):
        self.indictorcanvas.itemconfig(self.indicator, fill="red")
        if not self.peaks:
            print ("you haven't imported any data")
            return
        self.masserror=float(self.MFentry.get())
        self.refmass=float(self.refmassentry.get())
        self.refjitter=float(self.refjitterentry.get())
        self.xyres=float(self.xyresentry.get())
        self.yaxisround=float(self.yaxisroundentry.get())
        self.lowmass=float(self.lowmassentry.get())
        self.highmass=float(self.highmassentry.get())
        self.indictorcanvas.itemconfig(self.indicator, fill="red")
        self.prelistpaired=[]
        self.prepeaks=[[],[]]
        self.sorted_prelist=[]
        self.index_peaks=[]
        self.isotopeindexlist=[]
        self.isotopelist=[]
        self.preplot=plt.figure()
        for item in self.premasslist.get_children():
            self.premasslist.delete(item)
        for item in self.selectpremass.get_children():
            self.selectpremass.delete(item)
        for item in self.selectisotope.get_children():
            self.selectisotope.delete(item)
        print("starting to sort")
        self.indictorcanvas.itemconfig(self.indicator, fill="red")
        self.pre_peaklist=peaksort.identify_precursors(self.peaks[0],self.peaks[1],self.peaks[2], self.refmass, self.refjitter, self.masserror ,self.xyres)
        print(str(len(self.pre_peaklist[0]))+" precursors identified")
        self.prepeaks[0]=self.pre_peaklist[0]
        self.prepeaks[1]=self.pre_peaklist[1]
        for i in range (0,len(self.prepeaks[0])):
            self.prelistpaired.append([(self.prepeaks[0][i]),(self.prepeaks[1][i])])
            self.sorted_prelist=sorted(self.prelistpaired, key=lambda x: x[0], reverse=True)
        for i in range (0,len(self.prepeaks[0])):
            self.premasslist.insert("","end",values=(self.sorted_prelist[i]))
            self.selectpremass.insert("","end",values=(self.sorted_prelist[i]))
            self.selectisotope.insert("","end",values=(self.sorted_prelist[i]))
        self.preplot=peaksort.precursor_plotter(self.pre_peaklist[0],self.pre_peaklist[1],self.lowmass,self.highmass,self.yaxisround)
        self.precanvas = FigureCanvasTkAgg(self.preplot, master=self.precursorplots)  
        self.precanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0)
        print("precursors plotted")
        self.assignnonassign=peaksort.peak_sort(self.pre_peaklist[0],self.peaks[0], self.peaks[1], self.peaks[2],self.pre_peaklist[2])
        self.data=self.assignnonassign[0]
        self.indictorcanvas.itemconfig(self.indicator, fill="red")
        for i in range (0,len(self.prepeaks[0])):
            #print(self.data[i][0][0])
            entry=(self.data[i][0][0])
            entry=float(entry)
            self.index_peaks.append(entry)
        print("sorted peaks")
        self.counter=0
        self.comparisonlist=[]
        plt.close()
        self.indictorcanvas.itemconfig(self.indicator, fill="green")
        return self.assignnonassign[0]
        
    def new_preplot_param(self):
        self.xyres=float(self.xyresentry.get())
        self.yaxisround=float(self.yaxisroundentry.get())
        self.lowmass=float(self.lowmassentry.get())
        self.highmass=float(self.highmassentry.get())
        if not self.prepeaks:
            print("you have to sort the peaks first")
            return
        self.preplot=peaksort.precursor_plotter(self.pre_peaklist[0],self.pre_peaklist[1],self.lowmass,self.highmass,self.yaxisround)
        self.precanvas = FigureCanvasTkAgg(self.preplot, master=self.precursorplots)  
        self.precanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0)
        plt.close()
        print("precursors plotted")
        
    def selected_precusor(self,event):
        self.msmslistpaired=[]
        self.sorted_msmslist=[]
        self.yaxisround_msms=float(self.roundtomsms.get())
        self.lowmass_msms=float(self.lowmassmsms.get())
        self.highmass_msms=float(self.highmassmsms.get())
        for item in self.msmspeaks.get_children():
            self.msmspeaks.delete(item)
        item = self.selectpremass.selection()
        selection=self.selectpremass.item(item,"values")[0]
        self.selection=float(selection)
        self.index=self.index_peaks.index(self.selection)
        self.msmsplotparam=peaksort.msms_peaklister(self.data,self.index)
        for i in range (0,len(self.msmsplotparam[0])):
            self.msmslistpaired.append([(self.msmsplotparam[0][i]),(self.msmsplotparam[2][i]),(self.msmsplotparam[1][i])])
            self.sorted_msmslist=sorted(self.msmslistpaired, key=lambda x: x[0], reverse=True)
        for i in range (0,len(self.msmsplotparam[0])):
            self.msmspeaks.insert("","end",values=(self.sorted_msmslist[i]))
        self.msmsplot=peaksort.msms_plotter(self.msmsplotparam[0], self.msmsplotparam[1], self.msmsplotparam[2],self.selection, self.lowmass_msms, self.highmass_msms, self.yaxisround_msms)
        self.msmscanvas = FigureCanvasTkAgg(self.msmsplot, master=self.msmsplots)  
        self.msmscanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0)
        self.minipreplot=peaksort.mini_precursor_plotter(self.index,self.pre_peaklist[0],self.pre_peaklist[1],self.lowmass,self.highmass,self.yaxisround)
        self.miniprecanvas = FigureCanvasTkAgg(self.minipreplot, master=self.interprecursor)  
        self.miniprecanvas.get_tk_widget().grid(row=4, column=0,padx=0,pady=0)
        plt.close()
        return self.sorted_msmslist
    
    def selected_isotope(self,event):
        for item in self.isotopepeaks.get_children():
            self.isotopepeaks.delete(item)
        item = self.selectisotope.selection()
        selectioniso=self.selectisotope.item(item,"values")[0]
        self.selectioniso=float(selectioniso)
        self.index=self.index_peaks.index(self.selectioniso)
        self.isotopeindexlist.append(self.index)
        self.isotopelist.append(self.selectioniso)
        for i in range (0,len(self.isotopelist)):
            self.isotopepeaks.insert("","end",values=(self.isotopelist[i]))
            
    def delete_isotope(self,event):
        item = self.isotopepeaks.selection()
        selectiondeliso=self.isotopepeaks.item(item,"values")[0]
        self.selectiondeliso=float(selectiondeliso)
        self.index=self.index_peaks.index(self.selectiondeliso)
        self.isotopeindexlist.remove(self.index)
        self.isotopelist.remove(self.selectiondeliso)
        for item in self.isotopepeaks.get_children():
            self.isotopepeaks.delete(item)
        for i in range (0,len(self.isotopelist)):
            self.isotopepeaks.insert("","end",values=(self.isotopelist[i]))
     
    def plot_isotopes(self):
        self.isotopepeak=[]
        self.yaxisround_iso=float(self.roundtoiso.get())
        self.lowmass_iso=float(self.lowmassiso.get())
        self.highmass_iso=float(self.highmassiso.get())
        for item in self.isotopepeaklist.get_children():
            self.isotopepeaklist.delete(item)
        self.isodata=peaksort.isotope_msms_peaklister(self.data,self.isotopeindexlist,self.pre_peaklist[2])
        self.isotopeplot=peaksort.isotope_msms_plotter(self.isodata[0], self.isodata[1], self.isodata[2], self.lowmass_iso, self.highmass_iso, self.isotopeindexlist, self.yaxisround_iso)
        self.isotopecanvas = FigureCanvasTkAgg(self.isotopeplot, master=self.isotopeplotframe)  
        self.isotopecanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0,columnspan=3) 
        for i in range (0,len(self.isodata[0])):
            self.isotopepeak.append([self.isodata[0][i],self.isodata[2][i],self.isodata[1][i]])
        self.sorted_isotopepeaks=sorted(self.isotopepeak, key=lambda x: x[0], reverse=True)
        for i in range (0,len(self.isodata[0])):
            self.isotopepeaklist.insert("","end",values=(self.sorted_isotopepeaks[i]))
        plt.close()
    
    def save_isotopefigure(self):
        print("you have clicked to save your combined isotope figure")
        if self.savetypeiso.get() =="SVG":
            self.isotopesvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
            self.isotopeplot.savefig(self.isotopesvgfile)
        if self.savetypeiso.get() =="PDF":
            self.isotopepdffile=filedialog.asksaveasfilename(defaultextension=".pdf",filetypes=[("PDF Files", "*.pdf")])
            self.isotopeplot.savefig(self.isotopepdffile)
        if self.savetypeiso.get() =="JPG":
            self.isotopejpgfile=filedialog.asksaveasfilename(defaultextension=".jpg",filetypes=[("Picture Files", "*.jpg")])
            self.isotopeplot.savefig(self.isotopejpgfile)
        if self.savetypeiso.get() =="PNG":
            self.isotopepngfile=filedialog.asksaveasfilename(defaultextension=".png",filetypes=[("Picture Files", "*.jpg")])
            self.isotopeplot.savefig(self.isotopepngfile)
        plt.close()    
    
    def save_isotopepeaklist(self):
        print("you have clicked to save your combined isotope peaklist")
        self.isotopecsvfile=filedialog.asksaveasfile(defaultextension=".csv",filetypes=[("CSV Files", "*.csv")])
        pd.DataFrame(self.sorted_isotopepeaks).to_csv(self.isotopecsvfile, header=['Fragment m/z','intensity',"precursor m/z"])
        self.isotopecsvfile.close()
        print("saved")
        
    def clearisotope(self):
        self.isotopepeak=[]
        self.isodata=[]
        self.isotopeindexlist=[]
        self.isotopelist=[]
        self.isotopepeak=[]
        for item in self.isotopepeaklist.get_children():
            self.isotopepeaklist.delete(item)
        self.isotopeplot=plt.figure(figsize=(12,8))
        self.isotopecanvas = FigureCanvasTkAgg(self.isotopeplot, master=self.isotopeplotframe)  
        self.isotopecanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0,columnspan=3)
        for item in self.isotopepeaks.get_children():
            self.isotopepeaks.delete(item)
 
    def save_prefigure(self):
        print("you have clicked to save your precursor figure")
        if self.savetypepre.get() =="SVG":
            self.presvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
            self.preplot.savefig(self.presvgfile)
        if self.savetypepre.get() =="PDF":
            self.prepdffile=filedialog.asksaveasfilename(defaultextension=".pdf",filetypes=[("PDF Files", "*.pdf")])
            self.preplot.savefig(self.prepdffile)
        if self.savetypepre.get() =="JPG":
            self.prejpgfile=filedialog.asksaveasfilename(defaultextension=".jpg",filetypes=[("Picture Files", "*.jpg")])
            self.preplot.savefig(self.prejpgfile)
        if self.savetypepre.get() =="PNG":
            self.prepngfile=filedialog.asksaveasfilename(defaultextension=".png",filetypes=[("Picture Files", "*.jpg")])
            self.preplot.savefig(self.prepngfile)
            
    def save_prelist(self):
        print("you have clicked to save your precursor masslist")
        self.precsvfile=filedialog.asksaveasfile(defaultextension=".csv",filetypes=[("CSV Files", "*.csv")])
        pd.DataFrame(self.sorted_prelist).to_csv(self.precsvfile, header=['precursor m/z','intensity'])
        self.precsvfile.close()
        print("saved")
        
    def save_msmsfigure(self):
        print("you have clicked to save your MSMS figure")
        if self.savetypemsms.get() =="SVG":
            self.msmssvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
            self.msmsplot.savefig(self.msmssvgfile)
        if self.savetypemsms.get() =="PDF":
            self.msmspdffile=filedialog.asksaveasfilename(defaultextension=".pdf",filetypes=[("PDF Files", "*.pdf")])
            self.msmsplot.savefig(self.msmspdffile)
        if self.savetypemsms.get() =="JPG":
            self.msmsjpgfile=filedialog.asksaveasfilename(defaultextension=".jpg",filetypes=[("Picture Files", "*.jpg")])
            self.msmsplot.savefig(self.msmsjpgfile)
        if self.savetypemsms.get() =="PNG":
            self.msmspngfile=filedialog.asksaveasfilename(defaultextension=".png",filetypes=[("Picture Files", "*.jpg")])
            self.msmsplot.savefig(self.msmspngfile)
        plt.close()
        
    def save_msmslist(self):
        self.sorted_msmslist1d=[]
        print("you have clicked to save your MSMS masslist")
        for i in range (0,len(self.sorted_msmslist)):
            self.sorted_msmslist1d.append([self.sorted_msmslist[i][0],self.sorted_msmslist[i][1]])
        self.msmscsvfile=filedialog.asksaveasfile(defaultextension=".csv",filetypes=[("CSV Files", "*.csv")])
        pd.DataFrame(self.sorted_msmslist1d).to_csv(self.msmscsvfile, header=['m/z','intensity'])
        self.msmscsvfile.close()
        print("saved")   
        
    def add_to_compare(self):
        self.counter+=1
        self.comparisonlist.append(self.index)
        self.plotheight=14/self.counter
        self.ymaxcompare=float(self.ymaxcompareentry.get())
        self.lowmasscompare=float(self.lowmasscompareentry.get())
        self.highmasscompare=float(self.highmasscompareentry.get())
        if self.counter != 1:
            fig, axs=plt.subplots(self.counter,dpi=75)
            fig.set_figheight(11)
            fig.set_figwidth(19.5)
            fig.tight_layout(pad=5.0)
            for i in range (0,self.counter):
                self.msmsplotparam=peaksort.msms_peaklister(self.data,self.comparisonlist[i])
                axs[i].vlines(x = self.msmsplotparam[0], ymin = 0, ymax = self.msmsplotparam[2], colors="black")
                axs[i].set_ylim([0, float(self.ymaxcompare)])
                axs[i].set_xlim([self.lowmasscompare, self.highmasscompare])
                axs[i].set_xlabel("m/z", size =20,style="italic")
                axs[i].set_ylabel("Intensity", size=20)
                axs[i].set_title(str(self.msmsplotparam[1][0]), size=20)
            self.compareplot=fig
            self.comparecanvas = FigureCanvasTkAgg(self.compareplot, master=self.comparisonframe)  
            self.comparecanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0)
            
    def clearcomparison(self):
        self.counter=0
        self.comparisonlist=[]
        self.comparecanvas.delete(self.comparecanvas)
        
    def replot_compare(self):
        self.ymaxcompare=float(self.ymaxcompareentry.get())
        self.lowmasscompare=float(self.lowmasscompareentry.get())
        self.highmasscompare=float(self.highmasscompareentry.get())
        if self.counter != 1:
            fig, axs=plt.subplots(self.counter,dpi=75)
            fig.set_figheight(11)
            fig.set_figwidth(19.5)
            fig.tight_layout(pad=5.0)
            for i in range (0,self.counter):
                self.msmsplotparam=peaksort.msms_peaklister(self.data,self.comparisonlist[i])
                axs[i].vlines(x = self.msmsplotparam[0], ymin = 0, ymax = self.msmsplotparam[2], colors="black")
                axs[i].set_ylim([0, float(self.ymaxcompare)])
                axs[i].set_xlim([self.lowmasscompare, self.highmasscompare])
                axs[i].set_xlabel("m/z", size =20,style="italic")
                axs[i].set_ylabel("Intensity", size=20)
                axs[i].set_title(str(self.msmsplotparam[1][0]), size=20)
            self.compareplot=fig
            self.comparecanvas = FigureCanvasTkAgg(self.compareplot, master=self.comparisonframe)  
            self.comparecanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0)
        plt.close()
        
    def save_comparefigure(self):
        print("you have clicked to save your comparison figure")
        if self.savetypecompare.get() =="SVG":
            self.comparesvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
            self.compareplot.savefig(self.comparesvgfile)
        if self.savetypecompare.get() =="PDF":
            self.comparepdffile=filedialog.asksaveasfilename(defaultextension=".pdf",filetypes=[("PDF Files", "*.pdf")])
            self.compareplot.savefig(self.comparepdffile)
        if self.savetypecompare.get() =="JPG":
            self.comparejpgfile=filedialog.asksaveasfilename(defaultextension=".jpg",filetypes=[("Picture Files", "*.jpg")])
            self.compareplot.savefig(self.comparejpgfile)
        if self.savetypecompare.get() =="PNG":
            self.comparepngfile=filedialog.asksaveasfilename(defaultextension=".png",filetypes=[("Picture Files", "*.jpg")])
            self.compareplot.savefig(self.comparepngfile)
        plt.close()
        
    def highlight_ridges(self):
        self.ridges=peaksort.identify_ridges(self.msmsplotparam[0], self.msmsplotparam[1], self.msmsplotparam[2], self.selection, self.pre_peaklist[0], self.refmass, self.refjitter, ppm=0)
        self.msmsplot=peaksort.highlight_msms_plotter(self.msmsplotparam[0], self.msmsplotparam[1], self.msmsplotparam[2],  self.selection, self.lowmass_msms, self.highmass_msms, self.ridges, self.yaxisround_msms)
        if self.ridges != []:
            self.msmscanvas = FigureCanvasTkAgg(self.msmsplot, master=self.msmsplots)  
            self.msmscanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0)
        plt.close()
        
    def remove_ridges(self):
        for item in self.msmspeaks.get_children():
            self.msmspeaks.delete(item)
        self.ridgeindex=[]
        self.ridgelesspeak=[]
        self.tmpx=self.msmsplotparam[0]
        self.tmpy=self.msmsplotparam[1]
        self.tmpz=self.msmsplotparam[2]
        for i in self.ridges[0]:
            tmpridge=self.msmsplotparam[0].index(i)
            if tmpridge not in self.ridgeindex:
                self.ridgeindex.append(tmpridge)
        self.sortridgeindex=sorted(self.ridgeindex, key=lambda x: x, reverse=True)
        for i in self.sortridgeindex:
            del self.tmpx[i]
            del self.tmpy[i]
            del self.tmpz[i]
        self.msmsplot=peaksort.msms_plotter(self.tmpx, self.tmpy, self.tmpz,self.selection, self.lowmass_msms, self.highmass_msms, self.yaxisround_msms)
        self.msmscanvas = FigureCanvasTkAgg(self.msmsplot, master=self.msmsplots)  
        self.msmscanvas.get_tk_widget().grid(row=0, column=0,padx=0,pady=0)
        for i in range (0,len(self.tmpx)):
            self.ridgelesspeak.append([self.tmpx[i],self.tmpz[i],self.tmpy[i]])
        self.sorted_ridgelesspeaks=sorted(self.ridgelesspeak, key=lambda x: x[0], reverse=True)
        for i in range (0,len(self.tmpx)):
            self.msmspeaks.insert("","end",values=(self.sorted_ridgelesspeaks[i]))
        plt.close()
            
    def sort_by_intensity(self):
        for item in self.selectpremass.get_children():
            self.selectpremass.delete(item)
        self.prelistpairedi=[]
        for i in range (0,len(self.prepeaks[0])):
            self.prelistpairedi.append([(self.prepeaks[0][i]),(self.prepeaks[1][i])])
            self.sortedbyi=sorted(self.prelistpairedi, key=lambda x: x[1], reverse=True)
        for i in range (0,len(self.prepeaks[0])):
            self.selectpremass.insert("","end",values=(self.sortedbyi[i]))
        plt.close()
            
    def sort_by_mz(self):
        for item in self.selectpremass.get_children():
            self.selectpremass.delete(item)
        for i in range (0,len(self.prepeaks[0])):
            self.selectpremass.insert("","end",values=(self.sorted_prelist[i]))
            
    def export_as(self):
        if self.exporttype.get() =="MASCOT":
            self.mgffile=filedialog.asksaveasfilename(defaultextension=".mgf",filetypes=[("MGF Files", "*.mgf")])
            self.mgftext=peaksort.mgf_text(self.assignnonassign[0],"Title")
            print("converted")
            with open(self.mgffile, 'w') as file:
                file.write(self.mgftext)
            print("saved")
            
        if self.exporttype.get() =="mzML":
            self.precursor_info=[[],[]]
            self.analyser_type=self.analysertype.get()
            self.source_type=self.sourcetype.get()
            self.Polarity=self.polaritytype.get()
            self.fragmentation_technique=self.fragmentation_techniquetype.get()
            self.fragmentation_power=self.export_fragmentation_power.get()
            self.Fragmentation_duration=self.export_Fragmentation_duration.get()
            self.MS_level=self.MS_leveltype.get()
            self.ExD_bias=self.export_ExD_bias.get()
            self.ExD_lens=self.export_ExD_lens.get()
            self.precursor_info=self.pre_peaklist
            #for i in range(0,len(self.precursor_info)):
                #print(self.precursor_info[1][i])
            self.mzmlfile=exportasmzml.write_to()
            self.exportmz=exportasmzml.writemzml(self.mzmlfile,self.assignnonassign[0],self.source_type,
                                                 self.analyser_type,self.Polarity,self.fragmentation_technique,
                                                 self.fragmentation_power,self.Fragmentation_duration,
                                                 self.MS_level,self.ExD_bias, self.ExD_lens,self.precursor_info)
            print(self.exportmz)
            if self.exp_text==1:
                self.expfile=filedialog.asksaveasfilename(defaultextension=".txt",filetypes=[("TEXT Files", "*.txt")])
                self.exptext=peaksort.experiment_txt_file(self.data,100,2000)
                with open(self.expfile, 'w') as file:
                    file.write(self.exptext)
                print("saved")
                
    def single_export(self):
        if self.singleexporttype.get() =="Sirius":
            self.compound_name=self.compound.get()
            self.charge=self.chargeinput.get()
            self.frag_type=self.frag_input.get()
            self.siriustext=peaksort.sirius_text(self.data,self.compound_name,self.charge,self.isotopeindexlist,
                                                 self.xyres,self.isodata[0], self.isodata[1], self.isodata[2],
                                                 self.frag_type)
            self.siriusfile=filedialog.asksaveasfilename(defaultextension=".ms",filetypes=[("MS files", "*.ms")])
            print("converted")
            with open(self.siriusfile, 'w') as file:
                file.write(self.siriustext)
            print("saved")
            
        if self.singleexporttype.get() =="mzML":
            print("Anna hasn't tested this yet :)")
            self.analyser_type=self.analysertype.get()
            self.source_type=self.sourcetype.get()
            self.Polarity=self.polaritytype.get()
            self.fragmentation_technique=self.fragmentation_techniquetype.get()
            self.fragmentation_power=self.export_fragmentation_power.get()
            self.Fragmentation_duration=self.export_Fragmentation_duration.get()
            self.MS_level=self.MS_leveltype.get()
            self.ExD_bias=self.export_ExD_bias.get()
            self.ExD_lens=self.export_ExD_lens.get()
            self.precursor_int=[0,0]
            self.mzmlfile=exportasmzml.write_to()
            self.exportmz=exportasmzml.writesinglemzml(self.mzmlfile,self.sorted_isotopepeaks,self.source_type,
                                                       self.analyser_type,self.Polarity,
                                                       self.fragmentation_technique,self.fragmentation_power,
                                                       self.Fragmentation_duration,self.MS_level,self.ExD_bias,
                                                       self.ExD_lens,self.precursor_int)
            
    def exp_file_LipidDial(self):
        self.expfile=filedialog.asksaveasfilename(defaultextension=".txt",filetypes=[("TEXT Files", "*.txt")])
        self.exptext=peaksort.experiment_txt_file(self.data,min(self.precursormass),max(self.precursormass))
        with open(self.expfile, 'w') as file:
            file.write(self.exptext)
        print("created a experiment file")
        
    def NLVLsortandplot(self):
        self.NL=float(self.NLEntry.get())
        self.VL=float(self.VLEntry.get())
        self.refyaxisjitter=float(self.refjitterentry.get())
        self.refdiagonaljitter=float(self.xyresentry.get())
        self.refmasserror=float(self.MFentry.get())
        self.refmass2=float(self.refmassentry.get())
        self.lowmassNLVL=float(self.lowmassNLVLentry.get())
        self.highmassNLVL=float(self.highmassNLVLentry.get())
        self.yaxisroundNLVL=float(self.yaxisroundNLVLentry.get())
        self.NLVLlist=[[],[],[]]
        for item in self.NLVLpeaks.get_children():
            self.NLVLpeaks.delete(item)
        if self.NLVLtype.get() =="VL":
            self.yaxisjitter=peaksort.mass_jitter(self.refmass2,self.refyaxisjitter,self.refmasserror,self.VL)
            print("Looking for a vertical line at "+str(self.VL))
            for i in range (0,len(self.peaks[0])):
                difference=self.peaks[0][i]-self.VL
                if abs(difference) <= self.yaxisjitter:
                    self.NLVLlist[0].append(self.peaks[0][i])
                    self.NLVLlist[1].append(self.peaks[1][i])
                    self.NLVLlist[2].append(self.peaks[2][i])
            #♣print(len(self.NLVLlist),self.NLVLlist[1])
            for i in range (0,len(self.NLVLlist[0])):
                self.NLVLpeaks.insert("","end",values=(self.NLVLlist[0][i],
                                                       self.NLVLlist[1][i],self.NLVLlist[2][i]))
            self.NLVLplot=peaksort.VL_plotter(self.NLVLlist[1], self.NLVLlist[0], self.NLVLlist[2],
                                              self.VL, self.lowmassNLVL, self.highmassNLVL, roundto=self.yaxisroundNLVL)
            self.NLVLcanvas = FigureCanvasTkAgg(self.NLVLplot, master=self.NLVLplotter)  
            self.NLVLcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0,columnspan=3) 
        elif self.NLVLtype.get() =="NL":
            print("Looking for a neutral loss line at "+str(self.NL))
            for i in range (0,len(self.peaks[0])):
                self.yaxisjitter=peaksort.mass_jitter(self.refmass2,self.refdiagonaljitter,self.refmasserror,self.peaks[0][i])
                difference=(self.peaks[1][i]-self.peaks[0][i])-self.NL
                if abs(difference) <= self.yaxisjitter:
                    self.NLVLlist[0].append(self.peaks[0][i])
                    self.NLVLlist[1].append(self.peaks[0][i]+self.NL)
                    self.NLVLlist[2].append(self.peaks[2][i])
            for i in range (0,len(self.NLVLlist[0])):
                self.NLVLpeaks.insert("","end",values=(self.NLVLlist[0][i],
                                                       self.NLVLlist[1][i],self.NLVLlist[2][i]))
            self.NLVLplot=peaksort.NL_plotter(self.NLVLlist[1], self.NLVLlist[0], self.NLVLlist[2],
                                              self.NL, self.lowmassNLVL, self.highmassNLVL, roundto=self.yaxisroundNLVL)
            self.NLVLcanvas = FigureCanvasTkAgg(self.NLVLplot, master=self.NLVLplotter)  
            self.NLVLcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0,columnspan=3) 
            #self.refyaxisjitter=peaksort.mass_jitter(self.refmass2,self.refyaxisjitter,self.refmasserror,)
            
    def save_NLVLfigure(self):
        print("you have clicked to save your Neutral loss or vertical line figure")
        if self.savetypeNLVL.get() =="SVG":
            self.NLVLsvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
            self.NLVLplot.savefig(self.NLVLsvgfile)
        if self.savetypeNLVL.get() =="PDF":
            self.NLVLpdffile=filedialog.asksaveasfilename(defaultextension=".pdf",filetypes=[("PDF Files", "*.pdf")])
            self.NLVLplot.savefig(self.NLVLpdffile)
        if self.savetypeNLVL.get() =="JPG":
            self.NLVLjpgfile=filedialog.asksaveasfilename(defaultextension=".jpg",filetypes=[("Picture Files", "*.jpg")])
            self.NLVLplot.savefig(self.NLVLjpgfile)
        if self.savetypeNLVL.get() =="PNG":
            self.NLVLpngfile=filedialog.asksaveasfilename(defaultextension=".png",filetypes=[("Picture Files", "*.jpg")])
            self.NLVLplot.savefig(self.NLVLpngfile)
        plt.close()
        
    def save_NLVLlist(self):
        self.sorted_NLVLlist1d=[]
        print("you have clicked to save your Neutral loss or vertical line masslist")
        for i in range (0,len(self.NLVLlist[0])):
            self.sorted_NLVLlist1d.append([self.NLVLlist[0][i],self.NLVLlist[1][i],self.NLVLlist[2][i]])
        self.NLVLcsvfile=filedialog.asksaveasfile(defaultextension=".csv",filetypes=[("CSV Files", "*.csv")])
        pd.DataFrame(self.sorted_NLVLlist1d).to_csv(self.NLVLcsvfile, header=['fragment m/z','precursor m/z','intensity'])
        self.NLVLcsvfile.close()
        print("saved")

if __name__=="__main__":
    r=TwoDpeaksortergui()