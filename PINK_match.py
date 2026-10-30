# -*- coding: utf-8 -*-
"""
Created on Thu May 9 13:21:02 2024

@author: lwg513
"""

"""Peak INterpretation via picKing (P.I.N.K) Match for spectral matching and interpretation of MS2 spectra extracted from 2DMS peak lists
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

import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter import filedialog
from scripts.peaksorter import peaksort
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
from pymzml import run #M Kösters, J Leufken, S Schulze, K Sugimoto, J Klein, R P Zahedi, M Hippler, S A Leidel, C Fufezan; pymzML v2.0: introducing a highly compressed and seekable gzip format, Bioinformatics, doi: https://doi.org/10.1093/bioinformatics/bty046
from scripts.import_DBs import import_DB
#import csv
pythonicon = "resources\\PINK_logo_transparent.gif"
from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import Draw #RDKit: Open-source cheminformatics. https://www.rdkit.org https://doi.org/10.5281/zenodo.591637

class read_mzML():
    def __init__(self):
        self.root=tk.Tk()
        self.root.title("P.I.N.K match")
        self.root.geometry("1600x850")
        
        photo_image = tk.PhotoImage(file=pythonicon)
        #print(type(photo_image))
        self.root.wm_iconphoto(False, photo_image)
        
        style=ttk.Style(self.root)
        style.configure("TFrame", background="#E69", foreground="white")
        
        nb = ttk.Notebook(self.root)
        self.mzML_viewer = ttk.Frame(nb)
        nb.add(self.mzML_viewer, text = "View MSn data")
        self.DB_results = ttk.Frame(nb)
        nb.add(self.DB_results, text = "View assignments")
        #self.small_batch = ttk.Frame(nb)
        #nb.add(self.small_batch, text = "small batch")
        self.results_viewer = ttk.Frame(nb)
        nb.add(self.results_viewer, text = "View results")
        
        nb.pack(expand=1, fill="both")
        self.importcounter=0
        self.import_mzML=tk.Button(self.mzML_viewer, text="Import mzML files", 
                                   command=self.open_mzML_file)
        self.import_mzML.grid(row=0, column=0, padx=5,pady=5)
        
        self.lowcutoffentry=tk.Entry(self.mzML_viewer)
        self.lowcutoffentry.insert(0, "0") 
        self.lowcutoffentry.grid(row=2, column=0, padx=5,pady=5)
        
        self.highcutoffentry=tk.Entry(self.mzML_viewer)
        self.highcutoffentry.insert(0, "3000") 
        self.highcutoffentry.grid(row=3, column=0, padx=5,pady=5)
        
        self.precursor_cutoff=tk.Button(self.mzML_viewer, text="Do mass cut off of precursors", 
                                   command=self.mzml_masscutoff)
        self.precursor_cutoff.grid(row=4, column=0, padx=5,pady=5)
        
        self.msmsplot=plt.figure(figsize=(12,8),dpi=75)
        self.msmscanvas = FigureCanvasTkAgg(self.msmsplot, master=self.mzML_viewer)  
        self.msmscanvas.get_tk_widget().grid(row=1,column=1,padx=0,pady=0,columnspan=3,rowspan=3)
        
        #self.msmsplotsb=plt.figure(figsize=(12,8),dpi=75)
        #self.msmscanvassb = FigureCanvasTkAgg(self.msmsplotsb, master=self.small_batch)  
        #self.msmscanvassb.get_tk_widget().grid(row=1,column=1,padx=0,pady=0,columnspan=3,rowspan=5)
        
        self.mirrorplot=plt.figure(figsize=(12,8),dpi=75)
        self.mirrorcanvas = FigureCanvasTkAgg(self.mirrorplot, master=self.DB_results)  
        self.mirrorcanvas.get_tk_widget().grid(row=1,column=1,padx=0,pady=0,columnspan=3,rowspan=3)
        
        self.mirrorplotresults=plt.figure(figsize=(12,8),dpi=75)
        self.mirrorcanvasresults = FigureCanvasTkAgg(self.mirrorplotresults, master=self.results_viewer)  
        self.mirrorcanvasresults.get_tk_widget().grid(row=1,column=1,padx=0,pady=0,columnspan=3,rowspan=3)
        
        self.yaxisroundlabel=tk.Label(self.mzML_viewer,text="what would you like y axis rounded to?")
        self.yaxisroundlabel.grid(row=5, column=1, padx=5,pady=5 )
        self.yaxisroundentry=tk.Entry(self.mzML_viewer)
        self.yaxisroundentry.insert(0, "1e8") 
        self.yaxisroundentry.grid(row=4, column=1, padx=5,pady=5 )

        self.lowmasslabel=tk.Label(self.mzML_viewer,text="lowmass in x")
        self.lowmasslabel.grid(row=5, column=2, padx=5,pady=5 )
        self.lowmassentry=tk.Entry(self.mzML_viewer)
        self.lowmassentry.insert(0, "100") 
        self.lowmassentry.grid(row=4, column=2, padx=5,pady=5 )

        self.highmasslabel=tk.Label(self.mzML_viewer,text="highmass")
        self.highmasslabel.grid(row=5, column=3, padx=5,pady=5 )
        self.highmassentry=tk.Entry(self.mzML_viewer)
        self.highmassentry.insert(0, "1500") 
        self.highmassentry.grid(row=4, column=3, padx=5,pady=5 )
        
        #self.yaxisroundlabelsb=tk.Label(self.small_batch,text="what would you like y axis rounded to?")
        #self.yaxisroundlabelsb.grid(row=5, column=1, padx=5,pady=5 )
        #self.yaxisroundentrysb=tk.Entry(self.small_batch)
        #self.yaxisroundentrysb.insert(0, "1e8") 
        #self.yaxisroundentrysb.grid(row=4, column=1, padx=5,pady=5 )

        #self.lowmasslabelsb=tk.Label(self.small_batch,text="lowmass in x")
        #self.lowmasslabelsb.grid(row=5, column=2, padx=5,pady=5 )
        #self.lowmassentrysb=tk.Entry(self.small_batch)
        #self.lowmassentrysb.insert(0, "100") 
        #self.lowmassentrysb.grid(row=4, column=2, padx=5,pady=5 )

        #self.highmasslabelsb=tk.Label(self.small_batch,text="highmass")
        #self.highmasslabelsb.grid(row=5, column=3, padx=5,pady=5 )
        #self.highmassentrysb=tk.Entry(self.small_batch)
        #self.highmassentrysb.insert(0, "1500") 
        #self.highmassentrysb.grid(row=4, column=3, padx=5,pady=5 )
        

        self.lowmasslabelmirror=tk.Label(self.DB_results,text="lowmass in x")
        self.lowmasslabelmirror.grid(row=5, column=2, padx=5,pady=5 )
        self.lowmassentrymirror=tk.Entry(self.DB_results)
        self.lowmassentrymirror.insert(0, "100") 
        self.lowmassentrymirror.grid(row=4, column=2, padx=5,pady=5 )

        self.highmasslabelmirror=tk.Label(self.DB_results,text="highmass")
        self.highmasslabelmirror.grid(row=5, column=3, padx=5,pady=5 )
        self.highmassentrymirror=tk.Entry(self.DB_results)
        self.highmassentrymirror.insert(0, "1500") 
        self.highmassentrymirror.grid(row=4, column=3, padx=5,pady=5 )
        
        self.displaycossimlabel=tk.Label(self.DB_results,text="Similarity score")
        self.displaycossimlabel.grid(row=0, column=1, padx=5,pady=5 )
        self.displaycossim=tk.Label(self.DB_results,text="0")
        self.displaycossim.grid(row=0, column=2, padx=5,pady=5 )
        self.indictorcanvas = tk.Canvas(self.DB_results,bg="#E69",height=50,width=50,highlightthickness=0)
        self.indictorcanvas.grid(row=0, column=3,padx=5,pady=5,rowspan=1)
        self.indicator = self.indictorcanvas.create_oval(20, 10, 10, 20)
        self.indictorcanvas.itemconfig(self.indicator, fill="red")
        
        self.smiles=plt.figure(figsize=(5,5),dpi=75)
        self.SMILEScanvas = FigureCanvasTkAgg(self.smiles,master=self.DB_results)
        self.SMILEScanvas.get_tk_widget().grid(row=2, column=6,padx=5,pady=5,rowspan=1)
        self.infocanvas = tk.Canvas(self.DB_results,bg="#E69",height=50,width=50,highlightthickness=0)
        self.infocanvas.grid(row=3, column=6,padx=5,pady=5,rowspan=1)
        
        self.importDBbutton=tk.Button(self.mzML_viewer, text="Import database files", 
                                   command=self.import_Database)
        self.importDBbutton.grid(row=2, column=6, padx=5,pady=5)
        self.saveDBbutton=tk.Button(self.mzML_viewer, text="Save database files as csv", 
                                   command=self.save_DB_as_csv)
        self.saveDBbutton.grid(row=3, column=6, padx=5,pady=5)
        self.clearDBs=tk.Button(self.mzML_viewer, text="Clear Database list", 
                                   command=self.clear_DB_list)
        self.clearDBs.grid(row=4, column=6, padx=5,pady=5)
        self.compareDBbutton=tk.Button(self.mzML_viewer, text="Search\ndatabases", 
                                   command=self.compare_to_DB,height= 3, width=10)
        self.compareDBbutton.grid(row=5, column=6, padx=5,pady=5)
        self.ppmerrorlabel=tk.Label(self.mzML_viewer,text="PPM error")
        self.ppmerrorlabel.grid(row=2, column=7, padx=5,pady=5)
        self.ppmerrorentry=tk.Entry(self.mzML_viewer)
        self.ppmerrorentry.insert(0, "2")
        self.ppmerrorentry.grid(row=3, column=7, padx=5,pady=5)
        self.progressbar = ttk.Progressbar(self.mzML_viewer)
        self.progressbar.grid(row=6, column=6, padx=5,pady=5,columnspan=2)
        self.removetrue = tk.IntVar()
        self.removesinglepeaks=tk.Checkbutton(self.mzML_viewer, text='Remove single fragment matches?',variable=self.removetrue, onvalue=1, offvalue=0)
        self.removesinglepeaks.grid(row=4, column=7, padx=5,pady=5,rowspan=2)
        
        #self.importDBbuttonsb=tk.Button(self.small_batch, text="Import database files", 
        #                           command=self.import_Database)
        #self.importDBbuttonsb.grid(row=4, column=6, padx=5,pady=5)
        #self.saveDBbuttonsb=tk.Button(self.small_batch, text="Save database files as csv", 
        #                           command=self.save_DB_as_csv)
        #self.saveDBbuttonsb.grid(row=5, column=6, padx=5,pady=5)
        #self.clearDBssb=tk.Button(self.small_batch, text="Clear Database list", 
        #                           command=self.clear_DB_list)
        #self.clearDBssb.grid(row=6, column=6, padx=5,pady=5)
        #self.compareDBbuttonsb=tk.Button(self.small_batch, text="Search\ndatabases", 
        #                           command=self.compare_to_DB_sb,height= 3, width=10)
        #self.compareDBbuttonsb.grid(row=7, column=6, padx=5,pady=5)
        #self.ppmerrorlabelsb=tk.Label(self.small_batch,text="PPM error")
        #self.ppmerrorlabelsb.grid(row=4, column=7, padx=5,pady=5)
        #self.ppmerrorentrysb=tk.Entry(self.small_batch)
        #self.ppmerrorentrysb.insert(0, "2")
        #self.ppmerrorentrysb.grid(row=5, column=7, padx=5,pady=5)
        #self.progressbarsb = ttk.Progressbar(self.small_batch)
        #self.progressbarsb.grid(row=8, column=6, padx=5,pady=5,columnspan=2)
        #self.removesinglepeakssb=tk.Checkbutton(self.small_batch, text='Remove single fragment matches?',variable=self.removetrue, onvalue=1, offvalue=0)
        #self.removesinglepeakssb.grid(row=6, column=7, padx=5,pady=5,rowspan=2)
        
        self.saveresults=tk.Button(self.results_viewer, text="Save results", 
                                   command=self.export_results)
        self.saveresults.grid(row=0,column=6)        
        self.precursors=ttk.Treeview(self.mzML_viewer,height=20)
        vsbp = ttk.Scrollbar(self.mzML_viewer, orient="vertical", 
                             command=self.precursors.yview)
        vsbp.grid(row=1, column =0, sticky=N+S+W)
        self.precursors.configure(yscrollcommand=vsbp.set)
        self.precursors["columns"]=("m/z")
        self.precursors.column("#0",width=0)
        self.precursors.column("m/z",width=120)
        self.precursors.heading("m/z",text="precursors m/z")
        self.precursors.grid(row=1, column=0, padx=5,pady=5)
        self.precursors.bind("<Double-1>", self.selected_precusor) 
        self.precursors.bind("<Double-3>", self.delete_msn)
        
        #self.precursorssb=ttk.Treeview(self.small_batch,height=20)
        #vsbpsb = ttk.Scrollbar(self.small_batch, orient="vertical", 
        #                     command=self.precursorssb.yview)
        #vsbpsb.grid(row=1, column =0,rowspan=3, sticky=N+S+W)
        #self.precursorssb.configure(yscrollcommand=vsbpsb.set)
        #self.precursorssb["columns"]=("m/z")
        #self.precursorssb.column("#0",width=0)
        #self.precursorssb.column("m/z",width=120)
        #self.precursorssb.heading("m/z",text="precursors m/z")
        #self.precursorssb.grid(row=1, column=0, padx=5,pady=5,rowspan=2)
        #self.precursorssb.bind("<Double-1>",self.selected_precusorsb)
        #self.precursorssb.bind("<Double-3>",self.move_precursor)
        
        self.DBprecursors=ttk.Treeview(self.DB_results,height=20)
        vsbpdb = ttk.Scrollbar(self.DB_results, orient="vertical", 
                             command=self.DBprecursors.yview)
        vsbpdb.grid(row=1, column =0, sticky=N+S+W,rowspan=3)
        self.DBprecursors.configure(yscrollcommand=vsbpdb.set)
        self.DBprecursors["columns"]=("m/z")
        self.DBprecursors.column("#0",width=0)
        self.DBprecursors.column("m/z",width=120)
        self.DBprecursors.heading("m/z",text="precursors m/z")
        self.DBprecursors.grid(row=1, column=0, padx=5,pady=5,rowspan=3)
        self.DBprecursors.bind("<Double-1>", self.selected_insert_precusor) 
        
        self.msmspeakslabel=tk.Label(self.mzML_viewer,text="tandem peaklist")
        self.msmspeakslabel.grid(row=0, column=5,padx=5,pady=5)
        self.msmspeaks=ttk.Treeview(self.mzML_viewer,height=30)
        vsbmsms = ttk.Scrollbar(self.mzML_viewer, orient="vertical", 
                                command=self.msmspeaks.yview)
        vsbmsms.grid(row=1, column =5, sticky=N+S+E,columnspan=1,rowspan=3)
        self.msmspeaks.configure(yscrollcommand=vsbmsms.set)
        self.msmspeaks["columns"]=("m/z","I")
        self.msmspeaks.column("#0",width=0)
        self.msmspeaks.column("m/z",width=100)
        self.msmspeaks.column("I",width=100)
        self.msmspeaks.heading("m/z",text="fragment m/z")
        self.msmspeaks.heading("I",text="intensity")
        self.msmspeaks.grid(row=1, column=5, padx=5,pady=5,columnspan=1,rowspan=3)
        
        self.DBlabel=tk.Label(self.mzML_viewer,text="DB list")
        self.DBlabel.grid(row=0, column=6,padx=5,pady=5)
        self.DBList=ttk.Treeview(self.mzML_viewer,height=10)
        self.DBList["columns"]=("name")
        self.DBList.column("#0",width=0)
        self.DBList.column("name",width=400)
        self.DBList.heading("name",text="Databases")
        self.DBList.grid(row=1, column=6, padx=5,pady=5,columnspan=2)
        
        self.savemirrorplotbutton=tk.Button(self.DB_results,text="save your figure",command=self.save_mirror_plot)
        self.savemirrorplotbutton.grid(row=0, column=7, padx=5,pady=5 )
        self.saveoptionsmirror=["SVG","PDF","JPG","PNG"]
        self.savetypemirror = StringVar(self.DB_results)
        self.savetypemirror.set(self.saveoptionsmirror[0])
        self.drop=OptionMenu(self.DB_results, self.savetypemirror, *self.saveoptionsmirror)
        self.drop.grid(row=1, column=7, padx=5,pady=5)
        
        #self.DBlabelsb=tk.Label(self.small_batch,text="DB list")
        #self.DBlabelsb.grid(row=0, column=6,padx=5,pady=5,columnspan=2)
        #self.DBListsb=ttk.Treeview(self.small_batch,height=10)
        #self.DBListsb["columns"]=("name")
        #self.DBListsb.column("#0",width=0)
        #self.DBListsb.column("name",width=400)
        #self.DBListsb.heading("name",text="Databases")
        #self.DBListsb.grid(row=1, column=6, padx=5,pady=5,columnspan=2)
        
        #self.sblabelsb=tk.Label(self.small_batch,text="list of precursors to run")
        #self.sblabelsb.grid(row=2, column=6)
        #self.sbList=ttk.Treeview(self.small_batch,height=10)
        #vsbsb = ttk.Scrollbar(self.small_batch, orient="vertical", 
        #                      command=self.sbList.yview)
        #vsbsb.grid(row=3, column =7, sticky=N+S+W)
        #self.sbList.configure(yscrollcommand=vsbsb.set)
        #self.sbList["columns"]=("name")
        #self.sbList.column("#0",width=0)
        #self.sbList.column("name",width=400)
        #self.sbList.heading("name",text="precursor")
        #self.sbList.grid(row=3, column=6, padx=5,pady=5,columnspan=2)
        #self.sbList.bind("<Double-1>",self.selected_precusorsb)
        #self.sbList.bind("<Double-3>",self.delete_precursorsb)
        
        self.assignmentslabel=tk.Label(self.DB_results,text="Potential matches")
        self.assignmentslabel.grid(row=0, column=6)
        self.assignments=ttk.Treeview(self.DB_results,height=10)
        vsbass = ttk.Scrollbar(self.DB_results, orient="vertical", 
                              command=self.assignments.yview)
        vsbass.grid(row=1, column =6, sticky=N+S+E, columnspan=1,rowspan=1)
        self.assignments.configure(yscrollcommand=vsbass.set)
        self.assignments["columns"]=("name","cossim")
        self.assignments.column("#0",width=0)
        self.assignments.column("name",width=200)
        self.assignments.column("cossim",width=100)
        self.assignments.heading("name",text="assignment")
        self.assignments.heading("cossim",text="Cosine\nsimilarity")
        self.assignments.grid(row=1, column=6)
        self.assignments.bind("<Double-1>", self.selected_DB_precusor)
        self.assignments.bind("<Double-3>", self.confirm_assignment)
        
        self.resultsTV=ttk.Treeview(self.results_viewer,height=20)
        vsbp = ttk.Scrollbar(self.results_viewer, orient="vertical", 
                             command=self.resultsTV.yview)
        vsbp.grid(row=1, column =0, sticky=N+S+W)
        self.resultsTV.configure(yscrollcommand=vsbp.set)
        self.resultsTV["columns"]=("m/z","cos","ass","in")
        self.resultsTV.column("#0",width=0)
        self.resultsTV.column("m/z",width=120)
        self.resultsTV.heading("m/z",text="precursors m/z")
        self.resultsTV.column("cos",width=120)
        self.resultsTV.heading("cos",text="Cosine similarity")
        self.resultsTV.column("ass",width=120)
        self.resultsTV.heading("ass",text="Asignment")
        self.resultsTV.column("in",width=120)
        self.resultsTV.heading("in",text="Y value")
        self.resultsTV.grid(row=1, column=0, padx=5,pady=5)
        self.resultsTV.bind("<Double-1>", self.show_results) 
        self.resultsTV.bind("<Double-3>", self.unconfirm_assignment)
        
        self.DBlist=[]
        self.names=[]
        self.precursormz=[]
        self.smiles=[]
        self.adducts=[]
        self.inchikey=[]
        self.formula=[]
        self.retentiontime=[]
        self.ccs=[]
        self.ionmode=[]
        self.compoundclass=[]
        self.comments=[]
        self.number_of_peaks=[]
        self.peaksindexlists=[]
        self.mzlist=[]
        self.ilist=[]
        self.sblist=[]
        
        self.postcomparisondatasetsorted=None
        self.DBinfolabel=tk.Label(self.DB_results,text="Name: \nPrecursor m/z:\nAdduct: \nINCHIKEY: \nFormula: \nIon mode:")
        
        self.root.mainloop()
        
    def dummy(self):
        print("you've run a dummy")
        
    def dummy2(self,event):
        print("you've run a dummy")
            
    def open_mzML_file(self):
        if self.importcounter!=0:
            self.displayinfile.destroy()
        self.allowedIDs=[]
        self._DB_mz_list=[]
        self._assignment_string_list=[]
        self._cossim_list=[]
        self._ID_list=[]
        self._input_mzlist=[]
        self.infile=filedialog.askopenfilename(defaultextension=".mzML", filetypes=[("mzML Files", "*.mzML")])
        self.displayinfile=tk.Label(self.mzML_viewer,text=self.infile)
        self.displayinfile.grid(row=0,column=1,columnspan=3)
        self.importcounter=self.importcounter+1
        with open(str(self.infile),"r") as f:
            self.data = f.read()
        f.close()
        for item in self.precursors.get_children():
            self.precursors.delete(item)
        for item in self.DBprecursors.get_children():
            self.DBprecursors.delete(item)
        self.run = run.Reader(self.infile)
        self.fragmented_precursors = {}
        for spectrum in self.run:
            if spectrum.ms_level >1 and len(spectrum.mz)>1:
                selected_precursors = spectrum.selected_precursors
                if spectrum.selected_precursors is not None:
                    for precursor_dict in selected_precursors:
                        precursor_mz = precursor_dict["mz"]
                        rounded_precursor_mz = round(precursor_mz, 5)
                        if rounded_precursor_mz not in self.fragmented_precursors.keys():
                            self.fragmented_precursors[rounded_precursor_mz] = []
                        self.fragmented_precursors[rounded_precursor_mz].append(spectrum.ID)
                        self.allowedIDs.append(spectrum.ID)
       # print (self.fragmented_precursors)
        for i in range (0,len(self.allowedIDs)):
            for spec in self.run:
                if self.allowedIDs[i] == spec.ID:
                    selected_precursors = spec.selected_precursors
                    for precursor_dict in selected_precursors:
                        self.precursor_mz = precursor_dict["mz"]
                        #self.precursor_int=precursor_dict["intensity"]
                    self.precursors.insert("","end",values=(round(self.precursor_mz,5)))
                    self.DBprecursors.insert("","end",values=(round(self.precursor_mz,5)))
                    #self.precursorssb.insert("","end",values=(round(self.precursor_mz,5)))
        for item in self.resultsTV.get_children():
            self.resultsTV.delete(item)
        for i in range (0,len(self._DB_mz_list)):
            self.resultsTV.insert("","end",values=((self._DB_mz_list[i]),(self._cossim_list[i]),(self._assignment_string_list[i])))
        
    
    def import_Database(self):
        self.progressbar = ttk.Progressbar(mode="indeterminate")
        self.progressbar.start()
        self.progressbar.update_idletasks()
        initial=self.mzlist
        self.DBfile=filedialog.askopenfilename(filetypes=[("MSP Files", "*.msp"),("CSV Files","*.csv"),("MSP Files", "*.MSP")])
        if (".msp") or (".MSP") in self.DBfile:
            self.DBdata=import_DB.read_MSP_2(self.DBfile,self.names,self.precursormz,self.smiles,self.adducts,self.inchikey,self.formula,self.retentiontime,self.ccs,self.ionmode,self.compoundclass,self.comments,self.number_of_peaks,self.peaksindexlists,self.mzlist,self.ilist)
            #print(self.DBdata.iloc[4,7])
            #self.DBdata["Name"]=self.DBdata["Name"].astype("category")#change object to category to save RAM
            self.DBdata["adduct"]=self.DBdata["adduct"].astype("category")
            #self.DBdata["smile"]=self.DBdata["smile"].astype("category")
            #self.DBdata["INCHIKEY"]=self.DBdata["INCHIKEY"].astype("category")
            self.DBdata["Formula"]=self.DBdata["Formula"].astype("category")
            self.DBdata["Ion mode"]=self.DBdata["Ion mode"].astype("category")
            self.DBdata["Precursor m/z"]=pd.to_numeric(self.DBdata["Precursor m/z"], downcast="float") #downcast from float64
        if (".csv") in self.DBfile:
            self.DBdata=import_DB.read_csv(self.DBfile,self.names,self.precursormz,self.smiles,self.adducts,self.inchikey,self.formula,self.retentiontime,self.ccs,self.ionmode,self.compoundclass,self.comments,self.number_of_peaks,self.peaksindexlists,self.mzlist,self.ilist)
            #print(self.DBdata.iloc[4,7])
        #self.DBfile=self.DBfile.replace("pinkpeak/resources"," ")
        self.DBList.insert("","end",values=self.DBfile)
        self.progressbar.stop()
        self.progressbar.update_idletasks()
               
    def save_DB_as_csv(self):
        self.DBcsvfile=filedialog.asksaveasfile(defaultextension=".csv",filetypes=[("CSV Files", "*.csv")])
        pd.DataFrame(self.DBdata).to_csv(self.DBcsvfile, header=["Precursor","Precursor m/z","Adduct","SMILE","INCHIKEY","Formula","Ion mode","m/z","Intensity"])
        self.DBcsvfile.close()
        print("Saved at "+str(self.DBcsvfile))
        
    def selected_precusor(self,event):
        plt.close()
        self.fragments=[]
        self.fragmentintensity=[]
        self.lowmass=float(self.lowmassentry.get())
        self.highmass=float(self.highmassentry.get())
        self.roundto=float(self.yaxisroundentry.get())
        for item in self.msmspeaks.get_children():
            self.msmspeaks.delete(item)
        item = self.precursors.selection()
        selection=self.precursors.item(item,"values")[0]
#        print(selection)
        self.selectionpre=float(selection)
        self.ID=self.fragmented_precursors.get(self.selectionpre)
        self.ID=float(self.ID[0])
        #print(type(self.ID))
        for spec in self.run:
            if self.ID == spec.ID:
                self.fragments=spec.mz
                self.fragmentintensity=spec.i
        self.msmsplot=mzML_plotter(self.fragments, self.fragmentintensity,
                                   self.selectionpre, self.lowmass, self.highmass,
                                   self.roundto)
        self.msmscanvas = FigureCanvasTkAgg(self.msmsplot, master=self.mzML_viewer)  
        self.msmscanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,columnspan=3,rowspan=5)
        plt.close()
        for i in range (0,len(self.fragments)):
            self.msmspeaks.insert("","end",values=(round(self.fragments[i],5),
                                                   self.fragmentintensity[i]))
        return self.ID
    
    def selected_precusorsb(self,event):
        self.fragments=[]
        self.fragmentintensity=[]
        self.lowmass=float(self.lowmassentrysb.get())
        self.highmass=float(self.highmassentrysb.get())
        self.roundto=float(self.yaxisroundentrysb.get())
        for item in self.msmspeaks.get_children():
            self.msmspeaks.delete(item)
        item = self.precursorssb.selection()
        selection=self.precursorssb.item(item,"values")[0]
#        print(selection)
        self.selectionpre=float(selection)
        self.ID=self.fragmented_precursors.get(self.selectionpre)
        self.ID=float(self.ID[0])
        #print(type(self.ID))
        for spec in self.run:
            if self.ID == spec.ID:
                self.fragments=spec.mz
                self.fragmentintensity=spec.i
        self.msmsplot=mzML_plotter(self.fragments, self.fragmentintensity,
                                   self.selectionpre, self.lowmass, self.highmass,
                                   self.roundto)
        self.msmscanvas = FigureCanvasTkAgg(self.msmsplot, master=self.mzML_viewer)  
        self.msmscanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,columnspan=3,rowspan=5)
        plt.close()
        for i in range (0,len(self.fragments)):
            self.msmspeaks.insert("","end",values=(round(self.fragments[i],5),
                                                   self.fragmentintensity[i]))
        return self.ID
    
    def move_precursor(self,event):
        self.allowedIDssb=[]
        for item in self.sbList.get_children():
            self.sbList.delete(item)
        for item in self.DBprecursors.get_children():
            self.DBprecursors.delete(item)
        item = self.precursorssb.selection()
        selectionsb=self.precursorssb.item(item,"values")[0]
        self.selectionsb=float(selectionsb)
        self.sblist.append(self.selectionsb)
        for i in range (0,len(self.sblist)):
            self.IDsb=self.fragmented_precursors.get(self.sblist[i])
            self.sbList.insert("","end",values=(round(self.sblist[i],5)))
            self.allowedIDssb.append(self.IDsb)
            self.DBprecursors.insert("","end",values=(round(self.sblist[i],5)))
    
    def delete_precursorsb(self,event):
        item = self.sblist.selection()
        selection=self.sblist.item(item,"values")[0]
#        print(selection)
        self.selectiondelpre=float(selection)
        self.ID=self.fragmented_precursors.get(self.selectiondelpre)
        self.ID=float(self.ID[0])
        for ids in self.allowedIDssb:
            if self.ID == ids:
                self.allowedIDssb.remove(ids)
                print("removed "+str(self.selectiondelpre))  
        for item in self.sblist.get_children():
            self.sblist.delete(item)
        for item in self.DBprecursors.get_children():
            self.DBprecursors.delete(item)
        for i in range (0,len(self.allowedIDssb)):
            for spec in self.run:
                if self.allowedIDssb[i] == spec.ID:
                    selected_precursors = spec.selected_precursors
                    for precursor_dict in selected_precursors:
                        self.precursor_mz = precursor_dict["mz"]
                    self.precursors.insert("","end",values=(round(self.precursor_mz,5)))
                    self.DBprecursors.insert("","end",values=(round(self.precursor_mz,5)))
        
    def delete_msn(self,event):
        item = self.precursors.selection()
        selection=self.precursors.item(item,"values")[0]
#        print(selection)
        self.selectiondelmsn=float(selection)
        self.ID=self.fragmented_precursors.get(self.selectiondelmsn)
        self.ID=float(self.ID[0])
        for ids in self.allowedIDs:
            if self.ID == ids:
                self.allowedIDs.remove(ids)
                print("removed "+str(self.selectiondelmsn))  
        for item in self.precursors.get_children():
            self.precursors.delete(item)
        for item in self.DBprecursors.get_children():
            self.DBprecursors.delete(item)
        for i in range (0,len(self.allowedIDs)):
            for spec in self.run:
                if self.allowedIDs[i] == spec.ID:
                    selected_precursors = spec.selected_precursors
                    for precursor_dict in selected_precursors:
                        self.precursor_mz = precursor_dict["mz"]
                    self.precursors.insert("","end",values=(round(self.precursor_mz,5)))
                    self.DBprecursors.insert("","end",values=(round(self.precursor_mz,5)))
                    
    def mzml_masscutoff(self):
        self.lowcutoff=float(self.lowcutoffentry.get())
        self.highcutoff=float(self.highcutoffentry.get())
        self.notallowedIDs=[]
        for item in self.precursors.get_children():
            self.precursors.delete(item)
        for item in self.DBprecursors.get_children():
            self.DBprecursors.delete(item)
        for i in range (0,len(self.allowedIDs)):
            for spec in self.run:
                if self.allowedIDs[i] == spec.ID:
                    selected_precursors = spec.selected_precursors
                    for precursor_dict in selected_precursors:
                        self.precursor_mz = precursor_dict["mz"]
                        if self.precursor_mz >= self.highcutoff or self.precursor_mz <= self.lowcutoff:
                            self.notallowedIDs.append(spec.ID)
                        else:
                            self.precursors.insert("","end",values=(round(self.precursor_mz,5)))
                            self.DBprecursors.insert("","end",values=(round(self.precursor_mz,5)))
        print ("removed "+str(len(self.notallowedIDs))+" peaks")
        for i in self.notallowedIDs:
            self.allowedIDs.remove(i)

    def compare_to_DB(self):
        if self.postcomparisondatasetsorted is not None:
            del self.postcomparisondatasetsorted
        self._DB_mz_list=[]
        self._assignment_string_list=[]
        self._cossim_list=[]
        self._ID_list=[]
        self._input_mzlist=[]
        self.precursorlist=[]
        self.datasetindsexlist=[]
        self.DBprecursormzlist=[]
        self.DBprecursornamelist=[]
        self.databaseindexlist=[]
        self.mzinputlist=[]
        self.intensityinputlist=[]
        self.mzdatabaselist=[]
        self.intensitydatabaselist=[]
        self.similaritylist=[]
        self.adductDBlist=[]
        self.smilesDBlist=[]
        self.INCHIKEYDBlist=[]
        self.FORMULADBlist=[]
        self.ionmodeDBlist=[]
        self.on=self.removetrue.get()
        self.ppmerror=float(self.ppmerrorentry.get())
        print("Checking "+str(len(self.allowedIDs)))
        #self.dataset=pd.DataFrame(data=None, index=self.allowedIDs, columns=["precursor","ID","m/z","i"])
        self.dataset1=[]
        self.dataset2=[]
        self.dataset3=[]
        self.dataset4n=[]
        self.dataset4=[]
        tmpinlist=[]
        for i in range (0,len(self.allowedIDs)):
            for spec in self.run:
                if spec.ID == self.allowedIDs[i]:
                    selected_precursors = spec.selected_precursors
                    for precursor_dict in selected_precursors:
                        self.precursor_mz = precursor_dict["mz"]
                    self.dataset1.append(self.precursor_mz)
                    self.dataset2.append(spec.ID)
                    self.dataset3.append(spec.mz)
                    self.dataset4.append(spec.i)
        for i in self.dataset4: #normalise intensity
            #print(type(i))
            i=list(i)
            #print(type(i))
            scalingfactor=100/max(i)
            tmpinlist=[]
            for j in i:
                normint=j*scalingfactor
                tmpinlist.append(normint)
            self.dataset4n.append(tmpinlist)
        #print((self.dataset4n))
            
        self.dataset=pd.DataFrame(data={"precursor":self.dataset1,
                                        "ID":self.dataset2,
                                        "m/z":self.dataset3,
                                        "i":self.dataset4n})
        self.dataset["precursor"]=pd.to_numeric(self.dataset["precursor"], downcast="float") #downcast from float64
        self.dataset["ID"]=pd.to_numeric(self.dataset["ID"], downcast="integer")
        self.progressbar = ttk.Progressbar(maximum=len(self.dataset1))
        print(self.DBdata.memory_usage())
        #iloc row(ie the precursors),column(ie the data about that precursor)
        for i in range(0,len(self.dataset1)):
            self.progressbar.step(1)
            self.progressbar.update_idletasks()
            print("searching "+str(i+1)+"/"+str(len(self.dataset1)))
            for j in range(0,len(self.DBdata)):
                self.precursor=self.dataset.iloc[i,0]
                self.DBprecursorname=self.DBdata.iloc[j,0]
                self.DBprecursormz=self.DBdata.iloc[j,1]
                self.mzinput=self.dataset.iloc[i,2]
                #print(type(self.mzinput))
                self.intensityinput=self.dataset.iloc[i,3]
                self.mzdatabase=self.DBdata.iloc[j,7]
                #print(type(self.mzdatabase))
                self.intensitydatabase=self.DBdata.iloc[j,8]
                self.adductDB=self.DBdata.iloc[j,2]
                self.smilesDB=self.DBdata.iloc[j,3]
                self.INCHIKEYDB=self.DBdata.iloc[j,4]
                self.FORMULADB=self.DBdata.iloc[j,5]
                self.ionmodeDB=self.DBdata.iloc[j,6]
                self.combinedpeaklist=[]
                for m in range(0,len(self.mzdatabase)):
                    self.combinedpeaklist.append([0,0,0])
                    #print(self.combinedpeaklist)
                    self.combinedpeaklist[m][0]=self.mzdatabase[m]
                    self.combinedpeaklist[m][1]=self.intensitydatabase[m]
                    for n in range(0,len(self.mzinput)):
                        float(self.mzinput[n])
                        float(self.mzdatabase[m])
                        #print(self.mzdatabase[m])
                        ppm=abs(1000000*((self.mzdatabase[m]-self.mzinput[n])/self.mzdatabase[m]))
                        if ppm < self.ppmerror:
                            self.combinedpeaklist[m][2]=self.intensityinput[n]
                #self.cosinesimilarity=cossim(self.combinedpeaklist,self.on)
                self.revsimilarity=revcordinersim(self.mzdatabase,self.mzinput,self.intensitydatabase,self.ppmerror,self.on)
                if self.revsimilarity==0:
                    self.similarity=np.NaN
                elif self.revsimilarity>=-0.2:
                    self.forsimilarity=cordinersim(self.mzdatabase,self.mzinput,self.intensityinput,self.ppmerror,self.on)
                    self.noprecossim=cossim_noprecursor(self.mzdatabase,self.mzinput,self.intensitydatabase,self.intensityinput,self.ppmerror,self.precursor)
                    self.similarity=((3*self.forsimilarity)+(6*self.revsimilarity)+self.noprecossim)/10
                else:
                    self.similarity=np.NaN
                if self.similarity<=0.0:
                    self.similarity=np.NaN
                #if np.isnan(self.cosinesimilarity) != True:
                    #print(self.cosinesimilarity)
                    #print(self.combinedpeaklist)
                self.precursorlist.append(self.precursor)
                self.datasetindsexlist.append(i)
                self.DBprecursormzlist.append(self.DBprecursormz)
                self.DBprecursornamelist.append(self.DBprecursorname)
                self.databaseindexlist.append(j)
                self.mzinputlist.append(self.mzinput)
                self.intensityinputlist.append(self.intensityinput)
                self.mzdatabaselist.append(self.mzdatabase)
                self.intensitydatabaselist.append(self.intensitydatabase)
                self.similaritylist.append(self.similarity)
                self.adductDBlist.append(self.adductDB)
                self.smilesDBlist.append(self.smilesDB)
                self.INCHIKEYDBlist.append(self.INCHIKEYDB)
                self.FORMULADBlist.append(self.FORMULADB)
                self.ionmodeDBlist.append(self.ionmodeDB)
        #del self.dataset
        self.postcomparisondataset=pd.DataFrame(data={"input precursor":self.precursorlist,
                                                "Input index":self.datasetindsexlist,
                                                "Database precursor m/z":self.DBprecursormzlist,
                                                "Database precursor name":self.DBprecursornamelist,
                                                "Database index":self.databaseindexlist,
                                                "input m/z list":self.mzinputlist,
                                                "input i list":self.intensityinputlist,
                                                "DB m/z list":self.mzdatabaselist,
                                                "DB i list":self.intensitydatabaselist,
                                                "Cosine Similarity":self.similaritylist,
                                                "Adduct":self.adductDBlist,
                                                "SMILES":self.smilesDBlist,
                                                "INCHIKEY":self.INCHIKEYDBlist,
                                                "Formula":self.FORMULADBlist,
                                                "Polarity":self.ionmodeDBlist})
        self._DB_mz_list=[]
        self._assignment_string_list=[]
        self._cossim_list=[]
        self._ID_list=[]
        self._input_mzlist=[]
        self.precursorlist=[]
        self.datasetindsexlist=[]
        self.DBprecursormzlist=[]
        self.DBprecursornamelist=[]
        self.databaseindexlist=[]
        self.mzinputlist=[]
        self.intensityinputlist=[]
        self.mzdatabaselist=[]
        self.intensitydatabaselist=[]
        self.similaritylist=[]
        self.adductDBlist=[]
        self.smilesDBlist=[]
        self.INCHIKEYDBlist=[]
        self.FORMULADBlist=[]
        self.ionmodeDBlist=[]
        print(len(self.postcomparisondataset))
        self.postcomparisondataset=self.postcomparisondataset.dropna(subset=['Cosine Similarity'])
        print(len(self.postcomparisondataset))
        self.postcomparisondatasetsorted = self.postcomparisondataset.sort_values(['input precursor',"Cosine Similarity"],ascending=False)
        #print(self.postcomparisondatasetsorted)
        print(self.postcomparisondatasetsorted.info(memory_usage="deep"))
        #del self.postcomparisondataset
        
        
    def compare_to_DB_sb(self):
        if self.postcomparisondatasetsorted is not None:
            del self.postcomparisondatasetsorted
        self.precursorlist=[]
        self.datasetindsexlist=[]
        self.DBprecursormzlist=[]
        self.DBprecursornamelist=[]
        self.databaseindexlist=[]
        self.mzinputlist=[]
        self.intensityinputlist=[]
        self.mzdatabaselist=[]
        self.intensitydatabaselist=[]
        self.similaritylist=[]
        self.adductDBlist=[]
        self.smilesDBlist=[]
        self.INCHIKEYDBlist=[]
        self.FORMULADBlist=[]
        self.ionmodeDBlist=[]
        self.on=self.removetrue.get()
        self.ppmerrorsb=float(self.ppmerrorentrysb.get())
        print("Checking "+str(len(self.allowedIDssb)))
            #self.dataset=pd.DataFrame(data=None, index=self.allowedIDs, columns=["precursor","ID","m/z","i"])
        self.dataset1=[]
        self.dataset2=[]
        self.dataset3=[]
        self.dataset4n=[]
        self.dataset4=[]
        tmpinlist=[]
        for i in range (0,len(self.allowedIDssb)):
            for spec in self.run:
                if spec.ID == self.allowedIDssb[i]:
                    selected_precursors = spec.selected_precursors
                    for precursor_dict in selected_precursors:
                        self.precursor_mz = precursor_dict["mz"]
                    self.dataset1.append(self.precursor_mz)
                    self.dataset2.append(spec.ID)
                    self.dataset3.append(spec.mz)
                    self.dataset4.append(spec.i)
        for i in self.dataset4: #normalise intensity
                #print(type(i))
            i=list(i)
                #print(type(i))
            scalingfactor=100/max(i)
            tmpinlist=[]
            for j in i:
                normint=j*scalingfactor
                tmpinlist.append(normint)
            self.dataset4n.append(tmpinlist)
            #print((self.dataset4n))
                
        self.dataset=pd.DataFrame(data={"precursor":self.dataset1,
                                            "ID":self.dataset2,
                                            "m/z":self.dataset3,
                                            "i":self.dataset4n})
        
        self.progressbarsb = ttk.Progressbar(maximum=len(self.dataset1))
            #iloc row(ie the precursors),column(ie the data about that precursor)
        for i in range(0,len(self.dataset1)):
            self.progressbarsb.step(1)
            self.progressbarsb.update_idletasks()
            print("starting to seach "+str(self.dataset.iloc[i,0]))
            for j in range(0,len(self.DBdata)):
                self.precursor=self.dataset.iloc[i,0]
                self.DBprecursorname=self.DBdata.iloc[j,0]
                self.DBprecursormz=self.DBdata.iloc[j,1]
                self.mzinput=self.dataset.iloc[i,2]
                #print(type(self.mzinput))
                self.intensityinput=self.dataset.iloc[i,3]
                self.mzdatabase=self.DBdata.iloc[j,7]
                #print(type(self.mzdatabase))
                self.intensitydatabase=self.DBdata.iloc[j,8]
                self.adductDB=self.DBdata.iloc[j,2]
                self.smilesDB=self.DBdata.iloc[j,3]
                self.INCHIKEYDB=self.DBdata.iloc[j,4]
                self.FORMULADB=self.DBdata.iloc[j,5]
                self.ionmodeDB=self.DBdata.iloc[j,6]
                self.combinedpeaklist=[]
                for m in range(0,len(self.mzdatabase)):
                    self.combinedpeaklist.append([0,0,0])
                    #print(self.combinedpeaklist)
                    self.combinedpeaklist[m][0]=self.mzdatabase[m]
                    self.combinedpeaklist[m][1]=self.intensitydatabase[m]
                    for n in range(0,len(self.mzinput)):
                        float(self.mzinput[n])
                        #print(self.mzinput[n])
                        
                        float(self.mzdatabase[m])
                        #print(self.mzdatabase[m])
                        ppm=abs(1000000*((self.mzdatabase[m]-self.mzinput[n])/self.mzdatabase[m]))
                        if ppm < self.ppmerror and self.combinedpeaklist[m][2]<self.intensityinput[n]:
                            self.combinedpeaklist[m][2]=self.intensityinput[n]
                self.cosinesimilarity=cossim(self.combinedpeaklist,self.on)
                #if np.isnan(self.cosinesimilarity) != True:
                    #print(self.cosinesimilarity)
                    #print(self.combinedpeaklist)
                print(self.combinedpeaklist)
                self.precursorlist.append(self.precursor)
                self.datasetindsexlist.append(i)
                self.DBprecursormzlist.append(self.DBprecursormz)
                self.DBprecursornamelist.append(self.DBprecursorname)
                self.databaseindexlist.append(j)
                self.mzinputlist.append(self.mzinput)
                self.intensityinputlist.append(self.intensityinput)
                self.mzdatabaselist.append(self.mzdatabase)
                self.intensitydatabaselist.append(self.intensitydatabase)
                self.similaritylist.append(self.cosinesimilarity)
                self.adductDBlist.append(self.adductDB)
                self.smilesDBlist.append(self.smilesDB)
                self.INCHIKEYDBlist.append(self.INCHIKEYDB)
                self.FORMULADBlist.append(self.FORMULADB)
                self.ionmodeDBlist.append(self.ionmodeDB)
        del self.dataset
        self.postcomparisondataset=pd.DataFrame(data={"input precursor":self.precursorlist,
                                                "Input index":self.datasetindsexlist,
                                                "Database precursor m/z":self.DBprecursormzlist,
                                                "Database precursor name":self.DBprecursornamelist,
                                                "Database index":self.databaseindexlist,
                                                "input m/z list":self.mzinputlist,
                                                "input i list":self.intensityinputlist,
                                                "DB m/z list":self.mzdatabaselist,
                                                "DB i list":self.intensitydatabaselist,
                                                "Cosine Similarity":self.similaritylist,
                                                "Adduct":self.adductDBlist,
                                                "SMILES":self.smilesDBlist,
                                                "INCHIKEY":self.INCHIKEYDBlist,
                                                "Formula":self.FORMULADBlist,
                                                "Polarity":self.ionmodeDBlist})
        print(len(self.postcomparisondataset))
        self.postcomparisondataset=self.postcomparisondataset.dropna(subset=['Cosine Similarity'])
        print(len(self.postcomparisondataset))
        self.postcomparisondatasetsorted = self.postcomparisondataset.sort_values(['input precursor',"Cosine Similarity"],ascending=False)
        print(len(self.postcomparisondatasetsorted[9]))
        del self.postcomparisondataset
        print(self.postcomparisondatasetsorted.memory_usage())
        
    def clear_DB_list(self):
        self.names=[]
        self.precursormz=[]
        self.smiles=[]
        self.adducts=[]
        self.inchikey=[]
        self.formula=[]
        self.retentiontime=[]
        self.ccs=[]
        self.ionmode=[]
        self.compoundclass=[]
        self.comments=[]
        self.number_of_peaks=[]
        self.peaksindexlists=[]
        self.mzlist=[]
        self.ilist=[]
        for item in self.DBList.get_children():
            self.DBList.delete(item)
    
    def selected_insert_precusor(self,event):
        self.resultcounter=0
        self.DBinfolabel.destroy()
        self.displaycossim.destroy()
        self.DBinfolabel=tk.Label(self.DB_results,text="Name: \nPrecursor m/z:\nAdduct: \nINCHIKEY: \nFormula: \nIon mode:")
        for item in self.assignments.get_children():
            self.assignments.delete(item)
        item = self.DBprecursors.selection()
        selection=self.DBprecursors.item(item,"values")[0]
        self.selection=float(selection)
        self.selecteddataframe = self.postcomparisondatasetsorted.loc[(round(self.postcomparisondatasetsorted["input precursor"],5) == self.selection)]
        #print(self.selecteddataframe)
        if self.selecteddataframe.empty:
            self.assignments.insert("","end",values=("no matches",0))
        else:
            for i in range(0,((len(self.selecteddataframe)))):
                self.assignments.insert("","end",values=((self.selecteddataframe.iloc[i,3]),
                                                         self.selecteddataframe.iloc[i,9]))
    
    def selected_DB_precusor(self,event):
        plt.close()
        if self.resultcounter!=0:
            self.DBinfolabel.destroy()
        self.DBinfolabel=tk.Label(self.DB_results,text="Name: \nPrecursor m/z:\nAdduct: \nINCHIKEY: \nFormula: \nIon mode:")
        self.lowmassmirror=float(self.lowmassentrymirror.get())
        self.highmassmirror=float(self.highmassentrymirror.get())
        item = self.assignments.selection()
        selection=self.assignments.item(item,"values")[0]
        self.selectionDBname=selection
        self.tempDBnamelist=self.selecteddataframe['Database precursor name'].tolist()
        self.index=self.tempDBnamelist.index(self.selectionDBname)
        cossimselection=self.selecteddataframe.iloc[self.index,9]
        self.cossimselection=float(cossimselection)
        self.mirrorplot=mirror_plotter(self.selecteddataframe.iloc[self.index,5],
                                       self.selecteddataframe.iloc[self.index,6],
                                       self.selecteddataframe.iloc[self.index,7],
                                       self.selecteddataframe.iloc[self.index,8], 
                                       self.selectionDBname, self.lowmassmirror, 
                                       self.highmassmirror)
        self.SMILE=self.selecteddataframe.iloc[self.index,11]
        self.mirrorcanvas = FigureCanvasTkAgg(self.mirrorplot, master=self.DB_results)  
        self.mirrorcanvas.get_tk_widget().grid(row=1,column=1,padx=0,pady=0,
                                               columnspan=3,rowspan=3)
        self.displaycossim=tk.Label(self.DB_results,text=str(self.cossimselection))
        self.displaycossim.grid(row=0, column=2, padx=5,pady=5 )
        self.indicator = self.indictorcanvas.create_oval(20, 10, 10, 20)
        if self.cossimselection <= 0.7:
            self.indictorcanvas.itemconfig(self.indicator, fill="red")
        elif self.cossimselection <= 0.8:
            self.indictorcanvas.itemconfig(self.indicator, fill="orange")
        elif self.cossimselection <= 0.85:
            self.indictorcanvas.itemconfig(self.indicator, fill="yellow")
        elif self.cossimselection <= 0.9:
            self.indictorcanvas.itemconfig(self.indicator, fill="limegreen")
        elif self.cossimselection <= 0.98:
            self.indictorcanvas.itemconfig(self.indicator, fill="green")
        else:
          self.indictorcanvas.itemconfig(self.indicator, fill="darkgreen")
        try:
            m = Chem.MolFromSmiles(str(self.SMILE))
            self.smile=Draw.MolToImage(m)
            ph = ImageTk.PhotoImage(self.smile)
        #self.smilesSVG=Draw.rdMolDraw2D.MolDraw2DSVG(m)
            self.smileslabel = tk.Label(self.DB_results, image=ph)
            self.smileslabel.image=ph
        except:
            self.smileslabel = tk.Label(self.DB_results, text="No smile found for this molecule")
        self.smileslabel.grid(row=2, column=6,padx=5,pady=5,rowspan=1)
        self.infocanvas = tk.Canvas(self.DB_results,bg="#E69",height=50,width=50,highlightthickness=0)
        self.infocanvas.grid(row=3, column=6,padx=5,pady=5,rowspan=1)
        self.DBinfolabel=tk.Label(self.DB_results,text="Name: "+str(self.selectionDBname)+"\nPrecursor m/z: "+str(self.selecteddataframe.iloc[self.index,2])+"\nAdduct: "+str(self.selecteddataframe.iloc[self.index,10])+"\nINCHIKEY: "+str(self.selecteddataframe.iloc[self.index,12])+"\nFormula: "+str(self.selecteddataframe.iloc[self.index,13])+"\nIon mode: "+str(self.selecteddataframe.iloc[self.index,14]))
        self.DBinfolabel.grid(row=3,column=6,padx=5,pady=5)
        self.resultcounter=self.resultcounter+1
    
    def confirm_assignment(self,event):
        item = self.assignments.selection()
        selection=self.assignments.item(item,"values")[0]
        self.selectionDBname=selection
        self.tempDBnamelist=self.selecteddataframe['Database precursor name'].tolist()
        self.index=self.tempDBnamelist.index(self.selectionDBname)
        cossimselection=self.selecteddataframe.iloc[self.index,9]
        self.adductselection=self.selecteddataframe.iloc[self.index,10]
        self.DBprecursorselection=self.selecteddataframe.iloc[self.index,2]
        self.IDselection=self.selecteddataframe.iloc[self.index,1]
        self.cossimselection=float(cossimselection)
        self.assignmentstringselection=self.selectionDBname
        self._DB_mz_list.append(self.DBprecursorselection)
        self._assignment_string_list.append(self.assignmentstringselection)
        self._cossim_list.append(self.cossimselection)
        self._ID_list.append(self.IDselection)
        self._input_mzlist.append(self.selection)
        for item in self.resultsTV.get_children():
            self.resultsTV.delete(item)
        for i in range (0,len(self._DB_mz_list)):
            self.resultsTV.insert("","end",values=((self._DB_mz_list[i]),(self._cossim_list[i]),(self._assignment_string_list[i]),self._input_mzlist[i]))
    
    def unconfirm_assignment(self,event):
        item = self.resultsTV.selection()
        selectiondel=self.resultsTV.item(item,"values")[0]
        self.selectiondel=float(selectiondel)
        self.index=self._DB_mz_list.index(self.selectiondel)
        self._DB_mz_list.remove(self.selectiondel)
        self._assignment_string_list.remove(self._assignment_string_list[self.index])
        self._cossim_list.remove(self._cossim_list[self.index])
        self._ID_list.remove(self._ID_list[self.index])
        self._input_mzlist.remove(self._input_mzlist[self.index])
        for item in self.resultsTV.get_children():
            self.resultsTV.delete(item)
        for i in range (0,len(self._DB_mz_list)):
            self.resultsTV.insert("","end",values=((self._DB_mz_list[i]),(self._cossim_list[i]),(self._assignment_string_list[i]),self._input_mzlist[i]))
        
    def export_results(self):
        self._precursor_mz_list=[]
        self._precursor_int_list=[]
        self.resultsDF=[]
        for i in range (0,len(self._ID_list)):#identified ids should be for whats in the list of results
            for spec in self.run:
                if self._ID_list[i] == spec.ID:
                    selected_precursors = spec.selected_precursors
                    for precursor_dict in selected_precursors:
                        self.precursor_mz = precursor_dict["mz"]
                        self.precursor_int = spec["peak intensity"]
                    self._precursor_mz_list.append(self.selection)
                    self._precursor_int_list.append(self.precursor_int)
        for i in range (0,len(self._DB_mz_list)):
            self.resultsDF.append([self._input_mzlist[i],self._DB_mz_list[i],self._precursor_int_list[i],self._assignment_string_list[i],self._cossim_list[i]])
        self.resultscsvfile=filedialog.asksaveasfile(defaultextension=".csv",filetypes=[("CSV Files", "*.csv")])
        pd.DataFrame(self.resultsDF).to_csv(self.resultscsvfile, header=['Y value','DB precursor m/z','precursor intensity','assignment','cos sim'])
        self.resultscsvfile.close()
        
    def save_mirror_plot(self):
        print("you have clicked to save your mirror plot figure")
        if self.savetypemirror.get() =="SVG":
            self.mirrorsvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
            self.mirrorplot.savefig(self.mirrorsvgfile)
        if self.savetypemirror.get() =="PDF":
            self.mirrorpdffile=filedialog.asksaveasfilename(defaultextension=".pdf",filetypes=[("PDF Files", "*.pdf")])
            self.mirrorplot.savefig(self.mirrorpdffile)
        if self.savetypemirror.get() =="JPG":
            self.mirrorjpgfile=filedialog.asksaveasfilename(defaultextension=".jpg",filetypes=[("Picture Files", "*.jpg")])
            self.mirrorplot.savefig(self.mirrorjpgfile)
        if self.savetypemirror.get() =="PNG":
            self.mirrorpngfile=filedialog.asksaveasfilename(defaultextension=".png",filetypes=[("Picture Files", "*.jpg")])
            self.mirrorplot.savefig(self.mirrorpngfile)
            
    def show_results(self,event):
        plt.close()
        self.lowmassmirror=float(self.lowmassentrymirror.get())
        self.highmassmirror=float(self.highmassentrymirror.get())
        item = self.resultsTV.selection()
        selectionDB=self.resultsTV.item(item,"values")[2]
        self.selectionDBname=selectionDB
        selectionpre=self.resultsTV.item(item,"values")[3]
        self.selectionpre=float(selectionpre)
        self.selecteddataframe = self.postcomparisondatasetsorted.loc[(round(self.postcomparisondatasetsorted["input precursor"],5) == self.selectionpre)]
        #print(len(self.selecteddataframe))
        self.tempDBnamelist=self.selecteddataframe['Database precursor name'].tolist()
        #print(self.tempDBnamelist)
        self.index=self.tempDBnamelist.index(self.selectionDBname)
        self.mirrorplotresults=mirror_plotter(self.selecteddataframe.iloc[self.index,5],
                                       self.selecteddataframe.iloc[self.index,6],
                                       self.selecteddataframe.iloc[self.index,7],
                                       self.selecteddataframe.iloc[self.index,8], 
                                       self.selectionDBname, self.lowmassmirror, 
                                       self.highmassmirror)
        self.mirrorcanvasresults = FigureCanvasTkAgg(self.mirrorplotresults, master=self.results_viewer)  
        self.mirrorcanvasresults.get_tk_widget().grid(row=1,column=1,padx=0,pady=0,
                                               columnspan=3,rowspan=3)
        
def mzML_plotter(xvalues, zvalues,header, lowmass, highmass, roundto=1e8):
        #header=round(header,6)
        y_max=peaksort.roundup(max(zvalues),roundto)
        #x_max=peaksort.roundup(max(xvalues),50)
        fig, ax=plt.subplots(dpi=75)
        ax.vlines(x = xvalues, ymin = 0, ymax = zvalues, colors="black")
        ax.set_xlabel("m/z", size =20,style="italic")
        ax.set_ylabel("Intensity", size=20)
        ax.set_title(str(header), size=20)
        ax.set_ylim([0, y_max])
        ax.set_xlim([lowmass, highmass])
        fig.set_figheight(8)
        fig.set_figwidth(12)
        return fig

def mirror_plotter(xvalues1, zvalues1, xvalues2, zvalues2, header, lowmass, 
                   highmass, roundto=1e8):
    #y_max=peaksort.roundup(max(zvalues1),roundto)
    #x_max=peaksort.roundup(max(xvalues1),50)
    plottinglist=[]
    #plottinglistint=[]
    for i in zvalues2:
        newint=(-i)
        plottinglist.append(newint)
    #y_min=peaksort.roundup(min(plottinglist),roundto)
    fig, ax=plt.subplots(dpi=75)
    ax.vlines(x = xvalues1, ymin = 0, ymax = zvalues1, colors="black")
    ax.vlines(x = xvalues2, ymax = 0, ymin = plottinglist, colors="#FF10F0")
    ax.set_xlabel("m/z", size =20,style="italic")
    ax.set_ylabel(" Relative Intensity %", size=20)
    ax.set_title(str(header), size=20)
    ax.set_ylim([-110, 110])
    ax.set_xlim([lowmass, highmass])
    fig.set_figheight(8)
    fig.set_figwidth(12)
    return fig        
        
def cossim(combinedpeaklist,on=0):
    sum1_sq=0
    sum2_sq=0
    test_sum=0
    non_zero_count=0
    for i in range (0,len(combinedpeaklist)):
        test_sum=test_sum+combinedpeaklist[i][2]
        if combinedpeaklist[i][2]!=0:
            non_zero_count=non_zero_count+1 #counts non-zero intensties
    if test_sum == 0:
        cosine_similarity=np.nan 
    elif non_zero_count==1 and on==1: #if there 
    #is only one non-zero intensity in the combined peaklist for the imported 
    #data and the remove single fragments tick box is on, the cosine sim is set
    #to NaN
        cosine_similarity=np.nan
    else:
        for i in range (0,len(combinedpeaklist)):
            sum1_sq=sum1_sq+(combinedpeaklist[i][1]**2)
            sum2_sq=sum2_sq+(combinedpeaklist[i][2]**2)   
        peaklista=np.array(combinedpeaklist)
        #if peaklista[:,2]:
        #print(peaklista)
        dot_prod=peaklista[:,1]@peaklista[:,2]
        #print(peaklista[:,2])
        #print(dot_prod)
        cosine_similarity=dot_prod/(np.sqrt(sum1_sq)*np.sqrt(sum2_sq))
    return cosine_similarity

def jacsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
    AorB=len(combinedpeaklist)
    jaccard_similarity=AandB/AorB #currently no weighting based on intensity
    return jaccard_similarity

def dicesim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
    AorB=len(combinedpeaklist)
    dice_similarity=(2*AandB)/(AorB+AandB)
    return dice_similarity

def threeWjacsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
    AorB=len(combinedpeaklist)
    threeW_jac_similarity=(3*AandB)/(AorB+AandB+AandB)
    return threeW_jac_similarity

def SSsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    sokal_sneath_similarity=(BnotA+2*AnotB)/(BnotA+AnotB+AandB)
    return sokal_sneath_similarity

def mountsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    mountford_similarity=(2*AandB)/((AandB*(BnotA+AnotB))+(2*BnotA*AnotB))
    return mountford_similarity

def McCsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    mcconnaughey_similarity=((AandB)**2-(AnotB*BnotA))/((AnotB+AandB)*(BnotA+AandB))
    return mcconnaughey_similarity

def DKsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    a=AnotB
    b=BnotA
    c=AandB
    driver_kroeber_similarity=(c*(a+b+2*c))/(2*(a+c)*(b+c))
    return driver_kroeber_similarity

def simpsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    a=AnotB
    b=BnotA
    c=AandB
    simpson_similarity=(c)/(min(a+c),(b+c))
    return simpson_similarity

def BBsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    a=AnotB
    b=BnotA
    c=AandB
    Braun_Banquet_similarity=(c)/(max(a+c),(b+c))
    return Braun_Banquet_similarity

def FMsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    a=AnotB
    b=BnotA
    c=AandB
    fager_McGowan_similarity=((c)/(np.sqrt((a+c)*(b+c))))-((1)/(2*np.sqrt(max(a+c),(b+c))))
    return fager_McGowan_similarity

def hamsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    a=AnotB
    b=BnotA
    c=AandB
    Hamming_similarity=1/(a+b)
    return Hamming_similarity

def hellsim(DBmz,DBint,Acmz,Acint,ppmerror,on=0):
    combinedpeaklist=[]
    AandB=0
    AnotB=0
    BnotA=0
    for m in range(0,len(DBmz)):
        combinedpeaklist.append([0,0,0])
        combinedpeaklist[m][0]=DBmz[m]
        combinedpeaklist[m][1]=DBint[m]
        for n in range(0,len(Acmz)):
            float(Acmz[n])
            float(DBmz[m])
            ppm=abs(1000000*((DBmz[m]-Acmz[n])/DBmz[m]))
            if ppm < ppmerror and combinedpeaklist[m][2]<Acint[n]:
                combinedpeaklist[m][2]=Acint[n]   
    for n in range(0,len(Acmz)):
        ppmlist=[]
        for o in range (0,len(combinedpeaklist)):
            ppm=abs(1000000*((combinedpeaklist[o][0]-Acmz[n])/combinedpeaklist[o][0]))
            if ppm<ppmerror:
                ppmboo=True
            else:
                ppmboo=False
            ppmlist.append(ppmboo)
        if ppmlist.count(True) != 0:
            continue
        else:
            combinedpeaklist.append([0,0,0])
            combinedpeaklist[-1][0]=Acmz[n]
            combinedpeaklist[-1][1]=0
            combinedpeaklist[-1][2]=Acint[n]
    for i in range (0,len(combinedpeaklist)):
        if combinedpeaklist[i].count(0)==0:
            AandB=AandB+1
        elif combinedpeaklist[i][2] == 0:
            AnotB=AnotB+1
        elif combinedpeaklist[i][1] == 0:
            BnotA=BnotA+1
    AorB=len(combinedpeaklist)
    a=AnotB
    b=BnotA
    c=AandB
    Hamming_similarity=1-np.sqrt((1-c)/np.sqrt((a+c)*(b+c)))
    return Hamming_similarity

def cordinersim(DBmz,Acmz,Acint,ppmerror,on=0):
    cordinerpeaklist=[]
    ppmlist=[]
    for n in range(0,len(Acmz)):
        cordinerpeaklist.append([0,0,0,0])
        cordinerpeaklist[-1][0]=Acmz[n]
        cordinerpeaklist[-1][1]=Acint[n]
        for m in DBmz:
            ppm=abs(1000000*((m-Acmz[n])/m))
            if ppm<ppmerror:
                ppmboo=True
                cordinerpeaklist[-1][2]=1
            elif cordinerpeaklist[-1][2]==1:
                break
            else:
                ppmboo=False
                cordinerpeaklist[-1][2]=-1
            ppmlist.append(ppmboo)
            cordinerpeaklist[-1][3]=cordinerpeaklist[-1][2]*cordinerpeaklist[-1][1]
        cordinerpeaklist[-1][1]=Acint[n]
    #print([row[3] for row in cordinerpeaklist])
    fullscore=sum([row[3] for row in cordinerpeaklist])
    #print(fullscore)
    totalint=sum([row[1] for row in cordinerpeaklist])
    cordiner_similarity=fullscore/totalint
    #if on==1:
    #     [row[2] for row in cordinerpeaklist]
    if cordiner_similarity==-1.0:
        cordiner_similarity=np.nan
    return cordiner_similarity

def revcordinersim(DBmz,Acmz,DBint,ppmerror,on=0):
    cordinerpeaklist=[]
    ppmlist=[]
    for n in range(0,len(DBmz)):
        cordinerpeaklist.append([0,0,0,0])
        cordinerpeaklist[-1][0]=DBmz[n]
        cordinerpeaklist[-1][1]=DBint[n]
        for m in Acmz:
            ppm=abs(1000000*((m-DBmz[n])/m))
            if ppm<ppmerror:
                ppmboo=True
                cordinerpeaklist[-1][2]=1
            elif cordinerpeaklist[-1][2]==1:
                break
            else:
                ppmboo=False
                cordinerpeaklist[-1][2]=-1
            ppmlist.append(ppmboo)
            cordinerpeaklist[-1][3]=cordinerpeaklist[-1][2]*cordinerpeaklist[-1][1]
        cordinerpeaklist[-1][1]=DBint[n]
    #print([row[3] for row in cordinerpeaklist])
    fullscore=sum([row[3] for row in cordinerpeaklist])
    #print(fullscore)
    totalint=sum([row[1] for row in cordinerpeaklist])
    revcordiner_similarity=fullscore/totalint
    #if on==1:
    #     [row[2] for row in cordinerpeaklist]
    if revcordiner_similarity==-1.0:
        revcordiner_similarity=np.nan
    return revcordiner_similarity

def cossim_noprecursor(DBmz,Acmz,DBint,Acint,ppmerror,precursor):
    if len(DBmz) or len(Acmz) <= 1:
        cos_sim=0
        return cos_sim
    DBmz1=list(DBmz)
    Acmz1=list(Acmz)
    DBint1=list(DBint)
    Acint1=list(Acint)
    combinedpeaklist=[]
    lowmass=precursor-20
    bad_index_DB=[]
    bad_index_Ac=[]
    for i in range(0,len(DBmz1)):
        if DBmz1[i]>lowmass:
            bad_index_DB.append(i)
    for j in range(0,len(Acmz1)):
        if Acmz1[j]>lowmass:
            bad_index_Ac.append(j)
    for index in sorted(bad_index_DB, reverse=True):
        del DBmz1[index]
        del DBint1[index]
    for index in sorted(bad_index_Ac, reverse=True):
        del Acmz1[index]
        del Acint1[index]
    DBscalingfactor=100/max(DBint)
    tmpinlist=[]
    DBint1n=[]
    for j in DBint1:
        normint=j*DBscalingfactor
        tmpinlist.append(normint)
    DBint1n=tmpinlist
    #print(len(DBmz),len(DBmz1),len(DBint1n))
    Acscalingfactor=100/max(Acint)
    tmpinlist=[]
    Acint1n=[]
    for j in Acint1:
        normint=j*Acscalingfactor
        tmpinlist.append(normint)
    Acint1n=tmpinlist
    for m in range(0,len(DBmz1)):
        combinedpeaklist.append([0,0,0])
        #print(self.combinedpeaklist)
        combinedpeaklist[m][0]=DBmz1[m]
        combinedpeaklist[m][1]=DBint1n[m]
        for n in range(0,len(Acmz1)):
            float(Acmz1[n])
            float(DBmz1[m])
            #print(self.mzdatabase[m])
            ppm=abs(1000000*((DBmz1[m]-Acmz1[n])/DBmz1[m]))
            if ppm < ppmerror:
                combinedpeaklist[m][2]=Acint1n[n]
    cos_sim=cossim(combinedpeaklist,on=0)
    #print(combinedpeaklist)
    return cos_sim

if __name__=="__main__":
    r=read_mzML()
    