# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 10:55:47 2024

@author: lwg513
"""
import pandas as pd
from tkinter import ttk
import csv
import time
import requests
import numpy as np

class import_DB():
    def get_smiles_from_inchikey(inchikey):
        r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON').json()
        return r['PropertyTable']['Properties'][0]['ConnectivitySMILES']
    def read_MSP(MSP_file,names,precursors,smiles,adducts,inchikey,formula,retentiontime,ccs,ionmode,compoundclass,comments,number_of_peaks,peaksindexlists,mzlist,ilist):
        counter=0
        linecounter=0
        peaksindexlists=[]
        with open(str(MSP_file),"r") as f:
            lines=f.readlines()
            #progressbar = ttk.Progressbar(maximum=len(lines))
            #progressbar.grid(row=5, column=6, padx=5,pady=5)
            for line in lines:
                linecounter=linecounter+1
                if line.find("NAME:") != -1:
                    line=line.strip("NAME: \n")
                    names.append(line)
                    #print("I am importing the infomration for "+str(line))
                    counter=counter+1
                    #print(counter)
                    #progressbar.step(line)
                    continue
                if line.find("PRECURSORMZ:") != -1:
                    line=line.strip("PRECURSORMZ: \n")
                    line=float(line)
                    precursors.append(line)
                    continue
                if line.find("PRECURSORTYPE:") != -1:
                    line=line.strip("PRECURSORTYPE: \n")
                    adducts.append(line)
                    continue
                if line.find("SMILES:") != -1:
                    line=line.strip("SMILES: \n")
                    smiles.append(line)
                    continue
                if line.find("INCHIKEY:") != -1:
                    line=line.strip("INCHIKEY: ")
                    inchikey.append(line)
                    continue
                if line.find("FORMULA:") != -1:
                    line=line.strip("FORMULA: \n")
                    formula.append(line)
                    continue
                if line.find("RETENTIONTIME:") != -1:
                    line=line.strip("RETENTIONTIME: \n")
                    retentiontime.append(line)
                    continue
                if line.find("CCS:") != -1:
                    line=line.strip("CCS: \n")
                    ccs.append(line)
                    continue
                if line.find("IONMODE:") != -1:
                    line=line.strip("IONMODE: \n")
                    ionmode.append(line)
                    continue
                if line.find("COMPOUNDCLASS:") != -1:
                    line=line.strip("COMPOUNDCLASS: \n")
                    compoundclass.append(line)
                    continue
                if line.find("ONTOLOGY:") != -1:
                    line=line.strip("ONTOLOGY: \n")
                    compoundclass.append(line)
                    continue
                if line.find("Comment:") != -1:
                    comments.append(line)
                    continue
                if line.find("Num Peaks:") != -1:
                    nopeaks=line.strip("Num Peaks:\t\n")
                    nopeaks=int(nopeaks)
                    #print(type(nopeaks))
                    number_of_peaks.append(nopeaks)
                    tmp=[]
                    for i in range(1,(nopeaks+1)):
                        lineno=linecounter+i
                        tmp.append(lineno)
                    peaksindexlists.append(tmp)
                    #print(tmp)
                    continue
        for i in range(0,len(peaksindexlists)):
            peaklist=[]
            tmpmzlist=[]
            tmpilist=[]
            ilistn=[]
            
            for j in peaksindexlists[i]:
                index=j-1
                #print(index)
                mzandi=lines[index]
                peaklist.append(mzandi)
            peaklistreader=csv.reader(peaklist, delimiter="\t")
            #print(peaklist)
            for row in peaklistreader:
                (mz, intensity)=row
                mz=float(mz)
                intensity=float(intensity)
                tmpmzlist.append(mz)
                tmpilist.append(intensity)
            #print(tmpilist[0])
            if len(tmpmzlist) == 0:
                tmpmzlist=[1]
                tmpilist=[1]
            mzlist.append(tmpmzlist)
            ilist.append(tmpilist)
            #print(ilist[i])
        for i in ilist: #normalise intensity
            #print(i)
            i=list(i)
            scalingfactor=100/max(i)
            ilistntmp=[]
            for j in i:
                normint=j*scalingfactor
                ilistntmp.append(normint)
            ilistn.append(ilistntmp)
        print(str(len(names))+" names\n"+str(len(precursors))+" precursors\n"+str(len(smiles))+" smiles\n"+str(len(adducts))+" adducts\n"+str(len(inchikey))+" inchikey\n"+str(len(formula))+" formula\n"+str(len(ionmode))+" ionmode\n"+str(len(mzlist))+" mass lists\n"+str(len(ilist))+" intensity lists")
        if len(smiles) != len(names):
            t0=time.time()
            failcounter=0
            smiles=[]
            for i in range (0,len(inchikey)):
                #print(i)
                try:
                    smiles.append(import_DB.get_smiles_from_inchikey(inchikey[i]))
                except:
                    smiles.append("No smile found for this compound")
                    failcounter=failcounter+1
            print(str(failcounter)+" smiles not found")
            t1=time.time()
            print(t1-t0)
        database=pd.DataFrame(data={"Name":names,"Precursor m/z":precursors,"adduct":adducts,"smile":smiles,"INCHIKEY":inchikey,"Formula":formula,"Ion mode":ionmode,"mz list":mzlist,"Intensity list":ilistn})
        #print(type(mzlist[1]))
        print("finished importing")
        return database
    def read_MSP_2(MSP_file,names,precursors,smiles,adducts,inchikey,formula,retentiontime,ccs,ionmode,compoundclass,comments,number_of_peaks,peaksindexlists,mzlist,ilist):
        counter=0
        counter2=0
        linecounter=0
        peaksindexlists=[]
        with open(str(MSP_file),"r") as f:
            lines=f.readlines()
            #progressbar = ttk.Progressbar(maximum=len(lines))
            #progressbar.grid(row=5, column=6, padx=5,pady=5)
            for line in lines:
                linecounter=linecounter+1
                if line.find("NAME:") != -1 or line.find("Name:") != -1:
                    if len(names)>len(precursors):
                        print("check how precursor m/z is saved")
                        return
                    if len(names)>len(adducts):
                        adducts.append("N/A")
                    if len(names)>len(smiles):
                        smiles.append("N/A")
                    if len(names)>len(inchikey):
                        inchikey.append("N/A")
                    if len(names)>len(formula):
                        formula.append("N/A")
                    if len(names)>len(retentiontime):
                        retentiontime.append("N/A")
                    if len(names)>len(ccs):
                        ccs.append("N/A")
                    if len(names)>len(ionmode):
                        ionmode.append("N/A")
                    if len(names)>len(compoundclass):
                        compoundclass.append("N/A")
                    if len(names)>len(comments):
                        comments.append("N/A")
                    if len(peaksindexlists)!=0:
                        #print (len(peaksindexlists))
                        if peaksindexlists[-1]==0:#if number of peaks is zero, remove last entry from each column to not import that data
                            #print(names[-1])
                            del names[-1]
                            del precursors[-1]
                            del adducts[-1]
                            del smiles[-1]
                            del inchikey[-1]
                            del formula[-1]
                            del retentiontime[-1]
                            del ccs[-1]
                            del ionmode[-1]
                            del compoundclass[-1]
                            del comments[-1]
                            del peaksindexlists[-1]
                            del number_of_peaks[-1]
                    if line.find("NAME:") != -1:
                        line=line.strip("NAME: \n")
                    elif line.find("Name:") != -1:
                        line=line.strip("Name: \n")
                    names.append(line)
                    #print("I am importing the infomration for "+str(line))
                    counter=counter+1
                    #print(counter)
                    #progressbar.step(line)
                    continue
                if line.find("PrecursorMz:") != -1:
                    line=line.strip("PrecursorMz: \n")
                    line=float(line)
                    precursors.append(line)
                    continue
                if line.find("PrecursorMZ:") != -1:
                    line=line.strip("PrecursorMZ: \n")
                    line=float(line)
                    precursors.append(line)
                    continue
                if line.find("PRECURSORMZ:") != -1:
                    line=line.strip("PRECURSORMZ: \n")
                    line=float(line)
                    precursors.append(line)
                    continue
                if line.find("Selected Ion m/z:") != -1:
                    line=line.strip("Selected Ion m/z: \n")
                    line=float(line)
                    precursors.append(line)
                    continue
                if line.find("PRECURSORTYPE:") != -1:
                    line=line.strip("PRECURSORTYPE: \n")
                    adducts.append(line)
                    continue
                if line.find("Precursor_type:") != -1:
                    line=line.strip("Precursor_type: \n")
                    adducts.append(line)
                    continue
                if line.find("SMILES:",0,10) != -1:
                    line=line.strip("SMILES: \n")
                    smiles.append(line)
                    continue
                if line.find("INCHIKEY:",0,10) != -1:
                    line=line.strip("INCHIKEY: ")
                    inchikey.append(line)
                    counter2=counter2+1
                    continue
                if line.find("FORMULA:",0,10) != -1:
                    line=line.strip("FORMULA: \n")
                    formula.append(line)
                    continue
                if line.find("RETENTIONTIME:") != -1:
                    line=line.strip("RETENTIONTIME: \n")
                    retentiontime.append(line)
                    continue
                if line.find("CCS:") != -1:
                    line=line.strip("CCS: \n")
                    ccs.append(line)
                    continue
                if line.find("IONMODE:",0,10) != -1:
                    line=line.strip("IONMODE: \n")
                    ionmode.append(line)
                    continue
                if line.find("COMPOUNDCLASS:") != -1:
                    line=line.strip("COMPOUNDCLASS: \n")
                    compoundclass.append(line)
                    continue
                if line.find("Comment:") != -1:
                    comments.append(line)
                    continue
                if line.find("Num Peaks:") != -1:
                    nopeaks=line.strip("Num Peaks: \n")
                    nopeaks=int(nopeaks)
                    #print(type(nopeaks))
                    number_of_peaks.append(nopeaks)
                    tmp=[]
                    for i in range(1,(nopeaks+1)):
                        lineno=linecounter+i
                        tmp.append(lineno)
                    peaksindexlists.append(tmp)
                    #print(tmp)
                    continue
                if line.find("Num peaks:") != -1:
                    nopeaks=line.strip("Num peaks: \n")
                    nopeaks=int(nopeaks)
                    #print(type(nopeaks))
                    number_of_peaks.append(nopeaks)
                    tmp=[]
                    for i in range(1,(nopeaks+1)):
                        lineno=linecounter+i
                        tmp.append(lineno)
                    peaksindexlists.append(tmp)
                    #print(tmp)
                    continue
                if line.find("NUM PEAKS:") != -1:
                    nopeaks=line.strip("NUM PEAKS:\t\n")
                    nopeaks=int(nopeaks)
                    #print(type(nopeaks))
                    number_of_peaks.append(nopeaks)
                    tmp=[]
                    for i in range(1,(nopeaks+1)):
                        lineno=linecounter+i
                        tmp.append(lineno)
                    peaksindexlists.append(tmp)
                    #print(tmp)
                    continue
                if line.find("number of peaks:") != -1:
                    nopeaks=line.strip("number of peaks:\t\n")
                    nopeaks=int(nopeaks)
                    #print(type(nopeaks))
                    number_of_peaks.append(nopeaks)
                    tmp=[]
                    for i in range(1,(nopeaks+1)):
                        lineno=linecounter+i
                        tmp.append(lineno)
                    peaksindexlists.append(tmp)
                    #print(tmp)
                    continue
        for i in range(0,len(peaksindexlists)):
            #print(i)
            peaklist=[]
            tmpmzlist=[]
            tmpilist=[]
            ilistn=[]
            for j in peaksindexlists[i]:
                index=j-1
                #print(index)
                mzandi=lines[index]
                peaklist.append(mzandi)
                if "	" in mzandi:
                    delim="tab"
                    #print("tab")
                elif "," in mzandi:
                    delim="comma"
                    #print("comma")
                else:
                    delim="space"
                    #print("space")
            #print(len(initialmzlist))
            if delim=="tab":
                peaklistreader=csv.reader(peaklist, delimiter="\t")
            if delim=="comma":
                peaklistreader=csv.reader(peaklist, delimiter=",")
            if delim=="space":
                peaklistreader=csv.reader(peaklist, delimiter=" ")
            for row in peaklistreader:
                (mz, intensity)=row
                mz=float(mz)
                intensity=float(intensity)
                tmpmzlist.append(mz)
                tmpilist.append(intensity)
             #print(tmpilist[0])
            if len(tmpmzlist) == 0:
                tmpmzlist=[1]
                tmpilist=[1]
            mzlist.append(tmpmzlist)
            ilist.append(tmpilist)
            
        for i in ilist: #normalise intensity
            #print(i)
            i=list(i)
            scalingfactor=100/max(i)
            ilistntmp=[]
            for j in i:
                normint=j*scalingfactor
                ilistntmp.append(normint)
            ilistn.append(ilistntmp)
        if len(names)>len(precursors):
            precursors.append(0)
        if len(names)>len(adducts):
            adducts.append("N/A")
        if len(names)>len(smiles):
            smiles.append("N/A")
        if len(names)>len(inchikey):
            inchikey.append("N/A")
        if len(names)>len(formula):
            formula.append("N/A")
        if len(names)>len(retentiontime):
            retentiontime.append("N/A")
        if len(names)>len(ccs):
            ccs.append("N/A")
        if len(names)>len(ionmode):
            ionmode.append("N/A")
        if len(names)>len(compoundclass):
            compoundclass.append("N/A")
        if len(names)>len(comments):
            comments.append("N/A")
        print(str(len(names))+" names\n"+str(len(precursors))+" precursors\n"+str(len(smiles))+" smiles\n"+str(len(adducts))+" adducts\n"+str(len(inchikey))+" inchikey\n"+str(len(formula))+" formula\n"+str(len(ionmode))+" ionmode\n"+str(len(mzlist))+" mass lists\n"+str(len(ilist))+" intensity lists")
        if len(smiles) != len(names):
            t0=time.time()
            failcounter=0
            smiles=[]
            for i in range (0,len(inchikey)):
                #print(i)
                try:
                    smiles.append(import_DB.get_smiles_from_inchikey(inchikey[i]))
                except:
                    smiles.append("No smile found for this compound")
                    failcounter=failcounter+1
            print(str(failcounter)+" smiles not found")
            t1=time.time()
            print(t1-t0)
        database=pd.DataFrame(data={"Name":names,"Precursor m/z":precursors,"adduct":adducts,"smile":smiles,"INCHIKEY":inchikey,"Formula":formula,"Ion mode":ionmode,"mz list":mzlist,"Intensity list":ilistn})
        #print(type(mzlist[1]))
        print("finished importing")
        return database
    def read_csv(csv_file,names,precursors,smiles,adducts,inchikey,formula,retentiontime,ccs,ionmode,compoundclass,comments,number_of_peaks,peaksindexlists,mzlist,ilist):
        with open(str(csv_file),"r") as f:
            reader=pd.read_csv(f)
            reader.head()
            precursorname=reader["Precursor"]
            precursormass=reader["Precursor m/z"]
            smilelist=reader["SMILE"]
            adductlist=reader["Adduct"]
            keylist=reader["INCHIKEY"]
            formulalist=reader["Formula"]
            modelist=reader["Ion mode"]
            mzlists=reader["m/z"]
            ilists=reader["Intensity"]
        for i in range(0,len(precursorname)):
            names.append(precursorname[i])
            precursors.append(precursormass[i])
            smiles.append(smilelist[i])
            adducts.append(adductlist[i])
            inchikey.append(keylist[i])
            formula.append(formulalist[i])
            ionmode.append(modelist[i])
            mzlist.append(mzlists[i])
            ilist.append(ilists[i])
        print(str(len(names))+" names\n"+str(len(precursors))+" precursors\n"+str(len(smiles))+" smiles\n"+str(len(adducts))+" adducts\n"+str(len(inchikey))+" inchikey\n"+str(len(formula))+" formula\n"+str(len(ionmode))+" ionmode\n"+str(len(mzlist))+" mass lists\n"+str(len(ilist))+" intensity lists")
        database=pd.DataFrame(data={"Name":names,"Precursor m/z":precursors,"adduct":adducts,"smile":smiles,"INCHIKEY":inchikey,"Formula":formula,"Ion mode":ionmode,"mz list":mzlist,"Intensity list":ilist})
        print("finished importing")
        return database