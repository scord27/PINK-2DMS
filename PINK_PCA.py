# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 11:22:55 2025

@author: lwg513
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
pythonicon = "resources/PINK_logo_transparent.gif"
from PIL import Image, ImageTk
from sklearn.decomposition import PCA
import time

class whole_2d_PCA():
    def __init__(self):
        self.root=tk.Tk()
        self.root.title("P.I.N.K PCA")
        self.root.geometry("1600x850")
        
        photo_image = tk.PhotoImage(file=pythonicon)
        print(type(photo_image))
        self.root.wm_iconphoto(False, photo_image)
        
        style=ttk.Style(self.root)
        style.configure("TFrame", background="#E69", foreground="white")
        
        nb = ttk.Notebook(self.root)
        self.importtab = ttk.Frame(nb)
        nb.add(self.importtab, text = "Import")
        self.resultstab = ttk.Frame(nb)
        nb.add(self.resultstab, text = "Results")
        
        nb.pack(expand=1, fill="both")
        
        self.contourplot=plt.figure()
        self.contourcanvas = FigureCanvasTkAgg(self.contourplot, master=self.importtab)
        self.contourcanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,rowspan=6)
        self.contourtoolbar = NavigationToolbar2Tk(self.contourcanvas,self.importtab, pack_toolbar=False)
        self.contourtoolbar.update()
        self.contourtoolbar.grid(row=7,column=1,padx=0,pady=0)
        
        self.PCAplot=plt.figure()
        self.twoDPCAcanvas = FigureCanvasTkAgg(self.PCAplot, master=self.resultstab)
        self.twoDPCAcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0)
        self.twoDPCAtoolbar = NavigationToolbar2Tk(self.twoDPCAcanvas,self.resultstab, pack_toolbar=False)
        self.twoDPCAtoolbar.update()
        self.twoDPCAtoolbar.grid(row=2,column=0,padx=0,pady=0)
        self.scree_plot=plt.figure()
        self.screecanvas = FigureCanvasTkAgg(self.scree_plot, master=self.resultstab)  
        self.screecanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0)
        self.threePCAplot=plt.figure()
        self.threeDPCAcanvas = FigureCanvasTkAgg(self.threePCAplot, master=self.resultstab)
        self.threeDPCAcanvas.get_tk_widget().grid(row=1,column=1,padx=0,pady=0)
        self.threeDPCAtoolbar = NavigationToolbar2Tk(self.threeDPCAcanvas,self.resultstab, pack_toolbar=False)
        self.threeDPCAtoolbar.update()
        self.threeDPCAtoolbar.grid(row=2,column=1,padx=0,pady=0)
        self.loadplot=plt.figure()
        self.loadcanvas = FigureCanvasTkAgg(self.loadplot, master=self.resultstab)
        self.loadcanvas.get_tk_widget().grid(row=0,column=2,padx=0,pady=0,rowspan=2)
        self.loadtoolbar = NavigationToolbar2Tk(self.loadcanvas,self.resultstab, pack_toolbar=False)
        self.loadtoolbar.update()
        self.loadtoolbar.grid(row=2,column=2,padx=0,pady=0)
        
        self.scalelabel=tk.Label(self.resultstab,text="max loading:minimum loadings displayed")
        self.scalelabel.grid(row=4,column=3,padx=0,pady=5)
        self.scale=tk.Entry(self.resultstab)
        self.scale.insert(0, "50")
        self.scale.grid(row=5,column=3,padx=0,pady=5)
        self.reload=tk.Button(self.resultstab,text="replot loadings plot",command=self.replotloadings)
        self.reload.grid(row=6,column=3,padx=0,pady=5)
    
        self.savePCA=tk.Button(self.resultstab,text="save PCA plot",command=self.savePCAplot)
        self.savePCA.grid(row=0,column=3,padx=0,pady=0)
        self.savescree=tk.Button(self.resultstab,text="save Scree plot",command=self.savescreeplot)
        self.savescree.grid(row=1,column=3,padx=0,pady=0)
        self.save3dPCA=tk.Button(self.resultstab,text="save 3D PCA plot",command=self.save3dPCAplot)
        self.save3dPCA.grid(row=2,column=3,padx=0,pady=0)
        self.saveload=tk.Button(self.resultstab,text="save loading plot",command=self.saveloadplot)
        self.saveload.grid(row=3,column=3,padx=0,pady=0)
        
        self.importbutton=tk.Button(self.importtab,text="Pick your peaklist",command=self.importpeaklist)
        self.importbutton.grid(row=2, column=0, padx=5,pady=5)
        self.nicknamelabel=tk.Label(self.importtab,text="Name for the imported dataset\ndo before import")
        self.nicknamelabel.grid(row=0, column=0, padx=5,pady=5)
        self.nicknameentry=tk.Entry(self.importtab)
        self.nicknameentry.insert(0, "dataset 1")
        self.nicknameentry.grid(row=1, column=0, padx=5,pady=5)
        
        self.xjitterlabel=tk.Label(self.importtab,text="jitter in x dimension, mDa")
        self.xjitterlabel.grid(row=0, column=2, padx=5,pady=5)
        self.xjitterentry=tk.Entry(self.importtab)
        self.xjitterentry.insert(0,"50")
        self.xjitterentry.grid(row=1, column=2, padx=5,pady=5)
        self.yjitterlabel=tk.Label(self.importtab,text="jitter in y dimension, Da")
        self.yjitterlabel.grid(row=0, column=3, padx=5,pady=5)
        self.yjitterentry=tk.Entry(self.importtab)
        self.yjitterentry.insert(0,"1")
        self.yjitterentry.grid(row=1, column=3, padx=5,pady=5)
        self.masslabel=tk.Label(self.importtab,text="refernce jitter mass")
        self.masslabel.grid(row=2, column=2, padx=5,pady=5)
        self.massentry=tk.Entry(self.importtab)
        self.massentry.insert(0,"300")
        self.massentry.grid(row=2, column=3, padx=5,pady=5)
        self.PCAbutton=tk.Button(self.importtab,text="PCA analysis",command=self.PCA_go)
        self.PCAbutton.grid(row=3, column=2, padx=5,pady=5,columnspan=2)
        
        self.datasettreeview=ttk.Treeview(self.importtab,height=20)
        vsbds = ttk.Scrollbar(self.importtab, orient="vertical", 
                             command=self.datasettreeview.yview)
        vsbds.grid(row=3, column=0, sticky=N+S+W,rowspan=3)
        self.datasettreeview.configure(yscrollcommand=vsbds.set)
        self.datasettreeview["columns"]=("name","location")
        self.datasettreeview.column("#0",width=0)
        self.datasettreeview.column("name",width=150)
        self.datasettreeview.column("location",width=250)
        self.datasettreeview.heading("name",text="dataset name")
        self.datasettreeview.heading("location",text="dataset loc")
        self.datasettreeview.grid(row=3, column=0, padx=5,pady=5,rowspan=3)
        self.datasettreeview.bind("<Double-1>", self.plot_selected_2d) 
        self.datasettreeview.bind("<Double-3>", self.del_selected_2d) 
        
        self.data=[]
        self.nicknamelist=[]
        
        self.root.mainloop()
    
    def dummy(self):
        print("you have run a dummy")
        
    def dummy2(self,event):
        print("you have run a dummy")
        
    def importpeaklist(self):
        filterlist="resources\\filterpeaks.csv"
        with open(str(filterlist),"r") as f:
            reader=pd.read_csv(f)
            reader.head()
            self.bannedx=reader["x value"]
            self.bannedy=reader["y value"]
        self.bannedx=round(self.bannedx,3)
        self.bannedy=round(self.bannedx,0)
        t0=time.time()
        self.nickname=self.nicknameentry.get()
        self.infile=filedialog.askopenfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])
        if not self.infile:
            return
        if self.nickname not in self.nicknamelist:
            self.nicknamelist.append(self.nickname)
        else:
            print("please give unique nickname for dataset")
            self.unique_popup_window()
            return
        self.datasettreeview.insert("","end",values=[self.nickname,self.infile])
        with open(str(self.infile),"r") as f:
            self.peaks=[[],[],[]]
            reader=pd.read_csv(f)
            reader.head()
            self.fragmentmass=reader["X"]
            self.precursormass=reader["Y"]
            self.intensity=reader["Z"]
        self.totalintensity=sum(self.intensity)
        #self.normintesnity=[]
        #for i in range (0,len(self.intensity)):
            #tmp=1000*(self.intensity[i]/self.totalintensity)
        self.normintesnity=1000*(self.intensity/self.totalintensity)
        #self.normintesnity=pd.DataFrame(data=self.normintesnity)
        print("imported "+str(len(self.fragmentmass))+" peaks")
        print(self.bannedx)
        for i in range (0,len(self.fragmentmass)):
            notallowed=0
            for j in range (0,len(self.bannedx)):
                if float(round(self.fragmentmass[i],3)) ==self.bannedx[j]:
                    notallowed=1
                    break
                elif float(round(self.precursormass[i],0)) ==self.bannedy[j]:
                    notallowed=1
                    break
                else: continue
            if notallowed==0:
                self.peaks[0].append(self.fragmentmass[i])
                self.peaks[1].append(self.precursormass[i])
                self.peaks[2].append(self.normintesnity[i])
        print(len(self.peaks[0]))
        #print((self.peaks))
        #print(self.peaks[2])
        #self.peakssorted=sorted(self.peaks,key=lambda x: x[2])
        self.peakssorted=self.peaks
        #print(self.peakssorted[0])
        self.contourplot=TwoDMS_plot(self.peakssorted[0],self.peakssorted[1],self.peakssorted[2])
        self.contourcanvas = FigureCanvasTkAgg(self.contourplot, master=self.importtab)
        self.contourcanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,rowspan=6)
        self.contourtoolbar = NavigationToolbar2Tk(self.contourcanvas,self.importtab, pack_toolbar=False)
        self.contourtoolbar.update()
        self.contourtoolbar.grid(row=7,column=1,padx=0,pady=0)
        tmp=[[self.nickname]]
        for i in range (0,len (self.peakssorted[0])):
            tmp.append([self.peakssorted[0][i],self.peakssorted[1][i],self.peakssorted[2][i]])
        self.data.append(tmp)
        t1=time.time()
        print(str(t1-t0)+" s to import")
        plt.close()
        
    def unique_popup_window(self):
        window = tk.Toplevel()
        label = tk.Label(window, text="please give unique nickname for dataset")
        label.grid(row=0,column=0, padx=5, pady=5)
        button_close = tk.Button(window, text="Okay", command=window.destroy)
        button_close.grid(row=1,column=0, padx=5, pady=5)        
    
    def plot_selected_2d(self,event):
        item = self.datasettreeview.selection()
        selection=self.datasettreeview.item(item,"values")[0]
        self.selection=(selection)
        self.index=self.nicknamelist.index(self.selection)
        self.peaklist=self.data[self.index]
        print(self.peaklist[0])
        #del self.peaklist[0]
        self.peaklist=pd.DataFrame(self.peaklist[1:],columns=["X","Y","Z"])
        #print(self.peaklist)
        self.contourplot=TwoDMS_plot(self.peaklist["X"],self.peaklist["Y"],self.peaklist["Z"])
        self.contourcanvas = FigureCanvasTkAgg(self.contourplot, master=self.importtab)
        self.contourcanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0,rowspan=6)
        self.contourtoolbar = NavigationToolbar2Tk(self.contourcanvas,self.importtab, pack_toolbar=False)
        self.contourtoolbar.update()
        self.contourtoolbar.grid(row=7,column=1,padx=0,pady=0)
        plt.close()
     
    def del_selected_2d(self,event):
        item = self.datasettreeview.selection()
        selection=self.datasettreeview.item(item,"values")[0]
        self.selection=(selection)
        self.index=self.nicknamelist.index(self.selection)
        del self.data[self.index]
        del self.nicknamelist[self.index]
        self.datasettreeview.delete(item)
        #print(self.data)
        
    def PCA_go(self):
        t0=time.time()
        self.scaleloads=float(self.scale.get())
        self.xjitter=float(self.xjitterentry.get())*10**-3
        #print(self.xjitter)
        self.yjitter=float(self.yjitterentry.get())
        self.refmass=float(self.massentry.get())
        self.create3d=create_3d_dataset(self.data, self.nicknamelist, self.xjitter, self.yjitter, self.refmass)
        t1=time.time()
        print(str(t1-t0)+" s to make 3d dataset")
        t0=time.time()
        self.threeDlist=self.create3d[0]
        self.axis0=self.create3d[1]
        self.axis1=self.create3d[2]
        self.axis2=self.create3d[3]
        #print(self.threeDlist)
        self.flatten=flatten_x_and_y(self.threeDlist, self.axis0, self.axis1, self.axis2)
        t1=time.time()
        print(str(t1-t0)+" s to flatten")
        t0=time.time()
        self.flat_array=self.flatten[0]
        print(np.shape(self.flat_array))
        self.axis0=self.flatten[1]
        self.xyaxis=self.flatten[2]
        #print(self.axis0)
        self.PC=PCA_actuallydoingit(self.flat_array,self.axis0)
        self.PCdf=self.PC[0]
        self.explainedvarience=self.PC[1]
        self.singularvalues=self.PC[2]
        self.PClist=list(self.PCdf.columns.values)
        self.loadings=self.PC[3]
        print(self.PCdf)
        self.PCAplot=PCA_2d_plot(self.PCdf.loc[:,'principal component 1'],self.PCdf.loc[:,'principal component 2'],self.axis0,self.explainedvarience)
        self.twoDPCAcanvas = FigureCanvasTkAgg(self.PCAplot, master=self.resultstab)
        self.twoDPCAcanvas.get_tk_widget().grid(row=0,column=0,padx=0,pady=0,rowspan=2)
        self.twoDPCAtoolbar = NavigationToolbar2Tk(self.twoDPCAcanvas,self.resultstab, pack_toolbar=False)
        self.twoDPCAtoolbar.update()
        self.twoDPCAtoolbar.grid(row=2,column=0,padx=0,pady=0)
        self.scree_plot=scree_plot(self.explainedvarience,self.singularvalues,self.PClist)
        self.screecanvas = FigureCanvasTkAgg(self.scree_plot, master=self.resultstab)  
        self.screecanvas.get_tk_widget().grid(row=0,column=1,padx=0,pady=0)
        self.threePCAplot=PCA_3d_plot(self.PCdf.loc[:,'principal component 1'],self.PCdf.loc[:,'principal component 2'],self.PCdf.loc[:,'principal component 3'],self.axis0,self.explainedvarience)
        self.threeDPCAcanvas = FigureCanvasTkAgg(self.threePCAplot, master=self.resultstab)
        self.threeDPCAcanvas.get_tk_widget().grid(row=1,column=1,padx=0,pady=0)
        self.threeDPCAtoolbar = NavigationToolbar2Tk(self.threeDPCAcanvas,self.resultstab, pack_toolbar=False)
        self.threeDPCAtoolbar.update()
        self.threeDPCAtoolbar.grid(row=2,column=1,padx=0,pady=0)
        self.loadplot=loading_plot(self.loadings,self.explainedvarience,self.xyaxis,self.scaleloads)
        self.loadcanvas = FigureCanvasTkAgg(self.loadplot, master=self.resultstab)
        self.loadcanvas.get_tk_widget().grid(row=0,column=2,padx=0,pady=0,rowspan=2)
        self.loadtoolbar = NavigationToolbar2Tk(self.loadcanvas,self.resultstab, pack_toolbar=False)
        self.loadtoolbar.update()
        self.loadtoolbar.grid(row=2,column=2,padx=0,pady=0)
        plt.close()
        t1=time.time()
        print(str(t1-t0)+" s to do PCA analysis")
        
    def savePCAplot(self):
        self.PCAsvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
        self.PCAplot.savefig(self.PCAsvgfile)
            
    def savescreeplot(self):
        self.screesvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
        self.scree_plot.savefig(self.screesvgfile)
        
    def save3dPCAplot(self):
        self.PCA3dsvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
        self.threePCAplot.savefig(self.PCA3dsvgfile)
            
    def saveloadplot(self):
        self.loadsvgfile=filedialog.asksaveasfilename(defaultextension=".svg",filetypes=[("Vector Files", "*.svg")])
        self.loadplot.savefig(self.loadsvgfile)
        
    def replotloadings(self):
        self.scaleloads=float(self.scale.get())
        self.loadplot=loading_plot(self.loadings,self.explainedvarience,self.xyaxis,self.scaleloads)
        self.loadcanvas = FigureCanvasTkAgg(self.loadplot, master=self.resultstab)
        self.loadcanvas.get_tk_widget().grid(row=0,column=2,padx=0,pady=0,rowspan=2)
        self.loadtoolbar = NavigationToolbar2Tk(self.loadcanvas,self.resultstab, pack_toolbar=False)
        self.loadtoolbar.update()
        self.loadtoolbar.grid(row=2,column=2,padx=0,pady=0)
        plt.close()

def mass_jitter(refmass,refjitter,mass,ppm=0):
    constant=(refjitter-ppm*refmass*1e-6)/refmass**2
    PDMT=(constant*mass**2)+(ppm*mass*1e-6)
    return PDMT
        
def TwoDMS_plot(X,Y,Z):
    #print()
    X=list(X)
    Y=list(Y)
    Z=list(Z)
    X.reverse()
    Y.reverse()
    Z.reverse()
    fig, ax=plt.subplots(dpi=75)
    cm=matplotlib.colormaps.get_cmap('RdPu')
    sc = ax.scatter(X, Y, c=np.log(Z), s=0.5, cmap=cm)
    ax.set_ylim([min(Y), max(Y)])
    ax.set_xlim([min(X), max(X)])
    ax.set_xlabel("Fragment m/z", size =20,style="italic")
    ax.set_ylabel("Precursor m/z", size=20, style="italic")
    plt.colorbar(sc, ax=ax, label="log10 intensity", spacing='proportional') #log10 to make lower intensiites more visible
    fig.set_figheight(10)
    fig.set_figwidth(12)
    return fig        

def create_3d_dataset(data,nicknamelist,xjitter,yjitter,refmass):
    unique_y_values=[]
    unique_x_values=[]
    #print(data[0])
    #print(data[0][1])
    #print(data[0][1][0])
    for i in range (0,len(data)):
        x_values=[]
        y_values=[]
        for n in range(1,len(data[i])):
            x_values.append(data[i][n][0])
            y_values.append(data[i][n][1])
        #print(x_values)
        for j in range (0,len(x_values)):
            xmatched=False 
            ymatched=False
            if len(unique_x_values)==0:
                unique_x_values.append(x_values[j])
                unique_y_values.append(y_values[j])
            elif x_values[j] not in unique_x_values:
                for k in range (0,len(unique_x_values)):
                    if abs(unique_x_values[k]-x_values[j]) <= mass_jitter(refmass,xjitter,x_values[j]):
                        xmatched=True
                for m in range (0,len(unique_y_values)):
                    if abs(unique_y_values[m]-y_values[j]) <= mass_jitter(refmass,yjitter,y_values[j]):
                        ymatched=True
                if xmatched == False:
                    unique_x_values.append(x_values[j])
                if ymatched==False:
                    unique_y_values.append(y_values[j])
    #return(unique_x_values,unique_y_values)
    axis0=nicknamelist
    axis1=unique_y_values
    axis2=unique_x_values
    #print(axis2)
    three_d_data=np.zeros((len(axis0),len(axis1),len(axis2)))
    print(np.shape(three_d_data))
    for i in range (0,len(axis0)):
        print("adding data for "+str(axis0[i]))
        for n in range(1,len(data[i])):
            for j in range (0,len(axis1)):
                for k in range (0,len(axis2)):
                    if abs(data[i][n][1]-axis1[j])<= mass_jitter(refmass,yjitter,data[i][n][1]):
                        if abs(data[i][n][0]-axis2[k])<= mass_jitter(refmass,xjitter,data[i][n][0]):
                            if three_d_data[i][j][k] == 0 or data[i][n][2] >= three_d_data[i][j][k]:
                                three_d_data[i][j][k] = data[i][n][2]
    #print(three_d_data[0][:][1])
    return three_d_data, axis0, axis1, axis2

def flatten_x_and_y(array_3d,axis0,axis1,axis2):
    print(axis1)
    flat_array=np.reshape(array_3d,(len(axis0),-1),order="C")
    print(np.shape(flat_array))
    comaxis=[]
    for i in axis1:
        for j in axis2:
            comaxis.append(str(round(i,4))+"\n"+str(round(j,4)))
    print(len(comaxis))
    return flat_array, axis0, comaxis

def PCA_actuallydoingit(ar,labels):
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(ar)
    explained_varience=pca.explained_variance_ratio_
    print(pca.explained_variance_ratio_)
    principalDf = pd.DataFrame(data = principalComponents , columns = ['principal component 1', 'principal component 2', 'principal component 3'])#, 'principal component 4', 'principal component 5', 'principal component 6', 'principal component 7', 'principal component 8', 'principal component 9', 'principal component 10'])
    principalDf.index=labels
    singularvalues=pca.singular_values_
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    print(len(loadings),len(loadings[0]))
    return principalDf,explained_varience,singularvalues,loadings

def PCA_2d_plot(PC1,PC2,Labels,explained_varience):
    fig,ax=plt.subplots()
    ax.scatter(PC1,PC2,color="black")
    for i, txt in enumerate(Labels):
        ax.annotate(txt, (PC1[i], PC2[i]))
    ax.set_xlabel("PC1"+str(round(explained_varience[0],2)*100)+"%")
    ax.set_ylabel("PC2"+str(round(explained_varience[1],2)*100)+"%")
    fig.set_figheight(12)
    fig.set_figwidth(12)
    return fig

def scree_plot(explained_varience,singular_values,PCA_list):
    cumulative_variance=[]
    for i in range(0,len(explained_varience)):
        if i==0:
            cumulative_variance.append(explained_varience[i]*100)
        else:
            cumvar=cumulative_variance[-1]+(explained_varience[i]*100)
            cumulative_variance.append(cumvar)
    fig, ax1 = plt.subplots()
    fig.set_figheight(6)
    fig.set_figwidth(6)
    ax2 = ax1.twinx()
    ax1.set_xlabel('Principal componants', size =20)
    ax1.set_ylabel('Eigenvalue', size =20)
    ax2.set_ylabel('Cumulative varience', size =20)
    ax2.set_ylim(ymin=0,ymax=100)
    ax1.bar(PCA_list,singular_values/100,color="#D09")
    ax2.plot(PCA_list,cumulative_variance,c="black")
    return fig

def PCA_3d_plot(PC1,PC2,PC3,Labels,explained_varience):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(PC1,PC2,PC3,c="red",marker="s")
    ax.set_xlabel("PC1 "+str(round(explained_varience[0],3)*100)+"%")
    ax.set_ylabel("PC2 "+str(round(explained_varience[1],3)*100)+"%")
    ax.set_zlabel("PC3 "+str(round(explained_varience[2],3)*100)+"%")
    for i, txt in enumerate(Labels):
        ax.text(PC1[i], PC2[i], PC3[i],txt)
    fig.set_figheight(6)
    fig.set_figwidth(6)
    return fig

def loading_plot(loadings,explained_varience,labels,scaling):
    #print(loadings[0][0],loadings[0][300])
    #print(loadings[0][0],loadings[300][0])
    list_of_0=[]
    addedloadings=[]
    loadpc1=[]
    loadpc2=[]
    #print(labels)
    for i in range(0,len(loadings)):
        list_of_0.append(0)
        addedloadings.append(np.sqrt((loadings[i][0])**2+abs(loadings[i][1])**2))
        loadpc1.append(round(loadings[i][0],4))
        #print(loadings[i][0])
        loadpc2.append(round(loadings[i][1],4))
    print(len(loadpc1),len(loadpc2),len(labels),len(addedloadings))
    #loadingsdf=pd.DataFrame(data=[loadings[:][0],loadings[:][1],labels,addedloadings],columns=["loadingsx","loadingsy","labels","sum"])
    loadingsdf=pd.DataFrame({"loadingsx":loadpc1,"loadingsy":loadpc2,"labels":labels,"sum":addedloadings})
    loadingsdf=loadingsdf.sort_values("sum")
    print("sorted loadings")
    fig, ax = plt.subplots()
    ax.set_xlabel("PC1 "+str(round(explained_varience[0],3)*100)+"%")
    ax.set_ylabel("PC2 "+str(round(explained_varience[1],3)*100)+"%")
    for i in range (0,len(loadingsdf.loc[:,"loadingsx"])):
        if abs(loadingsdf.loc[i,"loadingsx"])>(max(loadpc1)/scaling) or abs(loadingsdf.loc[i,"loadingsy"])>(max(loadpc2)/scaling):
            #print(loadingsdfs.iloc[i,0])
            xlist=[0,loadingsdf.loc[i,"loadingsx"]]
            ylist=[0,loadingsdf.loc[i,"loadingsy"]]
            ax.plot(xlist,ylist)
            if abs(loadingsdf.loc[i,"loadingsx"])>(5*max(loadpc1)/scaling) or abs(loadingsdf.loc[i,"loadingsy"])>(5*max(loadpc2)/scaling):
                txt=loadingsdf.loc[i,"labels"]
                ax.annotate(txt, (round(loadingsdf.loc[i,"loadingsx"],4), round(loadingsdf.loc[i,"loadingsy"],4)))
            #print (loadingsdf.iloc[i,0],loadings[i][0])
    return fig

    
if __name__=="__main__":
    r=whole_2d_PCA()
    