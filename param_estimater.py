# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 16:53:15 2025

@author: lwg513
"""

import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import *
pythonicon = "resources/PINK_logo_transparent.gif"
from PIL import Image, ImageTk

class scan_estimate():
    def __init__(self):
        self.root=tk.Tk()
        self.root.title("Parameter estimator")
        self.root.geometry("450x250")
        
        photo_image = tk.PhotoImage(file=pythonicon)
        print(type(photo_image))
        self.root.wm_iconphoto(False, photo_image)
        
        self.successful_scan=np.nan
        
        self.avilable_d26=[1e-5,7.5e-6,5e-6,4e-6,3e-6, 2.5e-6, 1.1e-6, 1e-6, 7.5e-7, 5e-7, 2.5e-7]
            
        style=ttk.Style(self.root)
        style.configure("TFrame", background="#E69", foreground="white")
        
        nb = ttk.Notebook(self.root)
        self.frame2 = ttk.Frame(nb)
        nb.add(self.frame2, text = "in26")
        self.frame1 = ttk.Frame(nb)
        nb.add(self.frame1, text = "Scans")
        nb.pack(expand=1, fill="both")
        
        self.mz1label=tk.Label(self.frame1,text="m/z of 1st ion")
        self.mz1label.grid(row=0, column=0, padx=5,pady=5)
        self.mz1entry=tk.Entry(self.frame1)
        self.mz1entry.insert(0, "0")
        self.mz1entry.grid(row=1, column=0, padx=5,pady=5)
        
        self.mz2label=tk.Label(self.frame1,text="m/z of 2nd ion")
        self.mz2label.grid(row=3, column=0, padx=5,pady=5)
        self.mz2entry=tk.Entry(self.frame1)
        self.mz2entry.insert(0, "0")
        self.mz2entry.grid(row=4, column=0, padx=5,pady=5)
        
        self.d26label=tk.Label(self.frame1,text="delay time increment")
        self.d26label.grid(row=0, column=1, padx=5,pady=5)
        self.d26entry=tk.Entry(self.frame1)
        self.d26entry.insert(0, "1E-6")
        self.d26entry.grid(row=1, column=1, padx=5,pady=5)
        
        self.maglabel=tk.Label(self.frame1,text="Magnet strength")
        self.maglabel.grid(row=3, column=1, padx=5,pady=5)
        self.magentry=tk.Entry(self.frame1)
        self.magentry.insert(0, "9.4")
        self.magentry.grid(row=4, column=1, padx=5,pady=5)
        
        self.lowprelabel=tk.Label(self.frame2,text="m/z of lowest precursor of interest")
        self.lowprelabel.grid(row=0, column=0, padx=5,pady=5)
        self.lowpre=tk.Entry(self.frame2)
        self.lowpre.insert(0, "0")
        self.lowpre.grid(row=0, column=1, padx=5,pady=5)
        
        self.magstrengthlabel=tk.Label(self.frame2,text="Magnet strength")
        self.magstrengthlabel.grid(row=1, column=0, padx=5,pady=5)
        self.magstrength=tk.Entry(self.frame2)
        self.magstrength.insert(0, "9.4")
        self.magstrength.grid(row=1, column=1, padx=5,pady=5)
        
        self.freqcorlabel=tk.Label(self.frame2, text="Frequency correction\nlow excitation frequency (EXC_Freq_Low)")
        self.freqcorlabel.grid(row=2, column=0, padx=5,pady=5)
        self.freqcor=tk.Entry(self.frame2)
        self.freqcor.insert(0, "120250.0")
        self.freqcor.grid(row=2, column=1, padx=5,pady=5)
        
        self.go_button=tk.Button(self.frame1, text="Go", 
                                   command=self.scan_calculator)
        self.go_button.grid(row=5, column=0, padx=5,pady=5)
        
        self.go_button2=tk.Button(self.frame2, text="Go", 
                                   command=self.d_26_calculator)
        self.go_button2.grid(row=5, column=0, padx=5,pady=5)
        
        self.answerlabel=tk.Label(self.frame1,text="Please fill in m/z")
        self.answerlabel.grid(row=6, column=0, padx=5,pady=5, columnspan=2)
        
        self.answerlabel2=tk.Label(self.frame2,text="Please fill in")
        self.answerlabel2.grid(row=6, column=0, padx=5,pady=5, columnspan=2)
        
        self.root.mainloop()
    
    def scan_calculator(self):
        self.cut_offs=[]
        self.B0=float(self.magstrength.get())
        self.ind=float(self.d26entry.get())
        if self.successful_scan:
            del self.successful_scan
        self.successful_scan=np.nan
        self.scan_numbers=[128,256,512,1024,2048,4096,8192,16348,32696,65392,130784]
        self.mz2=float(self.mz2entry.get())
        self.mz1=float(self.mz1entry.get())
        for i in self.scan_numbers:
            tmp_cut_off=model_sigmacutoff(i,self.mz2,self.B0,self.ind)
            self.cut_offs.append(tmp_cut_off)
        self.min_cut_off=abs(self.mz2-self.mz1)
        for i in range (0,len(self.scan_numbers)):
            if self.min_cut_off>self.cut_offs[i]:
                self.successful_scan=self.scan_numbers[i]
                break
            else:
                continue
        if not self.successful_scan:
            self.successful_scan=np.nan
        print(self.cut_offs)
        #print(self.successful_scan)
        self.answerlabel=tk.Label(self.frame1,text="Assuming equal heights \nI think you'll need "+str(self.successful_scan)+" scans to resolve these peaks with peak picking")
        self.answerlabel.grid(row=6, column=0, padx=5,pady=5,columnspan=2)
    
    def d_26_calculator(self):
        self.B0=float(self.magstrength.get())
        self.lowestprecursor=float(self.lowpre.get())
        self.frequency_correction=float(self.freqcor.get())
        self.cutoffmz=[]
        for i in self.avilable_d26:
            self.cutoffmz.append(in_26_mz(self.B0,i,self.frequency_correction))
        #print(self.cutoffmz)
        for i in range(0,len(self.cutoffmz)):
            if self.cutoffmz[i] < self.lowestprecursor:
                self.answer=self.avilable_d26[i]
                break
        if self.cutoffmz[-1] > self.lowestprecursor:
            self.answer=np.nan
        print(self.answer)
        self.answerlabel2=tk.Label(self.frame2,text=str(self.answer)+" is the maximum d_26 you can use")
        self.answerlabel2.grid(row=6, column=0, padx=5,pady=5, columnspan=2)

def model_sigmacutoff(scans,mz,b0,ind):
    res=(8.4237e6*((scans*b0*ind)/mz))-2.2623
    FWHM=mz/res
    sigma=FWHM/(2*np.sqrt(2*np.log(2)))
    return sigma


def in_26_mz(B0,in_26,freq_corr):
    sf=1/in_26
    Vc=0.5*sf+freq_corr
    #print(Vc)
    mzc=(1.535611E7*B0)/Vc
    #print(mzc)
    return mzc
    
if __name__=="__main__":
    r=scan_estimate()