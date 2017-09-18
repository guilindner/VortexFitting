import tkinter as tk
from tkinter import *
from PIL import Image, ImageTk
from tkinter.filedialog import askopenfilename

import os

import vortexfitting



class VortexDetection:

    def __init__(self, master):
        
        self.master = master
        master.title("Vortex Detection")

        self.button0 = Button(master, text="Choose file", command=lambda: self.OpenFile())
        self.button0.grid(columns=1, sticky=W)
        self.name = "../data/test_dataHIT.nc"
        self.label0 = Label(master, text=self.name).grid(row=0, column=1, sticky=W)
        #self.var0 = StringVar(master, value="../data/test_dataHIT.nc")
        #self.entry0 = Entry(master, width=40, textvariable=self.var0)
        #self.entry0.grid(row=0, column=1,sticky=W)

        self.label1 = Label(master, text="boxsize").grid(row=1, column=0, sticky=W)
        self.int1 = IntVar(master, value=6)
        self.entry1 = Entry(master, width=5, textvariable=self.int1)
        self.entry1.grid(row=1, column=1,sticky=W)
        
        self.label2 = Label(master, text="threshold").grid(row=2, column=0, sticky=W)
        self.double2 = DoubleVar(master, value=0.5)
        self.entry2 = Entry(master, width=5, textvariable=self.double2)
        self.entry2.grid(row=2, column=1,sticky=W)
        
        self.label3 = Label(master, text="Numerical Scheme").grid(row=3, column=0, sticky=W)
        self.var3 = StringVar(master)
        self.schemes = ("Second Order","Least Square","Fourth Order")
        self.var3.set(self.schemes[0])
        self.list3 = OptionMenu(master,self.var3,*self.schemes)
        self.list3.grid(row=3, column=1,sticky=W)
        
        self.label4 = Label(master, text="Detection").grid(row=4, column=0, sticky=W)
        self.var4 = StringVar(master)
        self.methods = ("Q criterion","Delta criterion", "Swirling Strength")
        self.var4.set(self.methods[2])
        self.list4 = OptionMenu(master,self.var4, *self.methods)
        self.list4.grid(row=4, column=1,sticky=W)
        
        self.button1 = Button(master, text="Run", command=lambda: self.run_detection())
        self.button1.grid(columns=1, sticky=W)
        
        self.cadre=Canvas(master,width=800,height=600)#,bg="white")
        self.cadre.grid(columnspan=2)
        self.pilImage=Image.open("../results/tk.png")
        print(self.pilImage.mode)
        self.im=ImageTk.PhotoImage(self.pilImage)
        self.cadre.create_image(400,300,image = self.im)
        
    def update_image(self): 
        self.pilImage=Image.open("../results/tk.png")
        self.im=ImageTk.PhotoImage(self.pilImage)
        self.cadre.create_image(400,300,image = self.im)   
        
    def OpenFile(self):
        self.name = askopenfilename(initialdir="../data/",
                               filetypes =(("NetCDF4 file", "*.nc"),("All Files","*.*")),
                               title = "Choose a file.")
        self.label0 = Label(self.master, text=self.name).grid(row=0, column=1, sticky=W)
        print (self.name)        
            
    def run_detection(self):
        scheme = 2
        if (self.var3.get()) == "Second Order":
            scheme = 2
        if (self.var3.get()) == "Least Square":
            scheme = 22 
        if (self.var3.get()) == "Fourth Order":
            scheme = 4
        method = "swirling"
        if (self.var4.get()) == "Q criterion":
            method = "Q"
        if (self.var4.get()) == "Delta criterion":
            method = "delta" 
        if (self.var4.get()) == "Swirling Strength":
            method = "swirling"
        print("python3 vortexfitting.py -i %s -b %s -t %s -s %s -d %s" %(self.name,self.int1.get(),self.double2.get(),scheme,method))
        os.system("python3 vortexfitting.py -i %s -b %s -t %s -s %s -d %s" %(self.name,self.int1.get(),self.double2.get(),scheme,method))
        
        self.update_image()
     
        
root = Tk()
app = VortexDetection(root)
root.mainloop()
