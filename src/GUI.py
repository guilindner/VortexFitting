import tkinter as tk
from tkinter import *
from PIL import Image, ImageTk

import os

import vortexfitting

class VortexDetection:

    def __init__(self, master):
        self.master = master
        master.title("Vortex Detection")

        self.label1 = Label(master, text="boxsize").grid(row=0, column=0, sticky=W)
        self.int1 = IntVar(master, value=6)
        self.entry1 = Entry(master, width=5, textvariable=self.int1)
        self.entry1.grid(row=0, column=1,sticky=W)
        
        self.label2 = Label(master, text="threshold").grid(row=1, column=0, sticky=W)
        self.double2 = DoubleVar(master, value=0.5)
        self.entry2 = Entry(master, width=5, textvariable=self.double2)
        self.entry2.grid(row=1, column=1,sticky=W)
        
        self.label3 = Label(master, text="Numerical Scheme").grid(row=2, column=0, sticky=W)
        self.var3 = StringVar(master)
        self.schemes = ("Second Order","Least Square","Fourth Order")
        self.var3.set(self.schemes[0])
        self.list3 = OptionMenu(master,self.var3,*self.schemes)
        self.list3.grid(row=2, column=1,sticky=W)
        
        self.label4 = Label(master, text="Detection").grid(row=3, column=0, sticky=W)
        self.var4 = StringVar(master)
        self.methods = ("Q criterion","Delta criterion", "Swirling Strenght")
        self.var4.set(self.methods[2])
        self.list4 = OptionMenu(master,self.var4, *self.methods)
        self.list4.grid(row=3, column=1,sticky=W)
        
        self.button1 = Button(master, text="Run", command=lambda: self.run_detection(self.int1,self.double2,self.var3,self.var4))
        self.button1.grid(columns=1, sticky=W)
        
        self.cadre=Canvas(master,width=800,height=600,bg="white")
        self.cadre.grid(columnspan=2)
        #self.pilImage=Image.open("../results/tk.png")
        #self.im=ImageTk.PhotoImage(self.pilImage)
        #self.cadre.create_image(400,300,image = self.im)
        
    def update_image(self): 
        self.pilImage=Image.open("../results/tk.png")
        self.im=ImageTk.PhotoImage(self.pilImage)
        self.cadre.create_image(400,300,image = self.im)   
        
        
            
    def run_detection(self,int1,double2,var3,var4):
        scheme = 2
        if (var3.get()) == "Second Order":
            scheme = 2
        if (var3.get()) == "Least Square":
            scheme = 22 
        if (var3.get()) == "Fourth Order":
            scheme = 4
        method = "swirling"
        if (var4.get()) == "Q criterion":
            method = "Q"
        if (var4.get()) == "Delta criterion":
            method = "delta" 
        if (var4.get()) == "Swirling Strenght":
            method = "swirling"
        print("python3 vortexfitting.py -b %s -t %s -s %s -d %s" %(int1.get(),double2.get(),scheme,method))
        os.system("python3 vortexfitting.py -b %s -t %s -s %s -d %s" %(int1.get(),double2.get(),scheme,method))
        
        self.update_image()
     
        
root = Tk()
app = VortexDetection(root)
root.mainloop()
