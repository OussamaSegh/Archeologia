# -*- coding: utf-8 -*-
""" 
@author: Seghaier Oussama
@Project: 51 
"""
from pyqtgraph.Qt import *
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import pyqtgraph.console
from pyqtgraph.dockarea import *
import numpy as np
import sys
import time
import subprocess
import csv
from math import pi
#from readFile2 import *  # module to read the file written by us 
#from dataTreatements import * # module with all the date treatement functions
############################### Sample Data Test ####################################

##############################Constants####################################################
Units = ['mm3','cm3','dm3','m3','mL','cL','dL','L'] # The Volume Unit Chosen , We Can Add More If necessary
CoeffToMm3 = [1,10**3,10**6,10**9,10**3,10**4,10**5,10**6] # Coeeficients to convert from mm3 to athor units
Unit = 'mm3'
ConesList = []
############################## Samples Vases ##########################3
#Vase 1
def generateVase1():

    Rays1 = [100+i*0.05 for i in range(51)]
    Rays1 += [Rays1[-1]+0.3*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]+0.1*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]-0.1*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]-0.3*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1] for i in range(50)]

    heig1 = [500-i for i in range(501)]
    Vase1 = []
    for i in range(min(len(Rays1),len(heig1))):
        Vase1.append((Rays1[i],heig1[i]))
    return Vase1
#Vase2
def generateVase2():

    Rays1 = [100 for i in range(51)]
    Rays1 += [Rays1[-1]+0.8*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]+0.1*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]-0.1*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]-0.8*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]+0.01*(i+1) for i in range(50)]

    heig1 = [500-i for i in range(501)]
    Vase1 = []
    for i in range(min(len(Rays1),len(heig1))):
        Vase1.append((Rays1[i],heig1[i]))
    return Vase1
#Vase3
def generateVase3():

    Rays1 = [100+i*0.05 for i in range(51)]
    Rays1 += [Rays1[-1]+0.3*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]+0.1*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]-0.1*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1]-0.3*(i+1) for i in range(100)]
    Rays1 += [Rays1[-1] for i in range(50)]

    heig1 = [500-i for i in range(501)]
    Vase1 = []
    for i in range(min(len(Rays1),len(heig1))):
        Vase1.append((Rays1[i],heig1[i]))
    return Vase1
    
    
def readData(datafile):
    table =[]
    with open(datafile,'r') as csvfile :
        csvreader = csv.reader(csvfile) 
        for row in csvreader :
            table.append(row)
    csvfile.close()
    i = len(table)/2
    print(table)
    while table.count([]) > 0:
      
      
        table.remove([])
    #print(table)"""
    for i in range(len(table)):
        table[i][0] = round(float(table[i][0]),3)
        table[i][1] = round(float(table[i][1]),3)
     
    return table
#################################Volume Conversion Function##################
def ConvToMm3(Volume,CurrentUnit):

    coef = CoeffToMm3[Units.index(CurrentUnit)] # we get right the coeffient for the conversion 
    return (round(Volume[0]*coef,2),'mm3')
    
def ConvFromMm3(Volume,NewUnit):

    coef =CoeffToMm3[Units.index(NewUnit)] # we get right the coeffient for the conversion
    return (round(Volume[0]/coef,2),NewUnit)
     
def Convert(Volume,CurrentUnit,NewUnit):

    Vol1 = ConvToMm3(Volume,CurrentUnit)
    Vol2 = ConvFromMm3(Vol1,NewUnit)
    return Vol2
    
############################# Function To Treat The DataFile #############################
def TreatData(FileName):

    listDH = readFile(FileName)
    a = correction(listDH)
    a = distancesToPoints(a[0], a[1])
    return pointsToCircles(a[0], a[1])[0]
    
def CalcVolume(rayList):

    """
    @author = Alban GOUGOUA & Oussama SEGHAIER

    This function figures out the volume in millimeter cube from a list of rays and heights.
    
    rayList = [(R,h).....] list of radius and associated heights that calculates the vase volume
    """
    
    def volume(r1, r2, h) :
        """calculates the volume of a cone truncated with parameters R = radius and H = height. 
        source : https://calculis.net/volume/cone-tronque
        """
        return (h*pi/3)*((r1**2)+(r2**2)+r1*r2)
    
    volum = 0
    #: We figure out the volume by sum of each cone's volume.
    for i in range(len(rayList)) :
        if(i != len(rayList) - 1) :
            h = rayList[i][1] - rayList[i+1][1]
            r1 = rayList[i][0]
            r2 = rayList[i+1][0]
            volum += volume(r1, r2, h)
    
    return (round(volum, 3),'mm3')


def VolVasePortion(ConesList,H):

    """ To Calculate the volume of a portion from the vase """
    if H > ConesList[0][1] or H < 0 :
        return (0,'mm3') 
    i = 0
    while ConesList[i][1] > H and i < len(ConesList):
        i +=1
    return CalcVolume(ConesList[i:])
    
################################Draw a 3D vase  ############################################

def drawWithoutGaps(ConesList,w5):

    """Draw the 3D modelisation of the vase
     ConesList : List = [(R1,H1),...,(Rn,Hn)] model of the vase by circles qnd their heights
     w5        : the window where we will draw the vase
     """
    w5.items = []
    g = gl.GLGridItem()
   
    g.scale(200,200,200)
    g.setSpacing(0.01,0.01,0.01)
    w5.addItem(g)
    w5.update() # the items list set to empty in case we change another vase
    for i in range(len(ConesList)-1) :
        # we draw the Cone number i
        md = gl.MeshData.cylinder(rows=10, cols=100, radius=[ConesList[i+1][0],ConesList[i][0]],length=(ConesList[i][1]-ConesList[i+1][1]))
        colors = np.ones((md.faceCount(), 4), dtype=float)
        colors[::2,0] = 0
        colors[:,1] = np.linspace(0, 1, colors.shape[0])
        md.setFaceColors(colors)
        m5 = gl.GLMeshItem(meshdata=md, smooth=True, shader='viewNormalColor', glOptions='opaque')
        # ze translate the cone number i to the right position
        m5.translate(0,0,ConesList[i+1][1])
        # wz then draw 
        w5.addItem(m5)
    #update the window to show the new vase if there was any before
    w5.update()
    
######################## The Graphical Interface ##########################################

def gui():
    global ConesList

    """ This is the principal function that will be enable the user to interact with the software functionnalities """
    
    app = QtGui.QApplication([]) # the application
    win = QtGui.QMainWindow() # the main window
    HelpWindow = QtGui.QWidget()
    area = DockArea()
    win.setCentralWidget(area)
    win.resize(1024,576)
    win.setWindowTitle('Archeo-Metrologia') # Window Title
    
    ## Create docks, place them into the window one at a time.
    ## The Docks Are the Windows That will contain the necessary information
    """ 
    I will create 6 Docks:
    d1 : this dock contains the necessary buttons like open file,.. and the ability to change the global variables
    d4 : this dock will show the volume details
    d5 : ______________ show the 3D Contents
    d6 : hier we will have a 2d representation of the vase and be able to select a part of the vase 
    """
    d1 = Dock("Parameters", size=(1024,30))     
    d4 = Dock("Volume", size=(244,273))
    d5 = Dock("3D Model",size = (405,273))
    d6 = Dock("2D Model",size = (385,273))
    d1.hideTitleBar()
    d5.hideTitleBar()
    d6.hideTitleBar()
    
 #### positionning docks ######   
    area.addDock(d1, 'top')
    area.addDock(d5, 'bottom',d1)
    area.addDock(d4, 'left',d5)
    area.addDock(d6, 'left',d4) 
 
    w1 = pg.LayoutWidget()
    
    w4 = pg.TableWidget()  
    w4.setWindowTitle("Volume Details")
    w4.setRowCount(2)
    w4.setColumnCount(2) 
    w4.setHorizontalHeaderLabels("Value;Unit;".split(";"))
    w4.setVerticalHeaderLabels("Selected;Total;".split(";"))
    
    w5 = gl.GLViewWidget()
    
    w6 = pg.GraphicsWindow(title="2D Vase")
    w = w6.addPlot()
## Adding widget to docks areas    
    d1.addWidget(w1)
    
    ouvrir = QtGui.QPushButton('Open File')
    
    afficherTrous = QtGui.QCheckBox('Show Gaps')
    
    apropos = QtGui.QPushButton('About US')
    
    Help = QtGui.QPushButton('Help')
    
    Combo = QtGui.QComboBox()
    Combo.setWindowTitle('Choose Volume Unit')
    Combo.addItems(Units)
    
    textbox = QtGui.QTableWidgetItem()
    
    filedialog = pg.FileDialog(None,"Load Vase file..",'C:/', "Vase File (*.csv)")
    
####################### Design #################    
    w1.addWidget(ouvrir, row=0, col=0)
    w1.addWidget(afficherTrous, row=0, col=5)
    w1.addWidget(Help, row=0, col=10)
    w1.addWidget(apropos, row=0, col=15)
    w1.addWidget(Combo,row = 0,col = 2)
    
    #ConesList = generateVase2()    
    def openfile():
        global ConesList
        datafile = 'data'
        FileName = unicode(filedialog.getOpenFileName(None,"Load Vase file..",'C:/', "Vase File (*.csv)"))
        #ConesList1 = TreatData(FileName)[:-1]
        data = open('data','w')
        python3Script = 'python3 read.py ' + FileName
        subprocess.call(['python3','read.py',FileName],stdout = data)        
        data.close()     
        ConesList = readData('data')

        print(ConesList)
        drawWithoutGaps(ConesList,w5)
        w.clear()
        update2D(ConesList)
        
    ouvrir.clicked.connect(openfile)
    
    w5.setBackgroundColor(10,10,10,255)
    w5.show()
    w5.setWindowTitle('Vase 3D Representation')
    w5.setCameraPosition(distance=800)
   
    d5.addWidget(w5)
    
    d4.addWidget(w4)
    d6.addWidget(w6)
    #apropos.clicked.connect(afficher3D)         
    layout = QtGui.QGridLayout()
    HelpWindow.setLayout(layout)
    
    def help():
        ConesList = generateVase2()
        drawWithoutGaps(ConesList,w5)
        w.clear()
        update2D(ConesList)
    
    def apropos():
        ConesList = generateVase3()
        drawWithoutGaps(ConesList,w5)
        w.clear()
        update2D(ConesList)                       
    Help.clicked.connect(lambda: help())                      
    #apropos.clicked.connect(lambda: apropos()) 

    region = pg.LinearRegionItem()
    def update():
        minX, maxX = region.getRegion()
        w.setXRange(minX, maxX, padding=0)

    region.sigRegionChanged.connect(update)

    def updateRegion(window, viewRange):
        rgn = viewRange[0]
        region.setRegion(rgn)

    w.sigRangeChanged.connect(updateRegion)

    region.setRegion([-20, 20])
    # We create the horizontal line
    hLine = pg.InfiniteLine(angle=0, movable=True)
    
    def update2D(ConesList):
        """Function to call every time we import a new file  """
        w.clear() # we clear the last data to draw a new vase
        w.addItem(hLine, ignoreBounds=True) # we add the horizontal line and we draw the vase as explained
        Rays1 = [ConesList[i][0] for i in range(len(ConesList))]
        Rays2 = [-Rays1[i] for i in range(len(Rays1))] # 
        heig1 = [ConesList[i][1] for i in range(len(ConesList))] # heights list
                 
        w.plot(Rays1,heig1)
        w.plot(Rays2,heig1)
        w.plot([Rays1[-1],Rays2[-1]],[heig1[-1],heig1[-1]])
        w.enableAutoScale()
        
    vb = w.vb
          
    def mouseMoved(evt): # every time we move the mose , the line will move with it and recover
        global ConesList 
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if w.sceneBoundingRect().contains(pos):
            mousePoint = vb.mapSceneToView(pos)
            index = int(mousePoint.x())
            hLine.setPos(mousePoint.y())
            Selected = VolVasePortion(ConesList,vb.mapSceneToView(evt[0]).y())
            Total = CalcVolume(ConesList)
            l = [[str(Convert(Selected,Selected[1],Combo.currentText())[0]),str(Combo.currentText())],[str(Convert(Total,Total[1],Combo.currentText())[0]),str(Combo.currentText())]]
            w4.setData(l)
            w4.setHorizontalHeaderLabels("Value;Unit;".split(";"))
            w4.setVerticalHeaderLabels("Selected;Total;".split(";"))
            
    proxy = pg.SignalProxy(w.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved) 
                                  
    win.show()
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
############################################################################################
if __name__ == '__main__':
    gui()    
        
    
   





