import Solvers
import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QWidget, QLabel, QSizePolicy, QFrame, QCheckBox, QSlider, QLineEdit, QTextEdit, QGridLayout, QVBoxLayout, QApplication, QPushButton, QComboBox)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class WidgetPlot(QWidget):
    def __init__(self, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.setLayout(QVBoxLayout())
        self.canvas = PlotCanvas(self, width=10, height=8)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout().addWidget(self.toolbar)
        self.layout().addWidget(self.canvas)

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        dataInit = [0 for i in range(100)]
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(dataInit, 'r-', linewidth = 0.5)
        #self.ax.set_title("Solution")
        self.draw()

    def UpdatePlot(self, CSV, nbrIter, traceLongueur, traceEspacement):
        self.ax.cla()
        schemaSolver = Solvers.Get_Solver_Ordre1_DecentreArriere
        tailleDomaine = 500
        iterationTrace = [int(nbrIter * (traceEspacement/100) + nbrIter * ((100-traceEspacement)/100) * (i/(traceLongueur-1))) for i in range(traceLongueur)]
        listSol, pos, duree = Solvers.Get_Multiple_Solution(schemaSolver(CSV, tailleDomaine), iterationTrace)

        for i in range(len(listSol)):
            self.ax.plot(pos, listSol[i], 'r-', linewidth = 0.5, alpha = 0.5)

        self.draw()

class MainWindow(QWidget):
    
    def __init__(self):
        super().__init__()
        self.initUI()
    
    # Update CSV, dureeIteration et dureeIteration puis les affiche
    def UpdateParameters(self):
        self.CSV = (1.0,-1.0)[self.signCSV.isChecked()] * self.sliderCSV.value()/100.0
        self.valueCSV.setText(str(self.CSV))

        self.dureeIteration = self.sliderDuree.value()
        self.valueDuree.setText(str(self.dureeIteration))

        self.dureeTrace = self.sliderTraceDuree.value()
        self.valueTraceDuree.setText(str(self.dureeTrace))

        self.espacementTrace = self.sliderTraceEspacement.value()
        self.valueTraceEspacement.setText(str(self.espacementTrace))


    def UpdateGraph(self):
        print("Graph updated")
        self.plotFrame.canvas.UpdatePlot(self.CSV, self.dureeIteration, self.dureeTrace, self.espacementTrace)


    def initUI(self):
        
        ############################################
        ## Definition de la fenetre de parametres ##
        ############################################

        self.parametersFrame = QFrame(self)
        self.parametersFrame.setFrameStyle(QFrame.Panel | QFrame.Raised)
        gridParameters = QGridLayout(self.parametersFrame)

        #Methode
        self.lineMethode = QLabel("Methode",self)
        gridParameters.addWidget(self.lineMethode, 0,0)

        self.schemaSelection = QComboBox(self)
        self.schemaSelection.addItem("1er ordre, decentre a gauche")
        self.schemaSelection.addItem("1er ordre, decentre a droite")
        self.schemaSelection.addItem("2nd ordre, decentre a gauche")
        self.schemaSelection.addItem("2nd ordre, centre")
        self.schemaSelection.addItem("Mac Cormack")
        self.schemaSelection.addItem("Lax-Friecrichs")
        self.schemaSelection.addItem("Warming Beam amont")
        gridParameters.addWidget(self.schemaSelection, 0,1)

        #CSV
        self.lineCSV = QLabel("CSV",self)
        gridParameters.addWidget(self.lineCSV, 1,0)

        self.sliderCSV = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderCSV, 1,1)
        self.sliderCSV.setMaximum(300)
        self.sliderCSV.valueChanged.connect(self.UpdateParameters)
        self.sliderCSV.sliderReleased.connect(self.UpdateGraph)

        self.valueCSV = QLabel("0",self)
        gridParameters.addWidget(self.valueCSV, 1,2)

        self.signCSV = QCheckBox("Negatif",self)
        gridParameters.addWidget(self.signCSV, 1,3)
        self.signCSV.stateChanged.connect(self.UpdateParameters)
        self.signCSV.stateChanged.connect(self.UpdateGraph)

        #Duree (iteration)
        self.lineDuree = QLabel("Nombre d'iteration",self)
        gridParameters.addWidget(self.lineDuree, 2,0)

        self.sliderDuree = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderDuree, 2,1)
        self.sliderDuree.setMaximum(1000)
        self.sliderDuree.valueChanged.connect(self.UpdateParameters)
        self.sliderDuree.sliderReleased.connect(self.UpdateGraph)

        self.valueDuree = QLabel("0",self)
        gridParameters.addWidget(self.valueDuree, 2,2)

        #Trace duree
        self.lineTraceDuree = QLabel("Duree de la trace",self)
        gridParameters.addWidget(self.lineTraceDuree, 3,0)

        self.sliderTraceDuree = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderTraceDuree, 3,1)
        self.sliderTraceDuree.setMaximum(20)
        self.sliderTraceDuree.setMinimum(2)
        self.sliderTraceDuree.valueChanged.connect(self.UpdateParameters)
        self.sliderTraceDuree.sliderReleased.connect(self.UpdateGraph)

        self.valueTraceDuree = QLabel("0",self)
        gridParameters.addWidget(self.valueTraceDuree, 3,2)

        #Trace espacement
        self.lineTraceEspacement = QLabel("Espacement de la trace",self)
        gridParameters.addWidget(self.lineTraceEspacement, 4,0)

        self.sliderTraceEspacement = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderTraceEspacement, 4,1)
        self.sliderTraceEspacement.valueChanged.connect(self.UpdateParameters)
        self.sliderTraceEspacement.sliderReleased.connect(self.UpdateGraph)

        self.valueTraceEspacement = QLabel("0",self)
        gridParameters.addWidget(self.valueTraceEspacement, 4,2)
        

        ####################################
        ## Definition des la page globale ##
        ####################################
        
        layoutMain = QVBoxLayout(self)
        self.plotFrame = WidgetPlot(self)
        layoutMain.addWidget(self.parametersFrame)
        layoutMain.addWidget(self.plotFrame)
        
        self.setLayout(layoutMain)         
        self.setGeometry(300, 300, 350, 300)
        self.setWindowTitle('Review')    
        self.show()
        

if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())