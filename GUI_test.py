import Solvers
import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QWidget, QLabel, QSizePolicy, QFrame, QCheckBox, QSlider, QLineEdit, QTextEdit, QGridLayout, QVBoxLayout, QApplication, QPushButton, QComboBox)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
plt.style.use('ggplot')

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
    
    # GOAL # Actualise le graph du canvas
    # IN   # flaot CFL           : Nombre CFL du schema de resolution
    #      # int nbrIter         : Nombre d'iteration souhaité lors de la resolution
    #      # int traceLongueur   : Nombre de courbes composant la trace
    #      # int traceEspacement : Duree d'affichage de la solution dans le passe par la trace (pourcentage [0,100] de la duree totale)
    def UpdatePlot(self, CFL, schemaSolver, nbrIter, traceLongueur, traceEspacement):
        self.ax.cla()
        tailleDomaine = 1000
        if traceLongueur == 1:
            iterationTrace = [nbrIter]
        else :
            iterationTrace = [int(nbrIter * (traceEspacement/100) + nbrIter * ((100-traceEspacement)/100) * (i/(traceLongueur-1))) for i in range(traceLongueur)]
        listSol, pos, duree = Solvers.Get_Multiple_Solution(schemaSolver(CFL, tailleDomaine), iterationTrace)

        for i in range(len(listSol)-1):
            a = i / len(listSol)
            self.ax.plot(pos, listSol[i], 'k-', linewidth = 1.0, alpha = a/2.0)
        self.ax.plot(pos, listSol[len(listSol)-1], 'k-', linewidth = 2, alpha = 1.0)

        self.ax.set_title(r"$ \frac{\partial \Phi}{\partial t} + a \ \frac{\partial \Phi}{\partial x}$  %s  $a = 2 m.s^{-1}$  %s $ x \in [0,1] $  %s %s s" % ("\t","\t","\t",str(duree)), fontsize=20)
        self.draw()

class MainWindow(QWidget):
    
    def __init__(self):
        super().__init__()
        self.initUI()

        # Pour generer le premier graph aux positions par default
        self.UpdateParameters()
        self.UpdateGraph()
    
    # Update CFL, dureeIteration et dureeIteration puis les affiche
    def UpdateParameters(self):
        self.CFL = (1.0,-1.0)[self.signCFL.isChecked()] * self.sliderCFL.value()/100.0
        self.valueCFL.setText(str(self.CFL))

        self.dureeIteration = self.sliderDuree.value()
        self.valueDuree.setText(str(self.dureeIteration))

        self.dureeTrace = int(10 ** (self.sliderTraceDuree.value() / 10.0))
        self.valueTraceDuree.setText(str(self.dureeTrace))

        self.espacementTrace = 100 - self.sliderTraceEspacement.value()
        self.valueTraceEspacement.setText(str(100 - self.espacementTrace))

        solverIndex = self.schemaSelection.currentIndex()
        if solverIndex == 0:
            self.schemaSolver = Solvers.Get_Solver_Ordre1_DecentreArriere
        elif solverIndex == 1:
            self.schemaSolver = Solvers.Get_Solver_Ordre1_DecentreAvant
        elif solverIndex == 2:
            self.schemaSolver = Solvers.Get_Solver_Ordre2_DecentreArriere
        elif solverIndex == 3:
            self.schemaSolver = Solvers.Get_Solver_Ordre2_Centre
        elif solverIndex == 4:
            self.schemaSolver = Solvers.Get_Solver_McCormack
        elif solverIndex == 5:
            self.schemaSolver = Solvers.Get_Solver_LaxFriedrichs
        elif solverIndex == 6:
            self.schemaSolver = Solvers.Get_Solver_BeamWarmimg

    def UpdateGraph(self):
        print("Graph updated")
        self.plotFrame.canvas.UpdatePlot(self.CFL, self.schemaSolver, self.dureeIteration, self.dureeTrace, self.espacementTrace)
    
    def ShowRangeCFL(self):
        self.UpdateParameters()
        Solvers.Show_MaxCFL(self.schemaSolver)

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
        self.schemaSelection.currentIndexChanged.connect(self.UpdateParameters)
        self.schemaSelection.currentIndexChanged.connect(self.UpdateGraph)
        #print(self.schemaSelection.currentIndex())

        self.getCFLbutton = QPushButton("Voir intervalle CFL", self)
        gridParameters.addWidget(self.getCFLbutton, 0,2)
        self.getCFLbutton.clicked.connect(self.ShowRangeCFL)

        #CFL
        self.lineCFL = QLabel("CFL",self)
        gridParameters.addWidget(self.lineCFL, 1,0)

        self.sliderCFL = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderCFL, 1,1)
        self.sliderCFL.setMaximum(300)
        self.sliderCFL.setValue(25)
        self.sliderCFL.valueChanged.connect(self.UpdateParameters)
        self.sliderCFL.sliderReleased.connect(self.UpdateGraph)

        self.valueCFL = QLabel("0",self)
        gridParameters.addWidget(self.valueCFL, 1,2)

        self.signCFL = QCheckBox("Negatif",self)
        gridParameters.addWidget(self.signCFL, 1,3)
        self.signCFL.stateChanged.connect(self.UpdateParameters)
        self.signCFL.stateChanged.connect(self.UpdateGraph)

        #Duree (iteration)
        self.lineDuree = QLabel("Nombre d'iteration",self)
        gridParameters.addWidget(self.lineDuree, 2,0)

        self.sliderDuree = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderDuree, 2,1)
        self.sliderDuree.setMaximum(10000)
        self.sliderDuree.setValue(500)
        self.sliderDuree.valueChanged.connect(self.UpdateParameters)
        self.sliderDuree.sliderReleased.connect(self.UpdateGraph)

        self.valueDuree = QLabel("0",self)
        gridParameters.addWidget(self.valueDuree, 2,2)

        #Trace duree
        self.lineTraceDuree = QLabel("Duree de la trace",self)
        gridParameters.addWidget(self.lineTraceDuree, 3,0)

        self.sliderTraceDuree = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderTraceDuree, 3,1)
        self.sliderTraceDuree.setMaximum(30)
        self.sliderTraceDuree.setMinimum(0)
        self.sliderTraceDuree.setValue(10)
        self.sliderTraceDuree.valueChanged.connect(self.UpdateParameters)
        self.sliderTraceDuree.sliderReleased.connect(self.UpdateGraph)

        self.valueTraceDuree = QLabel("0",self)
        gridParameters.addWidget(self.valueTraceDuree, 3,2)

        #Trace espacement
        self.lineTraceEspacement = QLabel("Espacement de la trace",self)
        gridParameters.addWidget(self.lineTraceEspacement, 4,0)

        self.sliderTraceEspacement = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderTraceEspacement, 4,1)
        self.sliderTraceEspacement.setValue(80)
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
        self.setWindowTitle("Etude de l\'equation d\'advection")    
        self.show()


if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())