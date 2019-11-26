import Solvers
import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QWidget, QLabel, QFrame, QCheckBox, QSlider, QLineEdit, QTextEdit, QGridLayout, QApplication, QPushButton, QComboBox)


class Example(QWidget):
    
    def __init__(self):
        super().__init__()
        
        self.initUI()

    def UpdateGraph(self):
        CSV = (1.0,-1.0)[self.signCSV.isChecked()] * self.sliderCSV.value()/10.0
        self.valueCSV.setText(str(CSV))

        dureeIteration = self.sliderDuree.value()
        self.valueDuree.setText(str(dureeIteration))
        
        
    def button_clicked(self):
        print("clicked")
        print(Solvers.Get_Solution(Solvers.Get_Solver_Ordre1_DecentreArriere(0.5,10), 10))

    def initUI(self):
        
        gridMain = QGridLayout(self)
        
        ############################################
        ## Definition de la fenetre de parametres ##
        ############################################

        self.parametersFrame = QFrame(self)
        self.parametersFrame.setFrameStyle(QFrame.Panel | QFrame.Raised)
        gridParameters = QGridLayout(self.parametersFrame)
        
        self.b1 = QPushButton("Holla quetal",self)
        gridParameters.addWidget(self.b1, 0,2)
        self.b1.clicked.connect(self.button_clicked)

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
        self.sliderCSV.sliderReleased.connect(self.UpdateGraph)

        self.valueCSV = QLabel("0",self)
        gridParameters.addWidget(self.valueCSV, 1,2)

        self.signCSV = QCheckBox("Negatif",self)
        gridParameters.addWidget(self.signCSV, 1,3)
        self.signCSV.stateChanged.connect(self.UpdateGraph)

        #Duree (iteration)
        self.lineDuree = QLabel("Duree (iteration)",self)
        gridParameters.addWidget(self.lineDuree, 2,0)

        self.sliderDuree = QSlider(Qt.Horizontal)
        gridParameters.addWidget(self.sliderDuree, 2,1)
        self.sliderDuree.sliderReleased.connect(self.UpdateGraph)

        self.valueDuree = QLabel("0",self)
        gridParameters.addWidget(self.valueDuree, 2,2)

        ####################################
        ## Definition des la page globale ##
        ####################################
        
        gridMain.addWidget(self.parametersFrame, 2,0)

        self.b2 = QPushButton("test second layout")
        gridMain.addWidget(self.b2, 1,0)


        self.setLayout(gridMain)         
        self.setGeometry(300, 300, 350, 300)
        self.setWindowTitle('Review')    
        self.show()
        

if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())