import numpy as np
import matplotlib.pyplot as plt

## Methodes

def GetSolution_Ordre1_DecentreArriere(CFL, NbrIteration = 500, L = 10, N = 100, paramA = 2.0):
    
    deltaX = L/N
    deltaT = deltaX * (CFL / paramA)
    
    print(str(deltaT)[:6] + "*" + str(NbrIteration) + " = " + str(deltaT * NbrIteration)[:6] + " secondes")
    
    # Definition matrice
    phiCurrent = np.array(N*[0.0])
    phiNext = np.array(N*[0.0])
    matriceEDP = np.array(N*[N*[0.0]])
    
    # Initialisation matrice de calcul
    matriceEDP[0,N-1] = 1
    
    for i in range(1,N):
        matriceEDP[i,i-1] = CFL
        matriceEDP[i,i] = 1.0 - CFL
    
    
    # Initialisation conditions initiales
    for i in range(int(0.1 * N), int(0.2 * N)):
        phiCurrent[i] = 1.0
    phiInit = phiCurrent
    
    # Calcul de la solution
    
    for i in range(NbrIteration):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
    
    return phiCurrent 

def GetSolution_Ordre1_DecentreAvant(CFL, NbrIteration = 500, L = 10, N = 100, paramA = 2.0):
    
    deltaX = L/N
    deltaT = deltaX * (CFL / paramA)
    
    print(str(deltaT)[:6] + "*" + str(NbrIteration) + " = " + str(deltaT * NbrIteration)[:6] + " secondes")
    
    # Definition matrice
    phiCurrent = np.array(N*[0.0])
    phiNext = np.array(N*[0.0])
    matriceEDP = np.array(N*[N*[0.0]])
    
    # Initialisation matrice de calcul
    matriceEDP[N-1,0] = 1
    
    for i in range(0,N-1):
        matriceEDP[i,i+1] = -CFL
        matriceEDP[i,i] = 1.0 + CFL
    
    print(matriceEDP)
    
    # Initialisation conditions initiales
    for i in range(int(0.8 * N), int(0.9 * N)):
        phiCurrent[i] = 1.0
    phiInit = phiCurrent
    
    # Calcul de la solution
    
    for i in range(NbrIteration):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
    
    return phiCurrent 

def GetSolution_Ordre2_DecentreGauche(CFL, NbrIteration = 500, L = 10, N = 100, paramA = 2.0):
    
    deltaX = L/N
    deltaT = deltaX * (CFL / paramA)
    
    print(str(deltaT)[:6] + "*" + str(NbrIteration) + " = " + str(deltaT * NbrIteration)[:6] + " secondes")
    
    # Definition matrice
    phiCurrent = np.array(N*[0.0])
    phiNext = np.array(N*[0.0])
    matriceEDP = np.array(N*[N*[0.0]])
    
    # Initialisation matrice de calcul
    matriceEDP[0,N-2] = 1
    matriceEDP[1,N-1] = 1
    
    for i in range(2,N):
        matriceEDP[i,i-2] = - (CFL/2.0)
        matriceEDP[i,i] = 1.0 + (CFL/2.0)
    
    print(matriceEDP)
    
    # Initialisation conditions initiales
    for i in range(int(0.8 * N), int(0.9 * N)):
        phiCurrent[i] = 1.0
    phiInit = phiCurrent
    
    # Calcul de la solution
    
    for i in range(NbrIteration):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
    
    return phiCurrent 
    
## Etude CFL
    
def Etude_CFL(EDP_Function, CFL_graph, NbrIteration = 70):
    
    for CFL in CFL_graph:
        SolutionCFL = EDP_Function(CFL, NbrIteration)
        plt.plot(SolutionCFL, label = "CFL = " + str(CFL)[:5])
    
    plt.legend()
    plt.grid()
    plt.show()

## Etude temporelle

def Etude_Temporelle(EDP_Function, NbrIteration_graph, CFL):
    
    for Nbr in NbrIteration_graph:
        SolutionCFL = EDP_Function(CFL, Nbr)
        plt.plot(SolutionCFL, label = "Nombre d'itérations : " + str(Nbr))
    
    plt.legend()
    plt.grid()
    plt.show()

## Interface utilisateur

MODE_Selection = 4

if MODE_Selection == 1: # Etude temporelle : ordre 1 decentre arriere
    CFL = 0.2
    NbrIteration_graph = [0,50,100,150,200]
    Etude_Temporelle(GetSolution_Ordre1_DecentreArriere,NbrIteration_graph,CFL)
    
if MODE_Selection == 2: # Etude CFL : ordre 1 decentre arriere
    CFL_graph = [1.0, 0.8, 0.6, 0.4, 0.2]
    Etude_CFL(GetSolution_Ordre1_DecentreArriere,CFL_graph)
    
if MODE_Selection == 3: # Etude temporelle : ordre 1 decentre avant
    CFL = -0.2
    NbrIteration_graph = [0,50,100,150,200]
    Etude_Temporelle(GetSolution_Ordre1_DecentreAvant,NbrIteration_graph,CFL)
    
if MODE_Selection == 4: # Etude CFL : ordre 1 decentre avant
    CFL_graph = [-1.0, -0.8, -0.6, -0.4, -0.2]
    Etude_CFL(GetSolution_Ordre1_DecentreAvant,CFL_graph)

if MODE_Selection == 5: # Etude temporelle : ordre 2 decentre arriere
    CFL = 0.2
    NbrIteration_graph = [0,50,100,150,200]
    Etude_Temporelle(GetSolution_Ordre2_DecentreGauche,NbrIteration_graph,CFL)
    
if MODE_Selection == 6: # Etude CFL : ordre 2 decentre arriere
    CFL_graph = [-2.0]
    Etude_CFL(GetSolution_Ordre2_DecentreGauche,CFL_graph)