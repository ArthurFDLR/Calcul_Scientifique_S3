import numpy as np
import matplotlib.pyplot as plt

######################################
## Definitions matrices de calcules ##
######################################

# GOAL # Donne la matrice de calcul de l'equation d'advection par un schema decentre arriere d'ordre 1
# IN   # float CFL : Nombre CFL associe a la methode
#        int N     : Nombre d'intervalles de discretisation
# OUT  # np.array  : Matrice de calcul de taille N*N 
#      # float CFL 
def Get_Solver_Ordre1_DecentreArriere(CFL, N):
    out = np.array(N*[N*[0.0]])
    out[0,N-1] = 1
    
    for i in range(1,N):
        out[i,i-1] = CFL
        out[i,i] = 1.0 - CFL
    
    print(out)
    return out,CFL

# GOAL # Donne la matrice de calcul de l'equation d'advection par un schema decentre avant d'ordre 1
# IN   # float CFL : Nombre CFL associe a la methode
#        int N     : Nombre d'intervalles de discretisation
# OUT  # np.array  : Matrice de calcul de taille N*N 
#      # float CFL 
def Get_Solver_Ordre1_DecentreAvant(CFL, N):
    out = np.array(N*[N*[0.0]])
    out[N-1,0] = 1
    
    for i in range(0,N-1):
        out[i,i+1] = -CFL
        out[i,i] = 1.0 + CFL

    print(out)
    return out,CFL

# GOAL # Donne la matrice de calcul de l'equation d'advection par un schema decentre arriere d'ordre 2
# IN   # float CFL : Nombre CFL associe a la methode
#        int N     : Nombre d'intervalles de discretisation
# OUT  # np.array  : Matrice de calcul de taille N*N 
#      # float CFL 
def Get_Solver_Ordre2_DecentreArriere(CFL, N):
    out = np.array(N*[N*[0.0]])
    out[0,N-2] = 1
    out[1,N-1] = 1
    
    for i in range(2,N):
        out[i,i-2] = (CFL/2.0)
        out[i,i] = 1.0 - (CFL/2.0)

    print(out)
    return out,CFL

# GOAL # Donne la matrice de calcul de l'equation d'advection par un schema centre d'ordre 2
# IN   # float CFL : Nombre CFL associe a la methode
#        int N     : Nombre d'intervalles de discretisation
# OUT  # np.array  : Matrice de calcul de taille N*N 
#      # float CFL 
def Get_Solver_Ordre2_Centre(CFL, N):
    out = np.array(N*[N*[0.0]])
    out[0,N-1] = 1
    out[N-1,0] = 1
    
    for i in range(1,N-1):
        out[i,i-1] = (CFL/2.0)
        out[i,i] = 1.0
        out[i,i+1] = -(CFL/2.0)

    print(out)
    return out,CFL

# GOAL # Donne la matrice de calcul de l'equation d'advection par un schema de Mac Cormack
# IN   # float CFL : Nombre CFL associe a la methode
#        int N     : Nombre d'intervalles de discretisation
# OUT  # np.array  : Matrice de calcul de taille N*N
#      # float CFL 
def Get_Solver_McCormack(CFL, N):
    out = np.array(N*[N*[0.0]])
    out[0,N-1] = (CFL*CFL + CFL)/2.0
    out[N-1,0] = (CFL*CFL - CFL)/2.0
    
    for i in range(1,N):
        out[i,i-1] = (CFL*CFL + CFL)/2.0
    for i in range(0,N):
        out[i,i] = 1.0 - (CFL*CFL)
    for i in range(0,N-1):
        out[i,i+1] = (CFL*CFL - CFL)/2.0

    print(out)
    return out,CFL

# GOAL # Donne la matrice de calcul de l'equation d'advection par un schema de Lax Friedrichs
# IN   # float CFL : Nombre CFL associe a la methode
#        int N     : Nombre d'intervalles de discretisation
# OUT  # np.array  : Matrice de calcul de taille N*N 
#      # float CFL
def Get_Solver_LaxFriedrichs(CFL, N):
    out = np.array(N*[N*[0.0]])
    out[0,N-1] = (1.0 + CFL)/2.0
    out[N-1,0] = (1.0 - CFL)/2.0
    
    for i in range(1,N):
        out[i,i-1] = (1.0 + CFL)/2.0
    for i in range(0,N-1):
        out[i,i+1] = (1.0 - CFL)/2.0

    print(out)
    return out,CFL


########################
## Etude des methodes ##
########################

# GOAL # Donne la solution de l'equation associe au schema donne
# IN   # (np.ndarray(N,N), float) schemaSolver : Schema de la methode sous forme matricielle et nombre CFL associe
#        int dureeT                       : Nombre d'iteration de calcul de la solution
#      # float tailleDomaine                : Taille du domaine en unite SI
#      # float vitesseA                   : Vitesse d'advection de l'equation
# OUT  # np.ndarray(N)                    : Ordonnee (solution)
#        np.ndarray(N)                    : Abscisse (position)
#        float                            : Duree finale
def Get_Solution(schemaSolver, dureeT = 0, tailleDomaine = 1.0, vitesseA = 2.0):
    matriceEDP = schemaSolver[0]
    CFL = schemaSolver[1]
    N = matriceEDP.shape[0]

    deltaX = tailleDomaine / N
    deltaT = deltaX * (CFL / vitesseA)

    position = np.linspace(0.0,tailleDomaine,N)
    phiCurrent = np.array(N*[0.0]) # Condition initiales
    for i in range(int(0.1 * N), int(0.2 * N)): 
        phiCurrent[i] = 1.0

    for i in range(dureeT):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
    
    return phiCurrent, position, dureeT * deltaT

def Show_EvolutionTemporelle(schemaSolver):
    for i in range(10):
        sol, pos, duree = Get_Solution(schemaSolver(0.5, 1000), i*100)
        plt.plot(pos, sol, alpha = i * 0.1)
    plt.show()
## Methodes

def GetSolution_Ordre1_DecentreArriere(CFL, NbrIteration = 500, L = 10, N = 500, paramA = 2.0):
    
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

def GetSolution_Ordre1_DecentreAvant(CFL, NbrIteration = 500, L = 10, N = 500, paramA = 2.0):
    
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

def GetSolution_Ordre2_DecentreGauche(CFL, NbrIteration = 500, L = 10, N = 500, paramA = 2.0):
    
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
        matriceEDP[i,i-2] = (CFL/2.0)
        matriceEDP[i,i] = 1.0 - (CFL/2.0)
    
    print(matriceEDP)
    
    # Initialisation conditions initiales
    for i in range(int(0.1 * N), int(0.2 * N)):
        phiCurrent[i] = 1.0
    phiInit = phiCurrent
    
    # Calcul de la solution
    
    for i in range(NbrIteration):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
    
    return phiCurrent 

def GetSolution_Ordre2_Centre(CFL, NbrIteration = 500, L = 10, N = 500, paramA = 2.0):
    
    deltaX = L/N
    deltaT = deltaX * (CFL / paramA)
    
    print(str(deltaT)[:6] + "*" + str(NbrIteration) + " = " + str(deltaT * NbrIteration)[:6] + " secondes")
    
    # Definition matrice
    phiCurrent = np.array(N*[0.0])
    phiNext = np.array(N*[0.0])
    matriceEDP = np.array(N*[N*[0.0]])
    
    # Initialisation matrice de calcul
    matriceEDP[0,N-1] = 1
    matriceEDP[N-1,0] = 1
    
    for i in range(1,N-1):
        matriceEDP[i,i-1] = (CFL/2.0)
        matriceEDP[i,i] = 1.0
        matriceEDP[i,i+1] = -(CFL/2.0)
    
    print(matriceEDP)
    
    # Initialisation conditions initiales
    for i in range(int(0.1 * N), int(0.2 * N)):
        phiCurrent[i] = 1.0
    phiInit = phiCurrent
    
    # Calcul de la solution
    
    for i in range(NbrIteration):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
    
    return phiCurrent 

def GetSolution_McCormack(CFL, NbrIteration = 500, L = 10, N = 500, paramA = 2.0):
    
    deltaX = L/N
    deltaT = deltaX * (CFL / paramA)
    
    print(str(deltaT)[:6] + "*" + str(NbrIteration) + " = " + str(deltaT * NbrIteration)[:6] + " secondes")
    
    # Definition matrice
    phiCurrent = np.array(N*[0.0])
    phiNext = np.array(N*[0.0])
    matriceEDP = np.array(N*[N*[0.0]])
    
    # Initialisation matrice de calcul
    matriceEDP[0,N-1] = (CFL*CFL + CFL)/2.0
    matriceEDP[N-1,0] = (CFL*CFL - CFL)/2.0
    
    for i in range(1,N):
        matriceEDP[i,i-1] = (CFL*CFL + CFL)/2.0
    for i in range(0,N):
        matriceEDP[i,i] = 1.0 - (CFL*CFL)
    for i in range(0,N-1):
        matriceEDP[i,i+1] = (CFL*CFL - CFL)/2.0
    
    print(matriceEDP)
    
    # Initialisation conditions initiales
    for i in range(int(0.1 * N), int(0.2 * N)):
        phiCurrent[i] = 1.0
    phiInit = phiCurrent
    
    # Calcul de la solution
    
    for i in range(NbrIteration):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
    
    return phiCurrent

def GetSolution_LaxFriedrichs(CFL, NbrIteration = 500, L = 10, N = 500, paramA = 2.0):
    
    deltaX = L/N
    deltaT = deltaX * (CFL / paramA)
    
    print(str(deltaT)[:6] + "*" + str(NbrIteration) + " = " + str(deltaT * NbrIteration)[:6] + " secondes")
    
    # Definition matrice
    phiCurrent = np.array(N*[0.0])
    phiNext = np.array(N*[0.0])
    matriceEDP = np.array(N*[N*[0.0]])
    
    # Initialisation matrice de calcul
    matriceEDP[0,N-1] = (1.0 + CFL)/2.0
    matriceEDP[N-1,0] = (1.0 - CFL)/2.0
    
    for i in range(1,N):
        matriceEDP[i,i-1] = (1.0 + CFL)/2.0
    for i in range(0,N-1):
        matriceEDP[i,i+1] = (1.0 - CFL)/2.0
    
    print(matriceEDP)
    
    # Initialisation conditions initiales
    for i in range(int(0.1 * N), int(0.2 * N)):
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

MODE_Selection = 0

if MODE_Selection == 0:
    Show_EvolutionTemporelle(Get_Solver_Ordre1_DecentreArriere)

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
    CFL = 0.1
    NbrIteration_graph =  [0,100,400,1000]
    #NbrIteration_graph =  [0,1,2,3,4,5,6,7,8,9]
    Etude_Temporelle(GetSolution_Ordre2_DecentreGauche,NbrIteration_graph,CFL)
    
if MODE_Selection == 6: # Etude CFL : ordre 2 decentre arriere
    CFL_graph = [2.0, 1.6, 1.2, 0.8, 0.4]
    Etude_CFL(GetSolution_Ordre2_DecentreGauche,CFL_graph)

if MODE_Selection == 7: # Etude temporelle : ordre 2 centre
    CFL = 0.1
    NbrIteration_graph =  [90,80,70,60,50,40,0]
    Etude_Temporelle(GetSolution_Ordre2_Centre,NbrIteration_graph,CFL)
    
if MODE_Selection == 8: # Etude CFL : ordre 2 centre
    CFL_graph = [0.1]
    Etude_CFL(GetSolution_Ordre2_Centre,CFL_graph)

if MODE_Selection == 9: # Etude temporelle : McCormack
    CFL = 0.5
    NbrIteration_graph =  [1500,800,500,100,0]
    Etude_Temporelle(GetSolution_McCormack,NbrIteration_graph,CFL)

if MODE_Selection == 10: # Etude CFL : McCormack
    CFL_graph = [0.5]
    Etude_CFL(GetSolution_McCormack,CFL_graph)

if MODE_Selection == 11: # Etude temporelle : LaxFriedrichs
    CFL = 0.5
    NbrIteration_graph =  [500,400,300,200,100,50,0]
    Etude_Temporelle(GetSolution_LaxFriedrichs,NbrIteration_graph,CFL)

if MODE_Selection == 12: # Etude CFL : LaxFriedrichs
    CFL_graph = [0.5]
    Etude_CFL(GetSolution_LaxFriedrichs,CFL_graph)