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

    return out,CFL


########################
## Etude des methodes ##
########################

# GOAL # Donne la solution de l'equation associe au schema donne
# IN   # (np.ndarray(N,N), float) solver  : Schema de la methode sous forme matricielle et nombre CFL associe (cf. Get_Solver_...)
#        int dureeT                       : Nombre d'iteration de calcul de la solution
#      # float tailleDomaine              : Taille du domaine en unite SI
#      # float vitesseA                   : Vitesse d'advection de l'equation
# OUT  # np.ndarray(N)                    : Ordonnee (solution)
#        np.ndarray(N)                    : Abscisse (position)
#        float                            : Duree finale
def Get_Solution(solver, dureeT = 0, tailleDomaine = 1.0, vitesseA = 2.0):
    matriceEDP = solver[0]
    CFL = solver[1]
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

# GOAL # Donne les solution de l'equation associe au schema donne pour differents nombres d'iterations
# IN   # (np.ndarray(N,N), float) solver        : Schema de la methode sous forme matricielle et nombre CFL associe (cf. Get_Solver_...)
#      # list(int)                              : Iteration auquelles on souhaite le resultat
#      # float tailleDomaine                    : Taille du domaine en unite SI
#      # float vitesseA                         : Vitesse d'advection de l'equation
# OUT  # list(np.ndarray(N))                    : Abscisses (position)
#      # np.ndarray(N)                          : Ordonnee (solution)
#      # float                                  : Duree finale
def Get_Multiple_Solution(solver, iterationStop, tailleDomaine = 1.0, vitesseA = 2.0):
    matriceEDP = solver[0]
    CFL = solver[1]
    N = matriceEDP.shape[0]

    deltaX = tailleDomaine / N
    deltaT = deltaX * (CFL / vitesseA)

    position = np.linspace(0.0,tailleDomaine,N)

    phiCurrent = np.array(N*[0.0]) # Condition initiales
    for i in range(int(0.1 * N), int(0.2 * N)): 
        phiCurrent[i] = 1.0
    
    out = []

    limitIteration = max(iterationStop)
    for i in range(limitIteration+1):
        phiCurrent = np.dot(matriceEDP,phiCurrent)
        if i in iterationStop:
            out.append(phiCurrent)
    return out, position, limitIteration * deltaT

'''
schemaSolver = Get_Solver_Ordre1_DecentreArriere
listSol, pos, duree = Get_Multiple_Solution(schemaSolver(0.5, 1000), [50,100,150])
print(listSol)
'''