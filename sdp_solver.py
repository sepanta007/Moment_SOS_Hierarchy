import cvxpy as cp
import numpy as np


def generation_SDP_texte(c, filename):
    """
    Génère un fichier texte contenant la formulation du SDP.
    Paramètres :
    c : liste de contraintes sous forme de chaînes de caractères
    filename : nom du fichier de sortie (par défaut "sdp_problem.txt")
    """
    with open(filename, "w") as f:
        f.write("Maximize lambda\n")
        f.write("Subject to:\n")
        f.write("  Q0 >> 0\n")
        f.write("  Q1 >> 0\n")
        f.write("  Q2 >> 0\n")

        for l in c:
            for constraint in l:
                f.write(f"  {constraint}\n")
    print(f"Le fichier {filename} a été généré.")


# dico_monomials représente B0(x) (voir le rapport)
# help1 correspond a B1(x) fois x1
# help2 correspond a B1(x) fois x2
# help3 correspond a B1(x) fois x3
# help4 correspond a B1(x)
# d0 taille de Q0
# d1 taille de Q1
# d2 taille de Q2
def SDP(dico_monomials, help1, help2, help3, help4, d0, d1, d2, filename):
    # Définir les variables pour les matrices et lambda
    Q0 = cp.Variable((d0, d0), symmetric=True)
    Q1 = cp.Variable((d1, d1), symmetric=True)
    Q2 = cp.Variable((d2, d2), symmetric=True)
    lambda_var = cp.Variable()

    # C va correspondre à la somme des matrices F + Q(1)2 + Q(1)1 + ... + Q1 - Q2 qui sera égale a Q0 pour avoir ensuite les contraintes linéaires
    # C_str est la representation en chaine de caractères : on garde une trace pour avoir sur papier les contraintes
    C = np.zeros((d0, d0), dtype = object)
    C_str = [["" for _ in range(d0)] for _ in range(d0)]

    # Simule les deplacement de ligne et colonne et les rajout de zero (déjà fait dans l'initialisation de C), pour avoir les matrices voulu

    # Constructions de Q(1)2 + Q(1)2
    for i in range(d1):
        for j in range(d1):
                C[dico_monomials[help1[i]],dico_monomials[help1[j]]] += Q2[i,j] -1 * Q1[i,j]

                if C_str[dico_monomials[help1[i]]][dico_monomials[help1[j]]] == '':
                    C_str[dico_monomials[help1[i]]][dico_monomials[help1[j]]] += 'Q2['+str(i + 1)+','+str(j + 1)+'] - Q1['+str(i + 1)+','+str(j + 1)+']'
                else:
                    C_str[dico_monomials[help1[i]]][dico_monomials[help1[j]]] += ' + Q2['+str(i + 1)+','+str(j + 1)+'] - Q1['+str(i + 1)+','+str(j + 1)+']'

    # Construction de Q(1)2 - Q(1)2
    for i in range(d1):
        for j in range(d1):
                C[dico_monomials[help2[i]],dico_monomials[help2[j]]] += Q2[i,j] -1 * Q1[i,j]

                if C_str[dico_monomials[help2[i]]][dico_monomials[help2[j]]] == '' :
                    C_str[dico_monomials[help2[i]]][dico_monomials[help2[j]]] += 'Q2['+str(i + 1)+','+str(j + 1)+'] - Q1['+str(i + 1)+','+str(j + 1)+']'
                else:
                    C_str[dico_monomials[help2[i]]][dico_monomials[help2[j]]] += ' + Q2['+str(i + 1)+','+str(j + 1)+'] - Q1['+str(i + 1)+','+str(j + 1)+']'

    # Construction de Q(2)2 - Q(2)2 et additionne avec les matrices déjà construites
    for i in range(d1):
        for j in range(d1):
                C[dico_monomials[help3[i]],dico_monomials[help3[j]]] += Q2[i,j] -1 * Q1[i,j]

                if C_str[dico_monomials[help3[i]]][dico_monomials[help3[j]]] == '':
                    C_str[dico_monomials[help3[i]]][dico_monomials[help3[j]]] += 'Q2['+str(i + 1)+','+str(j + 1)+'] - Q1['+str(i + 1)+','+str(j + 1)+']'
                else:
                    C_str[dico_monomials[help3[i]]][dico_monomials[help3[j]]] += '+ Q2['+str(i + 1)+','+str(j + 1)+'] - Q1['+str(i + 1)+','+str(j + 1)+']'

    # Construction de Q(3)2 - Q(3)2 et additionne avec les matrices déjà construites
    for i in range(d1):
        for j in range(d1):
                C[dico_monomials[help4[i]],dico_monomials[help4[j]]] += -1 * Q2[i,j] + Q1[i,j]
                C_str[dico_monomials[help4[i]]][dico_monomials[help4[j]]] += '- Q2['+str(i + 1)+','+str(j + 1)+'] + Q1['+str(i + 1)+','+str(j + 1)+']'




    # Construction de F et additionne avec les matrices déjà construites
    C[0, 0] += -1 * lambda_var
    C_str[0][0] += ' - lambda'

    C[10, 10] += 1
    C_str[10][10] += ' + 1'

    C[10, 15] += -1
    C_str[10][15] += ' - 1'

    C[13, 10] += -1
    C_str[13][10] += ' - 1'

    C[14, 14] += 3
    C_str[14][14] += ' + 3'


    C[16, 11] += -1
    C_str[16][11] += ' - 1'

    C[17, 19] += -1
    C_str[17][19] += ' - 1'

    C[18, 16] += -1
    C_str[18][16] += ' - 1'

    C[19, 12] += -1
    C_str[19][12] += ' - 1'

    C[19, 19] += 1
    C_str[19][19] += ' + 1'

    C[16, 16] += 1
    C_str[16][16] += ' + 1'

    # Contrainte de positivité définie
    constraints = [
        Q0 >> 0,  # Q0 est semi-définie positive
        Q1 >> 0,  # Q1 est semi-définie positive
        Q2 >> 0,  # Q2 est semi-définie positive
    ]


    # Egalité des sommes des matrices reconstruites et Q0
    for i in range(d0):
        for j in range(d0):
            constraints += [C[i,j] == Q0[i,j]]

            if C_str[i][j] == '':
                C_str[i][j] = 'Q0['+str(i + 1)+','+str(j + 1)+'] = 0'
            else:
                C_str[i][j] = 'Q0['+str(i + 1)+','+str(j + 1)+'] = ' + C_str[i][j]


    # Définir l'objectif : maximiser lambda
    objective = cp.Maximize(lambda_var)

    # Créer et résoudre le problème
    problem = cp.Problem(objective, constraints)
    problem.solve()

    # Afficher les résultats
    print(f"Valeur optimale de lambda : {lambda_var.value}")
    print(f"Q0 = {Q0.value}")
    print(f"Q1 = {Q1.value}")
    print(f"Q2 = {Q2.value}")
    print(f"Statut du problème : {problem.status}")

    generation_SDP_texte(C_str, filename)



# ------------------RELAXATION (SOS)3------------------
print('RELAXATION (SOS)3')
# Représente B0(x) dans le rapport
dico_monomials = {'1':0,
                  'x1':1,'x2':2,'x3':3,
                  'x1**2':4,'x2**2':5,'x3**2':6,'x1*x2':7,'x1*x3':8,'x2*x3':9,
                  'x1**3':10,'(x1**2)*x2':11,'(x1**2)*x3':12,'x1*(x2**2)':13,'x1*x2*x3':14,'x1*(x3**2)':15,'x2**3':16,'(x2**2)*x3':17,'x2*(x3**2)':18,'x3**3':19}

# Dimensions des matrices carrés Q0, Q1 et Q2
d0, d1, d2 = 20, 10, 10
# help1 correspond a B1(x) fois x1
# help2 correspond a B1(x) fois x2
# help3 correspond a B1(x) fois x3
# help4 correspond a B1(x)
help1 = ['x1','x1**2','x1*x2','x1*x3','x1**3','x1*(x2**2)','x1*(x3**2)','(x1**2)*x2','(x1**2)*x3','x1*x2*x3']
help2 = ['x2','x1*x2','x2**2','x2*x3','(x1**2)*x2','x2**3','x2*(x3**2)','x1*(x2**2)','x1*x2*x3','(x2**2)*x3']
help3 = ['x3','x1*x3','x2*x3','x3**2','(x1**2)*x3','(x2**2)*x3','x3**3','x1*x2*x3','x1*(x3**2)','x2*(x3**2)']
help4 = ['1','x1','x2','x3','x1**2','x2**2','x3**2','x1*x2','x1*x3','x2*x3']

SDP(dico_monomials,help1,help2,help3,help4,d0,d1,d2,"sdp_problem_relaxation_t=3.txt")



# ------------------RELAXATION (SOS)4------------------
print('\n\nRELAXATION (SOS)4')

# Représente B0(x) dans le rapport
dico_monomials = {
    '1': 0,
    'x1': 1, 'x2': 2, 'x3': 3,
    'x1**2': 4, 'x2**2': 5, 'x3**2': 6, 'x1*x2': 7, 'x1*x3': 8, 'x2*x3': 9,
    'x1**3': 10, '(x1**2)*x2': 11, '(x1**2)*x3': 12, 'x1*(x2**2)': 13, 'x1*x2*x3': 14, 'x1*(x3**2)': 15,
    'x2**3': 16, '(x2**2)*x3': 17, 'x2*(x3**2)': 18, 'x3**3': 19,

    'x1**4': 20, 'x2**4': 21, 'x3**4': 22,
    '(x1**3)*x2': 23, '(x1**3)*x3': 24, 'x1*(x2**3)': 25, '(x2**3)*x3': 26, 'x1*(x3**3)': 27, 'x2*(x3**3)': 28,
    '(x1**2)*(x2**2)': 29, '(x1**2)*(x3**2)': 30, '(x2**2)*(x3**2)': 31,
    '(x1**2)*x2*x3': 32, 'x1*(x2**2)*x3': 33, 'x1*x2*(x3**2)': 34
}

# Dimensions des matrices carres Q0, Q1 et Q2
d0, d1, d2 = 35, 20, 20
# help1 correspond a B1(x) fois x1
# help2 correspond a B1(x) fois x2
# help3 correspond a B1(x) fois x3
# help4 correspond a B1(x)
help1 = [
    'x1', 'x1**2', 'x1*x2', 'x1*x3', 'x1**3', 'x1*(x2**2)', 'x1*(x3**2)',
    '(x1**2)*x2', '(x1**2)*x3', 'x1*x2*x3', 'x1**4', '(x1**3)*x2', '(x1**3)*x3',
    '(x1**2)*(x2**2)', '(x1**2)*x2*x3', '(x1**2)*(x3**2)', 'x1*(x2**3)', 'x1*(x2**2)*x3',
    'x1*x2*(x3**2)', 'x1*(x3**3)'
]

help2 = [
    'x2', 'x1*x2', 'x2**2', 'x2*x3', '(x1**2)*x2', 'x2**3', 'x2*(x3**2)', 'x1*(x2**2)', 'x1*x2*x3', '(x2**2)*x3',
    '(x1**3)*x2', '(x1**2)*(x2**2)', '(x1**2)*x2*x3', 'x1*(x2**3)', 'x1*(x2**2)*x3', 'x1*x2*(x3**2)',
    'x2**4', '(x2**3)*x3', '(x2**2)*(x3**2)', 'x2*(x3**3)'
]


help3 = [
    'x3', 'x1*x3', 'x2*x3', 'x3**2', '(x1**2)*x3', '(x2**2)*x3', 'x3**3',
    'x1*x2*x3', 'x1*(x3**2)', 'x2*(x3**2)', '(x1**3)*x3', '(x1**2)*x2*x3',
    '(x1**2)*(x3**2)', 'x1*(x2**2)*x3', 'x1*x2*(x3**2)', 'x1*(x3**3)',
    '(x2**3)*x3', '(x2**2)*(x3**2)', 'x2*(x3**3)', 'x3**4'
]
help4 = [
    '1', 'x1', 'x2', 'x3', 'x1**2', 'x2**2', 'x3**2', 'x1*x2', 'x1*x3', 'x2*x3',
    'x1**3', '(x1**2)*x2', '(x1**2)*x3', 'x1*(x2**2)', 'x1*x2*x3', 'x1*(x3**2)',
    'x2**3', '(x2**2)*x3', 'x2*(x3**2)', 'x3**3'
]
SDP(dico_monomials,help1,help2,help3,help4,d0,d1,d2,"sdp_problem_relaxation_t=4.txt")
