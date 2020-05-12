# DISTRIBUIÇÃO DA ENERGIA CINÉTICA APÓS EXCITAÇÃO DA MOLÉCULA PELO MÉTODO DE REFLEXÃO

#_____________________________________________________________________________________________________________________________

# OBJETIVOS:

# Gerar o sorteio segundo a distribuição gaussiana em um array
# Calcular valores de energia do estado fundamental da molécula de H2 e de curvas de potencial de estados excitados
# Selecionar os valores do array gerado que estão dentro do intervalor da região de Franck-Condom
# Com o novo array, calcular os valores de energia para cada curva de potencial
# Somar a componente de energia rotacional para N = 2 (N = 0 já é considerado com o objetivo anterior)
# Salvar os dados gerados das distância internuclear e as energias de cada curva de potencial em arquivos csv
# Método da reflexão: realizar o ajuste polinomial na função R(E) e cálculo das distribuições de energia cinética
# Ajuste gaussiana nas distribuições de energia cinéticas dos fragmentos
# Configurações dos histogramas e gráficos relativos a cada variável importante do programa

# OBS: os pesos de cada curva de potencial são para o programa principal, mas vale lembrar que as chances de excitação para cada curva de potencial não são iguais

#_____________________________________________________________________________________________________________________________

# BIBLIOTECAS

import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy as scipy
from scipy.optimize import curve_fit

#_____________________________________________________________________________________________________________________________

# SORTEIO SEGUNDO DISTRIBUIÇÃO GAUSSIANA

# Esse sorteio é utilizado na simulação para sorteio da distância internuclear para o cálculo de Delta E (vide tese de doutorado da Aline Medina), que servirá também para o código para exemplificar o método de reflexão para obtenção da energia cinética dos fragmentos.

R_min = 0.6306353832413667                                                     # Valor mínimo da região de Franck-Condom (Angstrom)
R_max = 0.8872245404202529                                                     # Valor máximo da região de Franck-Condom (Angstrom)
conv = 27.2113845                                                              # Valor de conversão entre Hartree e eV

mu = 0.74212                                                                   # Média do distribuição gaussiana (Angstrom)
sigma = 0.1233                                                                 # Desvio padrão da distribuição gaussiana (Angstrom)
g = np.random.normal (mu, sigma, 1000)                                         # Comando que gera o array segundo a distribuição gaussiana, utiliza a biblioteca numpy

#_____________________________________________________________________________________________________________________________

# TESTE PARA COLOCAR QUALQUER EQUAÇÃO DA CURVA POTENCIAL E EQUAÇÕES PARA AJUSTES NO TERMINAL

# Curvas de Potencial

#nome_da_curva = input("Digite o nome da curva: ")

#A0 = input("Digite a amplitude máxima da exponencial: ")
#B0 = input("Digite o expoente da exponencial: ")
#C0 = input("Digite o termo constante: ")

#def Exp_curva_potencial(x, A0, B0, C0):
#    return A0 * np.exp(-x*B0) + C0

# Equações para Ajustes

#nome_da_curva = input("Digite o nome da curva: ")

#A0 = input("Digite a amplitude máxima da exponencial: ")
#B0 = input("Digite o expoente da exponencial: ")
#C0 = input("Digite o termo constante: ")

#def Equacao_de_ajuste(x, A0, B0, C0):
#    return A0 * np.exp(-x*B0) + C0

#_____________________________________________________________________________________________________________________________

# VALORES DE ENERGIA DO ESTADO FUNDAMENTAL DA MOLÉCULA DE H2 E DE ALGUMAS CURVAS DE POTENCIAL DE ESTADOS EXCITADOS

# Aqui constam as equações do estado fundamental (em eV) e das curvas de potencial (em ua), bem como os valores máximos e mínimos em cada curva de potencial após excitação

E_ground = 66.15298 - 338.88776 * g + 702.36623 * (g**2) - 748.93254 * (g**3) + 413.7916 * (g**4) - 94.06366 * (g**5)       # Equação do estado fundamental (Realizar gráfico de E x R)

# Pi 1

E_pi1_min_ua = 2.5751 * np.exp(-R_min/0.3375) + 0.9879                         # Equação da curva de potencial Pi 1 (com valor mínimo da região de Franck-Condom) em unidades atômicas
E_pi1_min = conv * E_pi1_min_ua                                                # Conversão de Hartree para eV
E_pi1_max_ua  = 2.5751 * np.exp(-R_max/0.3375) + 0.9879                        # Equação da curva de potencial Pi 1 (com valor máximo da região de Franck-Condom) em unidades atômicas
E_pi1_max = conv * E_pi1_max_ua                                                # Conversão de Hartree para eV

E_pi1_mu_ua = 2.5751 * np.exp(-mu/0.3375) + 0.9879                             # Equação da curva de potencial Pi 1 (com valor média do distribuição gaussiana) em unidades atômicas
E_pi1_mu = conv * E_pi1_mu_ua                                                  # Conversão de Hartree para eV
print("Energia de excitação da curva Pi 1:", "%.2f" % E_pi1_mu)

# Pi 2

E_pi2_min_ua = 2.5479 * np.exp(-R_min/0.3476) + 1.0526                         # Equação da curva de potencial Pi 2 (com valor mínimo da região de Franck-Condom) em unidades atômicas
E_pi2_min = conv * E_pi2_min_ua                                                # Conversão de Hartree para eV
E_pi2_max_ua = 2.5479 * np.exp(-R_max/0.3476) + 1.0526                         # Equação da curva de potencial Pi 2 (com valor máximo da região de Franck-Condom) em unidades atômicas
E_pi2_max = conv * E_pi2_max_ua                                                # Conversão de Hartree para eV

E_pi2_mu_ua = 2.5479 * np.exp(-mu/0.3476) + 1.0526                             # Equação da curva de potencial Pi 2 (com valor média do distribuição gaussiana) em unidades atômicas
E_pi2_mu = conv * E_pi2_mu_ua                                                  # Conversão de Hartree para eV
print("Energia de excitação da curva Pi 2:", "%.2f" % E_pi2_mu)

# SigmaG 1

E_sigmag1_min_ua = 2.517 * np.exp(-R_min/0.335) + 0.957                        # Equação da curva de potencial SigmaG 1 (com valor mínimo da região de Franck-Condom) em unidades atômicas
E_sigmag1_min = conv * E_sigmag1_min_ua                                        # Conversão de Hartree para eV
E_sigmag1_max_ua = 2.517 * np.exp(-R_max/0.335) + 0.957                        # Equação da curva de potencial SigmaG 1 (com valor máximo da região de Franck-Condom) em unidades atômicas
E_sigmag1_max = conv * E_sigmag1_max_ua                                        # Conversão de Hartree para eV

E_sigmag1_mu_ua = 2.517 * np.exp(-mu/0.335) + 0.957                            # Equação da curva de potencial SigmaG 1 (com valor média do distribuição gaussiana) em unidades atômicas
E_sigmag1_mu = conv * E_sigmag1_mu_ua                                          # Conversão de Hartree para eV
print("Energia de excitação da curva SigmaG 1:", "%.2f" % E_sigmag1_mu)

E_sigmag1_mu_teste = conv * E_sigmag1_mu_ua - 26.04129446                      # Conversão de Hartree para eV

# H+ H+

E_pp_min_ua = 1.61973 * np.exp((-R_min)*0.76995) + 1.32999                     # Equação da curva de potencial Próton-próton (com valor mínimo da região de Franck-Condom) em unidades atômicas
E_pp_min = conv * E_pp_min_ua                                                  # Conversão de Hartree para eV
E_pp_max_ua = 1.61973 * np.exp((-R_max)*0.76995) + 1.32999                     # Equação da curva de potencial Próton-próton (com valor máximo da região de Franck-Condom) em unidades atômicas
E_pp_max = conv * E_pp_max_ua                                                  # Conversão de Hartree para eV

E_pp_mu_ua = 1.61973 * np.exp((-mu)*0.76995) + 1.32999                         # Equação da curva de potencial Próton-próton (com valor média do distribuição gaussiana) em unidades atômicas
E_pp_mu = conv * E_pp_mu_ua                                                    # Conversão de Hartree para eV
print("Energia de excitação da curva próton-próton:", "%.2f" % E_pp_mu)

E_pp_mu_teste = conv * E_pp_mu_ua - 36.1908693                                 # Conversão de Hartree para eV

#print("Energia cinética dos fragmentos da curva SigmaG 1:", "%.2f" % E_sigmag1_mu_teste)
#print("Energia cinética dos fragmentos da curva próton-próton:", "%.2f" % E_pp_mu_teste)

#_____________________________________________________________________________________________________________________________

# SELEÇÃO DOS VALORES QUE ESTÃO NA REGIÃO DE FRANCK-CONDOM

gn = ()

for i in g:                                                                    # Iteração dos elementos de g 
    if R_min < i < R_max:                                                      # Condição que determinar quais valores estão dentro do intervalo
        gn += (i,)                                                             # Acréscimo de cada valor da iteração a tupla gn

#_____________________________________________________________________________________________________________________________ 
    
# CURVAS DE POTENCIAIS DOS ESTADOS EXCITADOS

# Pi 1

E_pi1_ua = ()
E_pi1 = ()
E_pi1_teste = ()

for j1 in gn:                                                                  # Iteração dos elementos da tupla gn segundo a equação da curva de potencial Piu 1 
    u1 = 2.5751 * np.exp(-j1/0.3375) + 0.9879                                  # Pi 1 peso 0,419
    t1 = conv * u1 - 26.88212623                                               # Conversão de Hartree para eV e subtração do valor assintótico da curva de potencial
    v1 = t1 + 26.88212623
    E_pi1_ua += (u1,)                                                          # Acréscimo de cada valor da iteração a tupla E_pi1_ua
    E_pi1 += (t1,)                                                             # Acréscimo de cada valor da iteração a tupla E_pi1
    E_pi1_teste += (v1,)

# Pi 2

E_pi2_ua = ()
E_pi2 = ()
E_pi2_teste = ()

for j2 in gn:                                                                  # Iteração dos elementos da tupla gn segundo a equação da curva de potencial Piu 2
    u2 = 2.5479 * np.exp(-j2/0.3476) + 1.0526                                  # Pi 2 peso 0,379
    t2 = conv * u2 - 28.64270277                                               # Conversão de Hartree para eV e subtração do valor assintótico da curva de potencial
    v2 = t2 + 28.64270277
    E_pi2_ua += (u2,)                                                          # Acréscimo de cada valor da iteração a tupla E_pi2_ua
    E_pi2 += (t2,)                                                             # Acréscimo de cada valor da iteração a tupla E_pi2
    E_pi2_teste += (v2,)

# SigmaG 1

E_sigmag1_ua = ()
E_sigmag1 = ()
E_sigmag1_teste = ()

for j3 in gn:                                                                  # Iteração dos elementos da tupla gn segundo a equação da curva de potencial Sigma G 1
    u3 = 2.517 * np.exp(-j3/0.335) + 0.957                                     # Sigma G 1 peso 0,203
    t3 = conv * u3 - 26.04129446                                               # Conversão de Hartree para eV e subtração do valor assintótico da curva de potencial
    v3 = t3 + 26.04129446
    E_sigmag1_ua += (u3,)                                                      # Acréscimo de cada valor da iteração a tupla E_sigmag1_ua
    E_sigmag1 += (t3,)                                                         # Acréscimo de cada valor da iteração a tupla E_sigmag1
    E_sigmag1_teste += (v3,)

# H+ H+

E_pp_ua = ()
E_pp = ()
E_pp_teste = ()

for j4 in gn:                                                                  # Iteração dos elementos da tupla gn segundo a equação da curva de potencial Próton-próton
    u4 = 1.61973 * np.exp((-j4)*0.76995) + 1.32999                             # Próton - Próton não sei o peso dessa curva de potencial
    t4 = conv * u4 - 36.1908693                                                # Conversão de Hartree para eV e subtração do valor assintótico da curva de potencial
    v4 = t4 + 36.1908693
    E_pp_ua += (u4,)                                                           # Acréscimo de cada valor da iteração a tupla E_pp_ua
    E_pp += (t4,)                                                              # Acréscimo de cada valor da iteração a tupla E_pp
    E_pp_teste += (v4,)

#_____________________________________________________________________________________________________________________________
    
# SOMA DA ENERGIA ROTACIONAL COM A ENERGIA CINÉTICA DE CADA CURVA

N2 = 2
h_bar = 6.582119514 * (10**(-16))
mr = 0.83025 * (10**(-27))

# Energia rotacional para N = 2

E_rot = ()

for k0 in E_sigmag1:                                                           # Iteração dos elementos da tupla E_rot
    er = ((h_bar**2) * N2 * (N2 + 1))/(2 * mr * k0)
    E_rot += (er,)                                                             # Acréscimo de cada valor da iteração a tupla E_rot

# Pi 1 mais efeito da rotação da molécula

E_pi1_tot = ()

for k1, k0 in zip(E_pi1, E_rot):                                               # Iteração dos elementos, um a um, das tuplas E_pi1 e E_rot para obter a energia total da curva Pi 1
    e1 = k1 + k0
    E_pi1_tot += (e1,)                                                         # Acréscimo de cada valor da iteração a tupla E_pi1_tot
    
# Pi 2 mais efeito da rotação da molécula

E_pi2_tot = ()

for k2, k0 in zip(E_pi2, E_rot):                                               # Iteração dos elementos, um a um, das tuplas E_pi2 e E_rot para obter a energia total da curva Pi 2
    e2 = k2 + k0
    E_pi2_tot += (e2,)                                                         # Acréscimo de cada valor da iteração a tupla E_pi2_tot

# SigmaG 1 mais efeito da rotação da molécula

E_sigmag1_tot = ()

for k3, k0 in zip(E_sigmag1, E_rot):                                           # Iteração dos elementos, um a um, das tuplas E_sigmag1 e E_rot para obter a energia total da curva Sigma G 1
    e3 = k3 + k0
    E_sigmag1_tot += (e3,)                                                     # Acréscimo de cada valor da iteração a tupla E_sigmag1_tot
  
#_____________________________________________________________________________________________________________________________

# SALVANDO OS DADOS GERADOS DE DISTÂNCIA INTERNUCLEAR E AS ENERGIAS DE CADA CURVA DE POTENCIAL
    
# Pode ser utilizado para comparar os valores de ajustes polinomiais do Origin ou Mathematica (Os ajustes, pelos meus teste, são bem confiáveis)

# Distância Internuclear

#with open("dados_dist_intern.csv", 'w', newline='') as saida1:
#    escrever = csv.writer(saida1)
#    for ld in gn:
#        escrever.writerow([ld])

# Energia de Pi 1

#with open("dados_E_pi1.csv", 'w', newline='') as saida2:
#    escrever = csv.writer(saida2)
#    for l1 in E_pi1:
#        escrever.writerow([l1])

# Energia de Pi 2
        
#with open("dados_E_pi2.csv", 'w', newline='') as saida3:
#    escrever = csv.writer(saida3)
#    for l2 in E_pi2:
#        escrever.writerow([l2])
        
# Energia de SigmaG 1
        
#with open("dados_E_sigmag1.csv", 'w', newline='') as saida4:
#    escrever = csv.writer(saida4)
#    for l3 in E_sigmag1:
#        escrever.writerow([l3])

# Energia de H+ H+
        
#with open("dados_E_pp.csv", 'w', newline='') as saida5:
#    escrever = csv.writer(saida5)
#    for l4 in E_pp:
#        escrever.writerow([l4])

#_____________________________________________________________________________________________________________________________

# MÉTODO DA REFLEXÃO (TESTE)

# Inicialmente, os dados foram tratados no software Origin a fim de obter o ajuste de polinômio de terceiro grau para o gráfico R(E)

# Função do Ajuste Polinomial

def polinomial_3(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d                                           # Definição da função para o ajuste da curva de potencial com os valores da região de Franck-Condom

def exponencial(x, a, R_0, b):
    return a * np.exp(R_0*x) + b                                               # Definição da função para o ajuste da curva de potencial com os valores da região de Franck-Condom

# Pi 1

param_pi1, param_cov_pi1 = curve_fit(polinomial_3, E_pi1, gn)                  # Comando que realiza o ajuste com a função pré definida para Pi 1

print("Coeficientes da função polinomial de terceiro grau para Pi 1:")
print("a: ", param_pi1[0])
print("b: ", param_pi1[1])
print("c: ", param_pi1[2])
print("d: ", param_pi1[3])
#print("Covariância da função polinomial de terceiro grau:")
#print(param_cov_pi1)

# Pi 1 + E_rot

param_pi1_tot, param_cov_pi1_tot = curve_fit(polinomial_3, E_pi1_tot, gn)      # Comando que realiza o ajuste com a função pré definida para Pi 1 + E_rot

#print("Coeficientes da função polinomial de terceiro grau para Pi 1:")
#print("a: ", param_pi1_tot[0])
#print("b: ", param_pi1_tot[1])
#print("c: ", param_pi1_tot[2])
#print("d: ", param_pi1_tot[3])
#print("Covariância da função polinomial de terceiro grau:")
#print(param_cov_pi1_tot)

# Pi 2

param_pi2, param_cov_pi2 = curve_fit(polinomial_3, E_pi2, gn)                  # Comando que realiza o ajuste com a função pré definida para Pi 2

print("Coeficientes da função polinomial de terceiro grau para Pi 2:")
print("a: ", param_pi2[0])
print("b: ", param_pi2[1])
print("c: ", param_pi2[2])
print("d: ", param_pi2[3])
#print("Covariância da função polinomial de terceiro grau:")
#print(param_cov_pi2)

# Pi 2 + E_rot

param_pi2_tot, param_cov_pi2_tot = curve_fit(polinomial_3, E_pi2_tot, gn)      # Comando que realiza o ajuste com a função pré definida para Pi 2 + E_rot

#print("Coeficientes da função polinomial de terceiro grau para Pi 2 + E_rot:")
#print("a: ", param_pi2_tot[0])
#print("b: ", param_pi2_tot[1])
#print("c: ", param_pi2_tot[2])
#print("d: ", param_pi2_tot[3])
#print("Covariância da função polinomial de terceiro grau:")
#print(param_cov_pi2_tot)

# SigmaG 1

param_sigmag1, param_cov_sigmag1 = curve_fit(polinomial_3, E_sigmag1, gn)      # Comando que realiza o ajuste com a função pré definida para SigmaG 1

print("Coeficientes da função polinomial de terceiro grau para SigmaG 1:")
print("a: ", param_sigmag1[0])
print("b: ", param_sigmag1[1])
print("c: ", param_sigmag1[2])
print("d: ", param_sigmag1[3])
#print("Covariância da função polinomial de terceiro grau:")
#print(param_cov_sigmag1_tot)

# SigmaG 1 + E_rot

param_sigmag1_tot, param_cov_sigmag1_tot = curve_fit(polinomial_3, E_sigmag1_tot, gn)        # Comando que realiza o ajuste com a função pré definida para SigmaG 1

#print("Coeficientes da função polinomial de terceiro grau para SigmaG 1 + E_rot:")
#print("a: ", param_sigmag1_tot[0])
#print("b: ", param_sigmag1_tot[1])
#print("c: ", param_sigmag1_tot[2])
#print("d: ", param_sigmag1_tot[3])
#print("Covariância da função polinomial de terceiro grau:")
#print(param_cov_sigmag1_tot)

# H+ H+

param_pp, param_cov_pp = curve_fit(polinomial_3, E_pp, gn)                     # Comando que realiza o ajuste com a função pré definida para próton-próton

print("Coeficientes da função polinomial de terceiro grau para próton-próton:")
print("a: ", param_pp[0])
print("b: ", param_pp[1])
print("c: ", param_pp[2])
print("d: ", param_pp[3])
#print("Covariância da função polinomial de terceiro grau:")
#print(param_cov_pp)

# Método da Reflexão

E = np.linspace(0, 32, num=1000)                                               # Cria um vetor com valores igualmente espaçados dentro do intervalo escolhido

E_1 = np.linspace(0, 11, num=1000)                                             # Cria um vetor com valores igualmente espaçados dentro do intervalo escolhido

# Pi 1

# Para dois fragmentos

R0_pi1 = param_pi1[0]*(E**3) + param_pi1[1]*(E**2) + param_pi1[2]*(E) + param_pi1[3]        # Função inversa R(E)
dR0dE_pi1 = 3*param_pi1[0]*(E**2) + 2*param_pi1[1]*(E) + param_pi1[2]                       # Derivada de R(E)

f0_pi1 = 2.139 * np.exp(-32.882 * (R0_pi1 - 0.74212)**2)                                    # Função de onda
gn0_pi1 = (abs(f0_pi1))**2 *(abs(dR0dE_pi1))                                                # Distribuição de energia cinética


# Para um só fragmento

R0_pi1_1 = param_pi1[0]*(2*E_1**3) + param_pi1[1]*(2*E_1**2) + param_pi1[2]*(2*E_1) + param_pi1[3]        # Função inversa R(E)
dR0dE_pi1_1 = 3*param_pi1[0]*(2*E_1**2) + 2*param_pi1[1]*(2*E_1) + param_pi1[2]                           # Derivada de R(E)

f0_pi1_1 = 2.139 * np.exp(-32.882 * (R0_pi1_1 - 0.74212)**2)                                # Função de onda
gn0_pi1_1 = (abs(f0_pi1_1))**2 *(abs(dR0dE_pi1_1))                                          # Distribuição de energia cinética

# Pi 2

# Para dois fragmentos

R0_pi2 = param_pi2[0]*(E**3) + param_pi2[1]*(E**2) + param_pi2[2]*(E) + param_pi2[3]        # Função inversa R(E)
dR0dE_pi2 = 3*param_pi2[0]*(E**2) + 3*param_pi2[0]*(E) + param_pi2[2]                       # Derivada de R(E)

f0_pi2 = 2.139 * np.exp(-32.882 * (R0_pi2 - 0.74212)**2)                                    # Função de onda
gn0_pi2 = (abs(f0_pi2))**2 *(abs(dR0dE_pi2))                                                # Distribuição de energia cinética

# Para um só fragmento

R0_pi2_1 = param_pi2[0]*(2*E_1**3) + param_pi2[1]*(2*E_1**2) + param_pi2[2]*(2*E_1) + param_pi2[3]        # Função inversa R(E)
dR0dE_pi2_1 = 3*param_pi2[0]*(2*E_1**2) + 3*param_pi2[0]*(2*E_1) + param_pi2[2]                           # Derivada de R(E)

f0_pi2_1 = 2.139 * np.exp(-32.882 * (R0_pi2_1 - 0.74212)**2)                                # Função de onda
gn0_pi2_1 = (abs(f0_pi2_1))**2 *(abs(dR0dE_pi2_1))                                          # Distribuição de energia cinética

# SigmaG 1

# Para dois fragmentos

R0_sigmag1 = param_sigmag1[0]*(E**3) + param_sigmag1[1]*(E**2) + param_sigmag1[2]*(E) + param_sigmag1[3]        # Função inversa R(E)
dR0dE_sigmag1 = 3*param_sigmag1[0]*(E**2) + 2*param_sigmag1[1]*(E) + param_sigmag1[2]                           # Derivada de R(E)

f0_sigmag1 = 2.139 * np.exp(-32.882 * (R0_sigmag1 - 0.74212)**2)                            # Função de onda
gn0_sigmag1 = (abs(f0_sigmag1))**2 *(abs(dR0dE_sigmag1))                                    # Distribuição de energia cinética

# Para um só fragmento

R0_sigmag1_1 = param_sigmag1[0]*(2*E_1**3) + param_sigmag1[1]*(2*E_1**2) + param_sigmag1[2]*(2*E_1) + param_sigmag1[3]        # Função inversa R(E)
dR0dE_sigmag1_1 = 3*param_sigmag1[0]*(2*E_1**2) + 2*param_sigmag1[1]*(2*E_1) + param_sigmag1[2]                               # Derivada de R(E)

f0_sigmag1_1 = 2.139 * np.exp(-32.882 * (R0_sigmag1_1 - 0.74212)**2)                        # Função de onda
gn0_sigmag1_1 = (abs(f0_sigmag1_1))**2 *(abs(dR0dE_sigmag1_1))                              # Distribuição de energia cinética


# H+ H+

# Para dois fragmentos

R0_pp = param_pp[0]*(E**3) + param_pp[1]*(E**2) + param_pp[2]*(E) + param_pp[3]             # Função inversa R(E)
dR0dE_pp = 3*param_pp[0]*(E**2) + 2*param_pp[1]*(E) + param_pp[2]                           # Derivada de R(E)

f0_pp = 2.139 * np.exp(-32.882 * (R0_pp - 0.74212)**2)                                      # Função de onda
gn0_pp = (abs(f0_pp))**2 *(abs(dR0dE_pp))                                                   # Distribuição de energia cinética

# Para um só fragmento

R0_pp_1 = param_pp[0]*(2*E_1**3) + param_pp[1]*(2*E_1**2) + param_pp[2]*(2*E_1) + param_pp[3]        # Função inversa R(E)
dR0dE_pp_1 = 3*param_pp[0]*(2*E_1**2) + 2*param_pp[1]*(2*E_1) + param_pp[2]                          # Derivada de R(E)

f0_pp_1 = 2.139 * np.exp(-32.882 * (R0_pp_1 - 0.74212)**2)                                  # Função de onda
gn0_pp_1 = (abs(f0_pp_1))**2 *(abs(dR0dE_pp_1))                                             # Distribuição de energia cinética

#_____________________________________________________________________________________________________________________________

# AJUSTE GAUSSIANA NAS DISTRIBUIÇÕES DE ENERGIA CINÉTICA

#_____________________________________________________________________________________________________________________________

# CONFIGURAÇÃO DO HISTOGRAMA SEGUNDO A DISTRIBUIÇÃO GAUSSIANA (somente utilizado para compreensão do programa, na simulação não é necessária)

#count, bins, ignored = plt.hist(g, 40, density=True, histtype = 'stepfilled', facecolor = 'blue')                           # Comando que gera o histograma
#plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=2, color='red')      # Plota um ajuste segundo a equação colocada
#plt.axvline(R_min, color = 'red')
#plt.axvline(R_max, color = 'red')
#plt.title('Distribuição Gaussiana')                                                                                         # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('distribuicao_gaussiana.png')                                                                                   # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#_____________________________________________________________________________________________________________________________

# CONFIGURAÇÃO DO GRÁFICO DE ENERGIA DO ESTADO FUNDAMENTAL PELA DISTÂNCIA INTERNUCLEAR (somente utilizado para compreensão do programa, na simulação não é necessária)

#plt.scatter(g, E_ground, color='blue')                                                                                      # Plota pontos verdes no gráfico
#plt.title('Estado Fundamental da Molécula de H2')                                                                           # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('metodo_reflexao_1SIGMA.png')                                                                                   # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#_____________________________________________________________________________________________________________________________

# CONFIGURAÇÃO DOS GRÁFICOS ENERGIA CINÉTICA APÓS A EXCITAÇÃO DE CADA CURVA DE POTENCIAL SEM A ENERGIA ROTACIONAL (somente utilizado para compreensão do programa, na simulação não é necessária)

#plt.scatter(gn, E_pi1, color='blue')                                                                                        # Comando que gera o gráfico
#plt.title('E_PiU1(R)')                                                                                                      # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('1PIU_sem_rot.png')                                                                                             # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(gn, E_pi2, color='blue')                                                                                        # Comando que gera o gráfico
#plt.title('E_PiU2(R)')                                                                                                      # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('2PIU_sem_rot.png')                                                                                             # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(gn, E_sigmag1, color='blue')                                                                                    # Comando que gera o gráfico
#plt.title('E_SigmaG1(R)')                                                                                                   # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('1SIGMAG_sem_rot.png')                                                                                          # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(gn, E_pp, color='blue')                                                                                         # Comando que gera o gráfico
#plt.title('Próton-próton')                                                                                                  # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('pp_sem_rot.png')                                                                                               # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal


#_____________________________________________________________________________________________________________________________

# CONFIGURAÇÃO DOS GRÁFICOS ENERGIA CINÉTICA APÓS A EXCITAÇÃO DE CADA CURVA DE POTENCIAL COM A ENERGIA ROTACIONAL (somente utilizado para compreensão do programa, na simulação não é necessária)

#plt.scatter(gn, E_pi1_tot, color='blue')                                                                                    # Comando que gera o gráfico
#plt.title('E_PiU1(R) + E_rot(R)')                                                                                           # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('1PIU_com_rot.png')                                                                                             # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(gn, E_pi2_tot, color='blue')                                                                                    # Comando que gera o gráfico
#plt.title('E_PiU2(R) + E_rot(R)')                                                                                           # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('2PIU_com_rot.png')                                                                                             # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal
 
#plt.scatter(gn, E_sigmag1_tot, color='blue')                                                                                # Comando que gera o gráfico
#plt.title('E_SigmaG1(R) + E_rot(R)')                                                                                        # Definindo título
#plt.xlabel('Distância Internuclear (Angstrom)')                                                                             # Definindo nome do eixo X
#plt.ylabel('Energia (eV)')                                                                                                  # Definindo nome do eixo Y
#plt.savefig('1SIGMAG_com_rot.png')                                                                                          # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#_____________________________________________________________________________________________________________________________

# CONFIGURAÇÃO DO GRÁFICO DE TODAS AS VARIÁVEIS (somente utilizado para compreensão do programa, na simulação não é necessária)

plt.plot(g, E_ground, 'o', label='Fundamental')
plt.plot(gn, E_pi1_teste, '--', label='Pi1')
plt.plot(gn, E_pi2_teste, '--', label='Pi2')
plt.plot(gn, E_sigmag1_teste, '--', label='SigmaG1')
plt.plot(gn, E_pp_teste, '--', label='H+ H+')
plt.title('Gráfico de Energia vs Distância Internuclear')                                                                    # Definindo título
plt.xlabel('Distância Internuclear (Angstrom)')                                                                              # Definindo nome do eixo X
plt.ylabel('Energia (eV)')                                                                                                   # Definindo nome do eixo Y
plt.legend()                                                                                                                 # Imprime a legenda
plt.savefig('tudo.png')                                                                                                      # Salva imagem
plt.show()                                                                                                                   # Imprime no terminal

#_____________________________________________________________________________________________________________________________

# CONFIGURAÇÃO DO GRÁFICO DA DISTRIBUIÇÃO DE ENERGIA CINÉTICA (somente utilizado para compreensão do programa, na simulação não é necessária)

#plt.scatter(E, gn0_pi1, color='blue')
#plt.title('Distribuição de Energia Cinética para Pi 1')                                                                     # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_pi1.png')                                                                                           # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(E_1, gn0_pi1_1, color='blue')
#plt.title('Distribuição de Energia Cinética para Pi 1 um fragmento')                                                        # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_pi1_um_frag.png')                                                                                   # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(E, gn0_pi2, color='blue')
#plt.title('Distribuição de Energia Cinética para Pi 2')                                                                     # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_pi2.png')                                                                                           # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(E_1, gn0_pi2_1, color='blue')
#plt.title('Distribuição de Energia Cinética para Pi 2 um fragmento')                                                        # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_pi2_um_frag.png')                                                                                   # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(E, gn0_sigmag1, color='blue')
#plt.title('Distribuição de Energia Cinética para SigmaG 1')                                                                 # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_sigmag1.png')                                                                                       # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(E_1, gn0_sigmag1_1, color='blue')
#plt.title('Distribuição de Energia Cinética para SigmaG 1 um fragmento')                                                    # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_sigmag1_um_frag.png')                                                                               # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(E, gn0_pp, color='blue')
#plt.title('Distribuição de Energia Cinética para p-p')                                                                      # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_pp.png')                                                                                            # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

#plt.scatter(E_1, gn0_pp_1, color='blue', label='Pi1')
#plt.title('Distribuição de Energia Cinética para p-p um fragmento')                                                         # Definindo título
#plt.xlabel('Energia Cinética (eV)')                                                                                         # Definindo nome do eixo X
#plt.ylabel('Densidade de Probabilidade')                                                                                    # Definindo nome do eixo Y
#plt.savefig('dist_energ_pp_um_frag.png')                                                                                    # Salva imagem
#plt.show()                                                                                                                  # Imprime no terminal

plt.plot(E, gn0_pi1, 'ro', label='Pi1')
plt.plot(E, gn0_pi2, 'go', label='Pi2')
plt.plot(E, gn0_sigmag1, 'bo', label='SigmaG1')
plt.plot(E, gn0_pp, 'o', label='H+ H+')
plt.title('Pi1, Pi2, SigmaG1 e p-p')                                                                                         # Definindo título
plt.xlabel('Energia Cinética (eV)')                                                                                          # Definindo nome do eixo X
plt.ylabel('Densidade de Probabilidade')                                                                                     # Definindo nome do eixo Y
plt.legend()                                                                                                                 # Imprime a legenda
plt.savefig('pi1_pi2_sigmag1_H+H+_2frag.png')                                                                                                     # Salva imagem
plt.show()                                                                                                                   # Imprime no terminal

plt.plot(E, gn0_sigmag1, 'o', label='SigmaG1')
plt.plot(E, gn0_pp, 'o', label='H+ H+')
plt.title('Distribuição de Energia Cinética para SigmaG 1 e p-p')                                                            # Definindo título
plt.xlabel('Energia Cinética (eV)')                                                                                          # Definindo nome do eixo X
plt.ylabel('Densidade de Probabilidade')                                                                                     # Definindo nome do eixo Y
plt.legend()                                                                                                                 # Imprime a legenda
plt.savefig('sigmag1_H+H+_2frag.png')                                                                                                # Salva imagem
plt.show()                                                                                                                   # Imprime no terminal

plt.plot(E_1, gn0_pi1_1, 'ro', label='Pi1')
plt.plot(E_1, gn0_pi2_1, 'go', label='Pi2')
plt.plot(E_1, gn0_sigmag1_1, 'bo', label='SigmaG1')
plt.plot(E_1, gn0_pp_1, 'o', label='H+ H+')
plt.title('Pi1, Pi2 , SigmaG 1 e p-p um fragmento')                                                                          # Definindo título
plt.xlabel('Energia Cinética (eV)')                                                                                          # Definindo nome do eixo X
plt.ylabel('Densidade de Probabilidade')                                                                                     # Definindo nome do eixo Y
plt.legend()                                                                                                                 # Imprime a legenda
plt.savefig('pi1_pi2_sigmag1_H+H+_1frag.png')                                                                                             # Salva imagem
plt.show()                                                                                                                   # Imprime no terminal

plt.plot(E_1, gn0_sigmag1_1, 'o', label='SigmaG1')
plt.plot(E_1, gn0_pp_1, 'o', label='H+ H+')
plt.title('SigmaG 1 e p-p um fragmento')                                                                                     # Definindo título
plt.xlabel('Energia Cinética (eV)')                                                                                          # Definindo nome do eixo X
plt.ylabel('Densidade de Probabilidade')                                                                                     # Definindo nome do eixo Y
plt.legend()                                                                                                                 # Imprime a legenda
plt.savefig('sigmag1_H+H+_1frag.png')                                                                                        # Salva imagem
plt.show()                                                                                                                   # Imprime no terminal