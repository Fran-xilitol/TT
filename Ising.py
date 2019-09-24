from random import randint
from random import random
from math import e
import csv

def Ising(L):
    n = 10
    #n passos Monte Carlo

    data = []
    data.append(["T,u,m,c,chi"])

    lat = [] #linhas da rede
    i = 0
    while i < L:
        cell = [] #colunas da rede
        j = 0
        while j < L:
            cell.append(randint(0,1)*2-1) #momento magnético de cada célula
            j += 1
        lat.append(cell)
        i += 1
    #condição periódica de contorno:
    for c in lat:
        c.append(c[0]) #última coluna é igual à primeira
    lat.append(lat[0]) #última linha é igual à primeira
    
    #energia e momento magnético da rede
    M = 0
    E = 0
    j = 0
    while j < L:
        k = 0
        while k < L:
            M += lat[j][k]
            E -= lat[j][k]*lat[j+1][k] + lat[j][k]*lat[j][k+1]
            k += 1
        j += 1

    T = 1
    while T<4.1:
        Delta = 0 #variação de energia
        p = 0 #ativação térmica exp(1/kT)
        M1 = 0 #média da magnetização
        M2 = 0 #média dos quadrados de M
        E1 = 0 #média da energia
        E2 = 0 #média dos quadrados de E
        B = 1/T #para os próximos cálculos serem multiplicações em vez de divisões
        

        #minimização por Monte Carlo
        step = 0
        while step <= n:
            s = 0
            Ri = 0
            Rj= 0
            while s < L*L:
                Ri = randint(0,L-1)
                Rj = randint(0,L-1)
                #selecionando um sítio aleatório da rede
                #cálculos diferentes se o sítio está no interior ou na borda

                if Ri!=0 and Rj!=0:
                    Delta = 2*lat[Ri][Rj]*(lat[Ri-1][Rj] + lat[Ri][Rj-1] + lat[Ri][Rj+1] + lat[Ri+1][Rj])
                elif Ri!=0 and Rj==0:
                    Delta = 2*lat[Ri][Rj]*(lat[Ri-1][Rj] + lat[Ri][L-1] + lat[Ri+1][Rj] + lat[Ri][Rj+1])
                elif Ri==0 and Rj!=0:
                    Delta = 2*lat[Ri][Rj]*(lat[L-1][Rj] + lat[Ri][Rj-1]+ lat[Ri+1][Rj] + lat[Ri][Rj+1])
                elif Ri==0 and Rj==0:
                    Delta = 2*lat[Ri][Rj]*(lat[L-1][Rj] + lat[Ri][L-1]+ lat[Ri+1][Rj] + lat[Ri][Rj+1])

                if Delta <= 0:
                    lat[Ri][Rj] = -1*lat[Ri][Rj]
                    E = E + Delta
                    M = M + 2*lat[Ri][Rj]
                    
                    #não esquecer as condições de contorno
                    if Ri==0:
                        lat[L][Rj]=lat[0][Rj]
                    if Rj==0:
                        lat[Ri][L]=lat[Ri][0]

                else:
                    p = e**(-1*Delta*B)
                    chance = random()
                    if p > chance:
                        lat[Ri][Rj] = -1*lat[Ri][Rj]
                        E = E + Delta
                        M = M + 2*lat[Ri][Rj]
                    
                        #não esquecer as condições de contorno
                        if Ri==0:
                            lat[L][Rj]=lat[0][Rj]
                        if Rj==0:
                            lat[Ri][L]=lat[Ri][0]
                

                s += 1

            
            #descartando os primeiros 10% de passos para as médias
            if step > 0.1*n:
                E1 = E1 + E
                M1 = M1 + (M**2)**(0.5)
                E2 = E2 + E**2
                M2 = M2 + M**2

            step += 1

        E1 = E1/(0.9*n)
        M1 = M1/(0.9*n)
        E2 = E2/(0.9*n)
        M2 = M2/(0.9*n)

        data.append([T, E1/(L*L), M1/(L*L), (E2-(E1*E1))/(T*T)/(L*L), (M2-(M1*M1))/(T)/(L*L)])
        
        T += 0.1

    #exportando os dados para o arquivo
    with open('Rede de L = {}.csv'.format(L), 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(data)
    csvFile.close()
    return "Cálculos realizados com sucesso, ufa!"
print(Ising(50))