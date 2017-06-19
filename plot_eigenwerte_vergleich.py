import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp


Energien=np.genfromtxt('build/Durchlaufende_Energien.txt')
Potential=np.genfromtxt('build/Durchlaufende_Potentiale.txt')
frequenzzahler , frequenznenner =np.genfromtxt('build/Durchlaufende_Frequenzen.txt',unpack=True)
Gitterkonstante,Anzahl,Phasenverschiebung=np.genfromtxt('build/Parameter_fur_a='+str(int(Potential[0]))+'_E='+str(int(Energien[0]))+'_w='+ str(int(frequenzzahler[0]))+'%'+ str(int(frequenznenner[0]))+'a.txt',unpack=True)
Frequenz=np.array(frequenzzahler/frequenznenner)

#Eigenwerte=np.genfromtxt('build/Eigenwerte_fur_a=50_E=10_w=1%3a.txt')

def eig_erwartet_Matrix(Eigenwerte_H_0,Anzahl,Frequenz,Potential):
    m=np.size(Eigenwerte_H_0)*(1+Anzahl*2)
    Matrix=np.zeros((m,np.size(Frequenz)))
    Frequenz=Frequenz*Potential/100
    for i in range(np.size(Frequenz)):
        Erwartete_Eigenwerte=Eigenwerte_H_0
        for j in range(1,int(Anzahl+1)):
            Erwartete_Eigenwerte=np.append(Erwartete_Eigenwerte,Eigenwerte_H_0-Frequenz[i])
            Erwartete_Eigenwerte=np.append(Erwartete_Eigenwerte,Eigenwerte_H_0+Frequenz[i])
        Matrix[:,i]=np.sort(Erwartete_Eigenwerte)
    return Matrix


def Eigenwerte_Matrix(frequenzzahler,frequenznenner,Potential,Energie):
        Grosse=np.genfromtxt('build/Eigenwerte_fur_a='+str(int(Potential))+'_E='+str(int(Energie))+'_w='+ str(int(frequenzzahler[0]))+'%'+ str(int(frequenznenner[0]))+'a.txt')
        Matrix=np.zeros((np.size(Grosse),np.size(Frequenz)))
        for i in range(np.size(frequenznenner)):
            Matrix[:,i]=np.genfromtxt('build/Eigenwerte_fur_a='+str(int(Potential))+'_E='+str(int(Energie))+'_w='+ str(int(frequenzzahler[i]))+'%'+ str(int(frequenznenner[i]))+'a.txt')
        return Matrix

Eigenwerte_H_0=np.array([-2,0,0,2])
Figure_Zahler=1

for index_Potential , value_Potenial in enumerate(Potential):
    for index_Energie , value_Energie in enumerate(Energien):
        plt.figure(Figure_Zahler)
        Figure_Zahler=Figure_Zahler+1
        x=np.linspace(0,np.amax(Frequenz),100)
        Matrix_mit_Eigenwerten=Eigenwerte_Matrix(frequenzzahler,frequenznenner,value_Potenial,value_Energie)
        Matrix_mit_Erwarteteneigenwerten=eig_erwartet_Matrix(Eigenwerte_H_0,Anzahl,x,value_Potenial)
        n, m=np.shape(Matrix_mit_Eigenwerten)
        for j in range(n):
            plt.plot(Frequenz,Matrix_mit_Eigenwerten[j,:],'-b')#,label='Eigenwert'+ str(j) )
        n_2, m_2=np.shape(Matrix_mit_Eigenwerten)
        for i in range(n_2):
            plt.plot(x,Matrix_mit_Erwarteteneigenwerten[i,:],'--k')#,label='Eigenwert'+ str(j) )
        plt.legend(loc='best')
        plt.xlabel('Frequenz / a')
        plt.ylabel('E/J')
        plt.savefig('Plots/Plot_fur'+'_a='+str(value_Potenial/100)+'_E='+str(value_Energie/100)+'.pdf')
        plt.close()



        #     Erwartete_Eigenwerte=eig_war(Eigenwerte_H_0,Anzahl,Frequenz,value_Potenial)
        #     n=np.size(Eigenwerte)
        #     for i in range(n-1):
        #         plt.plot(Frequenz,Eigenwerte[i],'xb')
        #         plt.plot(Frequenz,Erwartete_Eigenwerte[i],'+r')
        #     plt.plot(Frequenz,Eigenwerte[n-1],'xb',label=r'E='+str(value_Energie/100))
        #     plt.plot(Frequenz,Erwartete_Eigenwerte[n-1],'+r',label=r'Erwartete_Eigenwerte')
        #     print('Potential',value_Potenial/100)
        #     print('value_Energie',value_Energie/100)
        #     print('Frequenz=', Frequenz)
        #     print('Erwartete_Eigenwerte',np.sort(Erwartete_Eigenwerte))
        #     print('Berechnete Eigenwerte',Eigenwerte)
