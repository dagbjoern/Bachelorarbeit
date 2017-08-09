import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from tqdm import tqdm


Energien=np.genfromtxt('build/Durchlaufende_Energien.txt')
Potential=np.genfromtxt('build/Durchlaufende_Potentiale.txt')
Frequenz_10000=np.genfromtxt('build/Durchlaufende_Frequenzen.txt',unpack=True)
Anzahl_N=np.genfromtxt('build/Durchlaufende_N.txt',unpack=True)

#Gitterkonstante,Anzahl,Phasenverschiebung=np.genfromtxt('Parameter/Parameter_fur_a='+str(int(Potential[0]))+'_w='+ str(int(Frequenz_10000[0]))+'.txt',unpack=True)





# if [(np.size(Frequenz_1000)==1)and(np.size(Potential)==1)]:
#     Gitterkonstante,Anzahl,Phasenverschiebung=np.genfromtxt('Parameter/Parameter_fur_a='+str(int(Potential))+'_w='+ str(int(Frequenz_1000))+'.txt',unpack=True)
# else:

Frequenz=Frequenz_10000/10000
#
# plt.figure(Figure_Zahler)
# Figure_Zahler=Figure_Zahler+1
# plt.title('Eigenwerte von a='+str(value_Potential/10000)+' w='+str(value_Frequenz)+' N='+str(Anzahl))
# plt.plot(x,x*0+value_Frequenz/2,'--k',alpha=0.5)
# plt.legend(loc='best')
# plt.xlabel('E')
# plt.ylabel(r'$\epsilon_\alpha $')
# plt.savefig('Plots/Plot_fur'+'_a='+str(value_Potential/10000)+'_w='+str(value_Frequenz)+'N='+str(Anzahl) +'.pdf')


#Eigenwerte=np.genfromtxt('build/Eigenwerte_fur_a=50_E=10_w=1%3a.txt')

def eig_erwartet_Matrix(Anzahl,Frequenz,Potential):
    Eigenwerte_H_0=np.genfromtxt('Parameter/eigenwerte_von_H_0_fur_a='+str(int(Potential)) +'a.txt')
    m=np.size(Eigenwerte_H_0)*(1+Anzahl*2)
    Matrix=np.zeros((m,np.size(Frequenz)))
    Frequenz=Frequenz
    for i in range(np.size(Frequenz)):
        Erwartete_Eigenwerte=Eigenwerte_H_0
        for j in range(1,int(Anzahl+1)):
            Erwartete_Eigenwerte=np.append(Erwartete_Eigenwerte,Eigenwerte_H_0-Frequenz[i])
            Erwartete_Eigenwerte=np.append(Erwartete_Eigenwerte,Eigenwerte_H_0+Frequenz[i])
        Matrix[:,i]=np.sort(Erwartete_Eigenwerte)
    return Matrix


def Eigenwerte_Matrix(Potential,Frequenz,Anzahl):
    print('test',Anzahl)
    Matrix =np.genfromtxt('build/Eigenwerte_fur_a='+str(int(Potential))+'_w='+str(int(Frequenz))+'_N='+str(int(Anzahl))+'.txt')
#        Matrix=np.zeros((np.size(Grosse),np.size(Frequenz_1000)))
    return Matrix

# Eigenwerte_H_0=np.array([-2,0,0,2])
Figure_Zahler=1

for index_Potential , value_Potential in enumerate(tqdm(Potential)):
    for index_Frequenz , value_Frequenz in enumerate(tqdm(Frequenz)):
        for index_N ,value_N in enumerate(tqdm(Anzahl_N)):
            print('N',value_N)
            plt.figure(Figure_Zahler)
            Figure_Zahler=Figure_Zahler+1
            x=np.linspace(np.amin(Energien)/10000,np.amax(Energien)/10000,100)
            Matrix_mit_Eigenwerten=Eigenwerte_Matrix(value_Potential,value_Frequenz*10000,value_N)
#           Matrix_mit_Erwarteteneigenwerten=eig_erwartet_Matrix(Anzahl,x,value_Potential)
            n, m=np.shape(Matrix_mit_Eigenwerten)
            plt.title('Eigenwerte von a='+str(value_Potential/10000)+' w='+str(value_Frequenz)+' N='+str(value_N))
            for j in tqdm(range(n-1)):
                plt.plot(Energien/10000,Matrix_mit_Eigenwerten[j,:],'-b',alpha=0.5)#,label='Eigenwert'+ str(j) )
                #plt.plot(Frequenz,Matrix_mit_Eigenwerten[j,:],'xr',alpha=0.5)#,label='Eigenwert'+ str(j) )
                #print(Frequenz)
            plt.plot(Energien/10000,Matrix_mit_Eigenwerten[n-1,:],'-b',alpha=0.5,label=r'$\epsilon_\alpha$ bei $\omega$='+ str(value_Frequenz) )
            plt.plot(x,x*0+value_Frequenz/2,'--k',alpha=0.5,linewidth=0.5)
            plt.plot(x, x*0-value_Frequenz/2,'--k',alpha=0.5,linewidth=0.5)#,label='Eigenwert'+ str(j) )
            plt.plot(x, x*0+value_Frequenz*3/2,'--k',alpha=0.5,linewidth=0.5)#,label='Eigenwert'+ str(j) )
            plt.plot(x, x*0-value_Frequenz*3/2,'--k',alpha=0.5,linewidth=0.5)#,label='Eigenwert'+ str(j) )
            plt.yticks( [-value_Frequenz*3/2, -value_Frequenz/2,0,  value_Frequenz/2, value_Frequenz*3/2],
               [r'$-\frac{3\omega}{2}$', r'$-\frac{\omega}{2}$',0,  r'$\frac{\omega}{2}$', r'$\frac{3\omega}{2}$'])
            plt.ylim(-value_Frequenz*2,value_Frequenz*2)
            plt.legend(loc='upper right')
            plt.xlabel(r'$E_0$')
            plt.ylabel(r'$\epsilon_\alpha $')
            plt.savefig('Plots/Plot_fur'+'_a='+str(Potential[index_Potential]/10000)+'_w='+str(value_Frequenz)+'N='+str(value_N) +'.pdf')
            plt.close()



        #     Erwartete_Eigenwerte=eig_war(Eigenwerte_H_0,Anzahl,Frequenz,value_Potential)
        #     n=np.size(Eigenwerte)
        #     for i in range(n-1):
        #         plt.plot(Frequenz,Eigenwerte[i],'xb')
        #         plt.plot(Frequenz,Erwartete_Eigenwerte[i],'+r')
        #     plt.plot(Frequenz,Eigenwerte[n-1],'xb',label=r'E='+str(value_Energie/100))
        #     plt.plot(Frequenz,Erwartete_Eigenwerte[n-1],'+r',label=r'Erwartete_Eigenwerte')
        #     print('Potential',value_Potential/100)
        #     print('value_Energie',value_Energie/100)
        #     print('Frequenz=', Frequenz)
        #     print('Erwartete_Eigenwerte',np.sort(Erwartete_Eigenwerte))
        #     print('Berechnete Eigenwerte',Eigenwerte)
