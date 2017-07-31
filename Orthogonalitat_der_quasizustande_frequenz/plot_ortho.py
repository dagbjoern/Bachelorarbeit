import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from tqdm import tqdm
import os





def skalarprodukt(v1,v2):
    ans=np.dot(np.conjugate(v1),v2)
    return ans



def phi_funktion(V,t,w):
    phi=np.zeros([4,4])
    phi=phi+1j*0
    for q in range(np.size(V[0,:])):
        n=-(np.size(V[:,0])/4-1)/2
        for r in range(int(np.size(V[:,0])/4)):
            for s in range(4):
                    phi[q,s] = phi[q,s] + V[r*4+s,q]*np.exp(1j*w*t)
            n=n+1
    return phi




Energien=np.genfromtxt('build/Durchlaufende_Energien.txt')
Potential=np.genfromtxt('build/Durchlaufende_Potentiale.txt')
Frequenz_10000=np.genfromtxt('build/Durchlaufende_Frequenzen.txt',unpack=True)
Anzahl_N=np.genfromtxt('build/Durchlaufende_N.txt',unpack=True)


str_Potential = Potential.astype(int)
str_Potential = str_Potential.astype(str)

str_Frequenz_10000 = Frequenz_10000.astype(int)
str_Frequenz_10000 = str_Frequenz_10000.astype(str)

str_Energien = Energien.astype(int)
str_Energien = str_Energien.astype(str)

str_Anzahl_N = Anzahl_N.astype(int)
str_Anzahl_N = str_Anzahl_N.astype(str)

Gitterkonstante, Anzahl, Phasenverschiebung = np.genfromtxt('Parameter/Parameter_fur_a='+str_Potential[0]+'_w='+ str_Frequenz_10000[0]+'_E='+str_Energien[0]+'_N='+str_Anzahl_N[0]+'.txt',unpack=True)

Frequenz = Frequenz_10000/10000

Figure_Zahler = 1
for a in tqdm(range(int(np.size(Potential)))):
    if not os.path.exists('Plots/Potential='+ str(Potential[a]/10000)):
        os.makedirs('Plots/Potential='+ str(Potential[a]/10000))
    for e in tqdm(range(int(np.size(Energien)))):
        if not os.path.exists('Plots/Potential='+ str(Potential[a]/10000) + '/Energie='+str(Energien[e]/10000)):
            os.makedirs('Plots/Potential='+ str(Potential[a]/10000)+ '/Energie='+str(Energien[e]/10000))
        for l in tqdm(range(np.size(Anzahl_N))):
            if not os.path.exists('Plots/Potential='+ str(Potential[a]/10000) + '/Energie='+str(Energien[e]/10000)+'/Anzahl='+ str_Anzahl_N[l]):
                os.makedirs('Plots/Potential='+ str(Potential[a]/10000)+ '/Energie='+str(Energien[e]/10000)+'/Anzahl=' + str_Anzahl_N[l])
            werte_1=np.zeros([4,np.size(Frequenz)])+1j*0
            werte_2=np.zeros([4,np.size(Frequenz)])+1j*0
            werte_3=np.zeros([4,np.size(Frequenz)])+1j*0
            werte_4=np.zeros([4,np.size(Frequenz)])+1j*0
            for f in tqdm(range(int(np.size(Frequenz)))):
                #Eigenwerte=np.genfromtxt('build/Eigenwerte_fur_a='+str_Potential[a]+'_w='+ str_Frequenz_1000[f]+'_E='+str_Energien[e]+'_N='+str_Anzahl_N[l]+'.txt',unpack=True)
                epsilon=np.genfromtxt('build/epsilon_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'a_N=' + str_Anzahl_N[l] +'.txt')
                Eigenzustande_realteil=np.genfromtxt('build/Realpart_Eigenzustande_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'a_N=' + str_Anzahl_N[l] +'.txt')
                Eigenzustande_imagteil=np.genfromtxt('build/Imagpart_Eigenzustande_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'a_N=' + str_Anzahl_N[l] +'.txt')
                Eigenzustande_realteil=Eigenzustande_realteil+1j*0
                Eigenzustande_imagteil=Eigenzustande_imagteil*1j
                V_phi=Eigenzustande_realteil+Eigenzustande_imagteil
                #print(Potential[a]/100, Energien[e]/10000)
                phi=phi_funktion(V_phi,0,Frequenz[f])
                #print(phi)
                for q in range(4):
                    werte_1[q,f]=skalarprodukt(phi[:,0],phi[:,q])
                    werte_2[q,f]=skalarprodukt(phi[:,1],phi[:,q])
                    werte_3[q,f]=skalarprodukt(phi[:,2],phi[:,q])
                    werte_4[q,f]=skalarprodukt(phi[:,3],phi[:,q])
                #print(werte_4[0,:])
            plt.figure(Figure_Zahler)
            Figure_Zahler=1+Figure_Zahler
            plt.title('Orthogonalität für n=' + str(Anzahl_N[l]) + '\n  E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/10000) )
            plt.plot(Frequenz, werte_1[0,:].real,'-y',alpha=0.5,label=r'phi1*phi1' )
            plt.plot(Frequenz, werte_1[1,:].real,'-r',alpha=0.5,label=r'phi1*phi2' )
            plt.plot(Frequenz, werte_1[2,:].real,'-g',alpha=0.5,label=r'phi1*phi3' )
            plt.plot(Frequenz, werte_1[3,:].real,'-b',alpha=0.5,label=r'phi1*phi4' )
            # plt.ylim(-Frequenz/2,Frequenz/2)
            plt.legend(loc='best')
            plt.xlabel('Frequenz')
            plt.ylabel(r'$phi*phi$')
            plt.savefig('Plots/Potential='+ str(Potential[a]/10000)+ '/Energie='+str(Energien[e]/10000) +'/Anzahl='+ str_Anzahl_N[l]+'/Plot_für_phi1_phi_i.pdf')
            plt.close()


            # plt.figure(Figure_Zahler)
            # Figure_Zahler=1+Figure_Zahler
            # plt.title('Orthogonalität für n=' + str(Anzahl_N[l]) + '\n w=' + str(Frequenz[f]) + ' E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/100) )
            # plt.plot(Anzahl_N, werte_2[0,:].real,'-y',alpha=0.5,label=r'phi2*phi1' )
            # plt.plot(Anzahl_N, werte_2[1,:].real,'-r',alpha=0.5,label=r'phi2*phi2' )
            # plt.plot(Anzahl_N, werte_2[2,:].real,'-g',alpha=0.5,label=r'phi2*phi3' )
            # plt.plot(Anzahl_N, werte_2[3,:].real,'-b',alpha=0.5,label=r'phi2*phi4' )
            # # plt.ylim(-Frequenz/2,Frequenz/2)
            # plt.legend(loc='best')
            # plt.xlabel('N')
            # plt.ylabel(r'$phi*phi$')
            # plt.savefig('Plots/Potential='+ str(Potential[a]/100)+ '/Energie='+str(Energien[e]/10000) +'/Frequenz='+ str(Frequenz[f])+'/Plot_für_phi2_phi_i.pdf')
            # plt.close()
            #
            # plt.figure(Figure_Zahler)
            # Figure_Zahler=1+Figure_Zahler
            # plt.title('Orthogonalität für n=' + str(Anzahl_N[l]) + '\n w=' + str(Frequenz[f]) + ' E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/100) )
            # plt.plot(Anzahl_N, werte_3[0,:].real,'-y',alpha=0.5,label=r'phi3*phi1' )
            # plt.plot(Anzahl_N, werte_3[1,:].real,'-r',alpha=0.5,label=r'phi3*phi2' )
            # plt.plot(Anzahl_N, werte_3[2,:].real,'-g',alpha=0.5,label=r'phi3*phi3' )
            # plt.plot(Anzahl_N, werte_3[3,:].real,'-b',alpha=0.5,label=r'phi3*phi4' )
            # # plt.ylim(-Frequenz/2,Frequenz/2)
            # plt.legend(loc='best')
            # plt.xlabel('N')
            # plt.ylabel(r'$phi*phi$')
            # plt.savefig('Plots/Potential='+ str(Potential[a]/100)+ '/Energie='+str(Energien[e]/10000) +'/Frequenz='+ str(Frequenz[f])+'/Plot_für_phi3_phi_i.pdf')
            # plt.close()
            #
            # plt.figure(Figure_Zahler)
            # Figure_Zahler=1+Figure_Zahler
            # plt.title('Orthogonalität für n=' + str(Anzahl_N[l]) + '\n w=' + str(Frequenz[f]) + ' E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/100) )
            # plt.plot(Anzahl_N, werte_4[0,:].real,'-y',alpha=0.5,label=r'phi4*phi1' )
            # plt.plot(Anzahl_N, werte_4[1,:].real,'-r',alpha=0.5,label=r'phi4*phi2' )
            # plt.plot(Anzahl_N, werte_4[2,:].real,'-g',alpha=0.5,label=r'phi4*phi3' )
            # plt.plot(Anzahl_N, werte_4[3,:].real,'-b',alpha=0.5,label=r'phi4*phi4' )
            # # plt.ylim(-Frequenz/2,Frequenz/2)
            # plt.legend(loc='best')
            # plt.xlabel('N')
            # plt.ylabel(r'$phi*phi$')
            # plt.savefig('Plots/Potential='+ str(Potential[a]/100)+ '/Energie='+str(Energien[e]/10000) +'/Frequenz='+ str(Frequenz[f])+'/Plot_für_phi4_phi_i.pdf')
            # plt.close()


#
#print(np.dot(phi[:,0],phi[:,1]))

#   if [(np.size(Frequenz_1000)==1)and(np.size(Potential)==1)]:
#   r i in range(np.size(phi[0,:])):
#     test(phi,i+1)


#print(Eigenvektoren_realteil,Eigenvektoren_imagteil)
#print(Eigenvektoren)






# def eig_erwartet_Matrix(Anzahl,Frequenz,Potential):
#     Eigenwerte_H_0=np.genfromtxt('Parameter/eigenwerte_von_H_0_fur_a='+str(int(Potential)) +'a.txt')
#     m=np.size(Eigenwerte_H_0)*(1+Anzahl*2)
#     Matrix=np.zeros((m,np.size(Frequenz)))
#     Frequenz=Frequenz*Potential/100
#     for i in range(np.size(Frequenz)):
#         Erwartete_Eigenwerte=Eigenwerte_H_0
#         for j in range(1,int(Anzahl+1)):
#             Erwartete_Eigenwerte=np.append(Erwartete_Eigenwerte,Eigenwerte_H_0-Frequenz[i])
#             Erwartete_Eigenwerte=np.append(Erwartete_Eigenwerte,Eigenwerte_H_0+Frequenz[i])
#         Matrix[:,i]=np.sort(Erwartete_Eigenwerte)
#     return Matrix

#
# def Eigenwerte_Matrix(Potential,Frequenz):
#     Matrix =np.genfromtxt('build/Eigenwerte_fur_a='+str(int(Potential))+'_w='+str(int(Frequenz))+'.txt')
# #        Matrix=np.zeros((np.size(Grosse),np.size(Frequenz_1000)))
#     return Matrix

# Eigenwerte_H_0=np.array([-2,0,0,2])
# Figure_Zahler=1
#
# #Figure_Zahler=Figure_Zahler+1
# plt.title('Eigenwerte von a='+str(Potential/100)+' E='+str(Frequenz) )
# plt.figure(Figure_Zahler)
# for index_N , value_N in enumerate(Anzahl_N):
#         Eigenwerte = np.genfromtxt('build/Eigenwerte_fur_a='+str(int(Potential))+'_w='+str(int(Frequenz_1000))+'_E='+str(int(Energien))+'_N='+str(int(value_N))+'.txt')
#         plt.plot(value_N+Eigenwerte*0,Eigenwerte,'xb',alpha=0.5)
#         print('N=',value_N,)
#         for index_e , value_e in enumerate(Eigenwerte):
#             if (value_e<(Frequenz/2) and (-Frequenz/2) < value_e):
#                 print(value_e)
#         # plt.plot(x, x*0-value_Frequenz/2,'--k',alpha=0.5)#,label='Eigenwert'+ str(j) )
#         # plt.plot(x, x*0-value_Frequenz/2,'--k',alpha=0.5)#,label='Eigenwert'+ str(j) )
#         # plt.plot(x, x*0+value_Frequenz*3/2,'--k',alpha=0.5)#,label='Eigenwert'+ str(j) )
#         # plt.plot(x, x*0-value_Frequenz*3/2,'--k',alpha=0.5)#,label='Eigenwert'+ str(j) )
# plt.ylim(-Frequenz/2,Frequenz/2)
# plt.legend(loc='best')
# plt.xlabel('N')
# plt.ylabel(r'$\epsilon_\alpha $')
# plt.savefig('Plots/Plot_fur'+'_a='+str(Potential/100)+'_w='+str(Frequenz) +'_E='+str(int(Energien))+'.pdf')
# plt.close()
#
        # for j in range(n):
        #     r=j%2
        #     b=(j+1)%2
        #
        #     plt.plot(Energien/100,Matrix_mit_Eigenwerten[j,:],'-b',color=(r,0,b),alpha=0.5)#,label='Eigenwert'+ str(j) )
            #plt.plot(Frequenz,Matrix_mit_Eigenwerten[j,:],'xr',alpha=0.5)#,label='Eigenwert'+ str(j) )



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
