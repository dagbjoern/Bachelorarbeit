import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from tqdm import tqdm


def norm(v):
    n=v/np.sqrt(np.dot(v,v))
    return n

def test(phi, a):
    vektor = np.zeros([4,1])+1j*0
    for i in range(np.size(phi[0,:])):
        vektor[i,0]= np.dot(norm(phi[:,a-1]),norm(phi[:,i]))
        print(vektor)
    return vektor

def skalarprodukt(v1,v2):
    ans=np.dot(np.conjugate(v1),v2)
    return ans


def konstanten(Startzustand, phi, epsilon):
    c = np.zeros(np.size(epsilon))+1j*0
    for i in range(np.size(c)):
            c[i] = skalarprodukt(phi[:, i], Startzustand)
    return c

def zeitentwicklung(Startzustand,V_phi,epsilon,frequenz,t):
    psi=np.zeros([np.size(Startzustand),np.size(t)])+1j*0
    phi_0=phi_funktion(V_phi,0,frequenz)
    for index_Zeit , value_Zeit in enumerate(tqdm(t)):
    #    print('fortschritt:',value_Zeit,'/',np.amax(t))
        phi_t=phi_funktion(V_phi,value_Zeit,frequenz)
        #print('t=',value_Zeit)
        for index_epsilon, value_epsilon in enumerate(epsilon):
            #print('c_',index_epsilon,'=',  skalarprodukt(phi_0[:,index_epsilon],Startzustand))
            c=konstanten(Startzustand,phi_0,epsilon)
#            print('Startzustand=',Startzustand)
#            print('Startzustand Berechnete',phi_0[:,0]*c[0]+phi_0[:,1]*c[1]+phi_0[:,2]*c[2]+phi_0[:,3]*c[3])
            #print('c=',c)

            psi[:,index_Zeit]=psi[:,index_Zeit]+phi_t[:,index_epsilon]*np.exp(-1j*value_epsilon*value_Zeit)*skalarprodukt(phi_0[:,index_epsilon],Startzustand)
#        print('psi',psi)
    return psi




def betragsquadrad(messwerte,psi_t):
    for index , value_Zeit in enumerate(t):
        if index==0:
            werte=np.abs(np.dot(messwerte,psi_t[:,index]))**2
        else:
            werte=np.append(werte,np.abs(np.dot(messwerte,psi_t[:,index]))**2)
    return werte

def phi_funktion(V,t,w):
    phi=np.zeros([4,4])
    phi=phi+1j*0
    for q in range(np.size(V[0,:])):
        n=-(np.size(V[:,0])/4-1)/2
#        print(np.size(V[:,0])/4)
        for r in range(int(np.size(V[:,0])/4)):
            for s in range(4):
    #                print('n',n,'\n s',s,'\n')
                    phi[s,q] = phi[s,q] + V[r*4+s,q]*np.exp(1j*w*t)
            n=n+1
    return phi



Energien=np.genfromtxt('build/Durchlaufende_Energien.txt')
Potential=np.genfromtxt('build/Durchlaufende_Potentiale.txt')
Frequenz_1000=np.genfromtxt('build/Durchlaufende_Frequenzen.txt',unpack=True)
Anzahl_N=np.genfromtxt('build/Durchlaufende_N.txt',unpack=True)
Gitterkonstante, Anzahl, Phasenverschiebung = np.genfromtxt('Parameter/Parameter_fur_a='+str(int(Potential))+'_w='+ str(int(Frequenz_1000))+'_E='+str(int(Energien))+'_N='+str(int(Anzahl_N[0]))+'.txt',unpack=True)
H_0_eigenvektoren=np.genfromtxt('Parameter/eigenvektoren_von_H_0_fur_a='+str(int(Potential))+'.txt')
H_0_Eigenwerte=np.genfromtxt('Parameter/eigenwerte_von_H_0_fur_a='+str(int(Potential)) +'.txt')


Figure_Zahler=1

for l in tqdm(range(np.size(Anzahl_N))):
    Eigenwerte=np.genfromtxt('build/Eigenwerte_fur_a='+str(int(Potential))+'_w='+ str(int(Frequenz_1000))+'_E='+str(int(Energien))+'_N='+str(int(Anzahl_N[0]))+'.txt',unpack=True)

    #print(Eigenwerte)
    #print(np.size(Eigenwerte)/2)
    epsilon=np.genfromtxt('build/epsilon_fur_a='+str(int(Potential))+'_E='+str(int(Energien))+'_w=' + str(int(Frequenz_1000)) +'a_N=' + str(int(Anzahl_N[l])) +'.txt')

    #print(epsilon)
    Frequenz=Frequenz_1000/1000
    Eigenzustande_realteil=np.genfromtxt('build/Realpart_Eigenzustande_fur_a='+str(int(Potential))+'_E='+str(int(Energien))+'_w=' + str(int(Frequenz_1000)) +'a_N=' + str(int(Anzahl_N[l])) +'.txt')
    Eigenzustande_imagteil=np.genfromtxt('build/Imagpart_Eigenzustande_fur_a='+str(int(Potential))+'_E='+str(int(Energien))+'_w=' + str(int(Frequenz_1000)) +'a_N=' + str(int(Anzahl_N[l])) +'.txt')
    #print(Eigenvektoren_realteil,Eigenvektoren_imagteil)
    Eigenzustande_realteil=Eigenzustande_realteil+1j*0
    Eigenzustande_imagteil=Eigenzustande_imagteil*1j
    V_phi=Eigenzustande_realteil+Eigenzustande_imagteil
    phi=phi_funktion(V_phi,0,Frequenz)
    plt.figure(Figure_Zahler)
    t=np.linspace(0,30,500)
    Startzustand=H_0_eigenvektoren[:,0]
    psi_t=zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t)

    Figure_Zahler=1+Figure_Zahler
    plt.title('zeitentwicklung für den Startzustand ' + str(phi[:,0]))
    plt.plot(t,np.real(psi_t[2]),'-b',label=r'3 real')
    plt.plot(t,np.imag(psi_t[2]),'-b',alpha=0.3,label=r'3 imag')
    plt.plot(t,np.real(psi_t[3]),'-r',label=r'4 real')
    plt.plot(t,np.imag(psi_t[3]),'-r',alpha=0.3,label=r'4 imag')
    plt.plot(t,np.real(psi_t[0]),'-g',label=r'1 real')
    plt.plot(t,np.imag(psi_t[0]),'-g',alpha=0.3,label=r'1 imag')
    plt.plot(t,np.real(psi_t[1]),'-y',label=r'2 real')
    plt.plot(t,np.imag(psi_t[1]),'-y',alpha=0.3,label=r'2 imag')
    #
    # #
    # plt.plot(t,zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t),'-g',label=r'2')
    #
    # plt.plot(t,zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t),'-g',label=r'3')
    # plt.plot(t,zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t)[0],'-g',label=r'4')
    #
    plt.legend(loc='best')
    plt.savefig('Plots/zeitentwicklung phi_1 N='+str(int(Anzahl_N[l]))+ 'a='+str(Potential/1000)+'.pdf')





#    print(betragsquadrad(psi_1,psi_1,phi,epsilon,t))
    print(Startzustand)
    t=np.linspace(0,100,1000)
    T = np.linspace(0, 1, 2)
    psi_t=zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t)
    Periodendauer = 2 * np.pi / Frequenz

    phi=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    plt.figure(Figure_Zahler )
    Figure_Zahler=1+Figure_Zahler

    #betragsquadrad(messwerte,Startzustand,V_phi,epsilon,frequenz,t):
    plt.title('zeitentwicklung für den Startzustand ' + str(np.round(Startzustand,3))+ '\n w=' + str(Frequenz) + 'E='+str(Energien/100) + 'a='+str(Potential/100)  )
    plt.plot(t,betragsquadrad(phi[:,0],psi_t),'-r',alpha=0.5,label=r'zustand 1')
    plt.plot(t,betragsquadrad(phi[:,1],psi_t),'-y',alpha=0.5,label=r'zustand 2')
    plt.plot(t,betragsquadrad(phi[:,2],psi_t),'-g',alpha=0.5,label=r'zustand 3')
    plt.plot(t,betragsquadrad(phi[:,3],psi_t),'-b',alpha=0.5,label=r'zustand 4')
    for n in range(4):
        plt.plot(T * 0 + Periodendauer * n, T, '--k',linewidth=0.5)
    plt.legend(loc='best')
    plt.savefig('Plots/zeitentwicklung_der_Besetzungen_N='+str(int(Anzahl_N[l]))+ 'a='+str(Potential/1000)+'.pdf')


    plt.figure(Figure_Zahler )
    Figure_Zahler=1+Figure_Zahler
    plt.title('zeitentwicklung für den Startzustand ' + str(np.round(Startzustand,3))+ '\n w=' + str(Frequenz) + 'E='+str(Energien/100) + 'a='+str(Potential/100)  )
    plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 0], psi_t),'-r',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[0],3)))
    plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 1], psi_t),'-y',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[1],3)))
    plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 2], psi_t),'-g',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[2],3)))
    plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 3], psi_t),'-b',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[3],3)))
    for n in range(4):
        plt.plot(T * 0 + Periodendauer * n, T, '--k')
    plt.legend(loc='best')
    plt.savefig('Plots/zeitentwicklung_der_Eigenzustande_N='+str(int(Anzahl_N[l]))+ 'a='+str(Potential/1000)+'.pdf')
#
#   print(betragsquadrad(psi_1,psi_1,phi,epsilon,t)+betragsquadrad(psi_2,psi_1,phi,epsilon,t)+betragsquadrad(psi_3,psi_1,phi,epsilon,t)+betragsquadrad(psi_4,psi_1,phi,epsilon,t))

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