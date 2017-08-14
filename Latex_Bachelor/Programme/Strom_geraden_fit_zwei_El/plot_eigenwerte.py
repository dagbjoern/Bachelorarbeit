import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from tqdm import tqdm
from scipy.optimize import curve_fit
import os
from scipy.optimize import curve_fit

def norm(v):
    n=v/np.sqrt(np.dot(v,v))
    return n


def  Strom(c):
    J = np.matrix([[0, -c, 0, c, ], [c, 0, -c, 0], [0, c, 0, -c], [-c, 0, c, 0]])
    return 1j/4*J

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
            # print('Startzustand=',Startzustand)
#            print('Startzustand Berechnete',phi_0[:,0]*c[0]+phi_0[:,1]*c[1]+phi_0[:,2]*c[2]+phi_0[:,3]*c[3])
            #print('c=',c)

            psi[:,index_Zeit]=psi[:,index_Zeit]+phi_t[:,index_epsilon]*np.exp(-1j*value_epsilon*value_Zeit)*skalarprodukt(phi_0[:,index_epsilon],Startzustand)
#        print('psi',psi)
    return psi


def Stromittelwert_Foquet(Startzustand, c, V_phi):
        phi_0 = phi_funktion(V_phi, 0, 0)
        J = Strom(c)
        Mittelwert = 0+1j*0
        c_a = np.array([0.0, 0.0, 0.0, 0.0])
        hilfsarray = np.zeros([4, int(np.size(V_phi[:,0])/4)])  #
        for i in range(np.size(phi_0[0, :])):
#            print('test',np.abs(np.dot(Startzustand, phi_0[:, i]))**2)
            c_a[i] =  np.abs(np.dot(Startzustand, phi_0[:, i]))**2  # konstanten berechnen
#            print(c_a[i])
            for n in range(int(np.size(hilfsarray[0, :]))):
                obergrenze=int(4*(1+n))
                untergrenze=int(0+(4*n))
                # print('max',np.size(hilfsarray[0, :]))
                # print(obergrenze,untergrenze)
                # print(V_phi[untergrenze:obergrenze,i])
                # print('skalarprodukt',np.dot(np.dot(J,V_phi[untergrenze:obergrenze,i]),np.conjugate(V_phi[untergrenze:obergrenze,i]))+1j*0)
                # print(V_phi[untergrenze:obergrenze,i])
                wert = np.dot(np.dot(J,V_phi[untergrenze:obergrenze,i]),np.conjugate(V_phi[untergrenze:obergrenze,i]))+1j*0
                hilfsarray[i, n]=np.amax(wert)
    #            print(0+(4*n), 4*(1+n))
#            print('c_a',c_a)
            Mittelwert = Mittelwert+c_a[i]*np.sum(hilfsarray[i, :])
        return Mittelwert


def betragsquadrad(messwerte,psi_t):
    for index , value_Zeit in enumerate(t):
        if index==0:
            werte=np.abs(np.dot(messwerte,psi_t[:,index]))**2
        else:
            werte=np.append(werte,np.abs(np.dot(messwerte,psi_t[:,index]))**2)
    return werte

def Erwartungswert(Operator,psi_t,t):
    for index , value_Zeit in enumerate(t):
        if index==0:
            werte=np.dot(np.dot(Operator,psi_t[:,index]),np.conjugate(psi_t[:,index]))+1j*0
        else:
            werte=np.append(werte,np.dot(np.dot(Operator,psi_t[:,index]),np.conjugate(psi_t[:,index])),axis=0)
    return werte

def gerade(x,m):
    return m*x


def phi_funktion(V,t,w):
    phi=np.zeros([4,4])
    phi=phi+1j*0
    N=(np.size(V[:,0])/4-1)/2
    n=np.linspace(-N,N,2*N+1)
#    print(n)
    for q in range(np.size(V[0,:])):
#        print(np.size(V[:,0])/4)
        for r in range(int(np.size(V[:,0])/4)):
            for s in range(4):
    #                print('n',n,'\n s',s,'\n')
                    phi[s,q] = phi[s,q] + V[r*4+s,q]*np.exp(1j*n[r]*w*t)
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

str_Anzahl_N = str(int(Anzahl_N))

Figure_Zahler = 1
a=np.linspace(0,np.amax(Potential)/10000,10000)
Frequenz = Frequenz_10000/10000

if not os.path.exists('Plots_mittelwerte/'):
    os.makedirs('Plots_mittelwerte/')


plt.figure(Figure_Zahler)
Figure_Zahler=Figure_Zahler+1
#plt.title('zeitentwicklung für den Startzustand ' + str(np.round(Startzustand,3))+'fur n=' + str(Anzahl_N[l]) + '\n w=' + str(Frequenz[f]) + ' E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/10000) )
plt.plot(a,np.sqrt(a**2+4*1**2)-a , '-r', alpha=0.75, label=r'erste Resonanz')
plt.plot(a,np.sqrt(a**2+4*1**2)+a , '-r', alpha=0.75, label=r'zweite Resonanz')
plt.plot(a,2*np.sqrt(a**2+4*1**2) , '-r', alpha=0.75, label=r'dritte Resonanz')
for i in range(np.size(Frequenz)):
    plt.plot(a,a*0+Frequenz[i],'--b')
plt.xlabel(r'Potential a')
plt.ylabel(r'Frequenz w')
#plt.xlim(np.amin(t),np.amax(t))
plt.legend(loc='best')
plt.savefig('Plots_mittelwerte/Resonanzen.pdf')
plt.close()




#Gitterkonstante, Anzahl, Phasenverschiebung = np.genfromtxt('Parameter/Parameter_fur_a='+str_Potential[0]+'_w='+ str_Frequenz_100000[0]+'_E='+str_Energien[0]+'_N='+str_Anzahl_N[0]+'.txt',unpack=True)

print('int',Anzahl_N )
print('str', str_Anzahl_N)
#print('cool' + str_Anzahl_N[0] + 'cool')
t_lsode=np.genfromtxt('build/Zeit_lsode.txt')
grenze=np.array([1,0.8,0.6])
I_bar=Energien*0 #array der länge von E erschaffen

farbe=np.array(['r','b','g','m'])


for a in tqdm(range(np.size(Potential))):
    H_0_eigenvektoren=np.genfromtxt('Parameter/eigenvektoren_von_H_0_fur_a='+str_Potential[a]+'.txt')
    H_0_Eigenwerte=np.genfromtxt('Parameter/eigenwerte_von_H_0_fur_a='+str_Potential[a] +'.txt')
    if not os.path.exists('Plots/Potential='+ str(Potential[a]/10000)):
        os.makedirs('Plots/Potential='+ str(Potential[a]/10000))
    for e in tqdm(range(np.size(Energien))):
        I_bar_1=Frequenz*0+1j*0 #array der länge von E erschaffen
        I_bar_2=Frequenz*0+1j*0 #array der länge von E erschaffen
        I_bar_lsode=Frequenz*0+1j*0 #array der länge von E erschaffen
        #               Numerische Werte Einlesen
        if not os.path.exists('Plots/Potential='+ str(Potential[a]/10000) + '/Energie='+str(Energien[e]/10000)):
            os.makedirs('Plots/Potential='+ str(Potential[a]/10000)+ '/Energie='+str(Energien[e]/10000))
        for f in tqdm(range(np.size(Frequenz))):
            f=f
            psi_real_lsode=np.genfromtxt('build/Realpart_Eigenvektoren_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'lsode.txt')
            psi_imag_lsode=np.genfromtxt('build/Imagpart_Eigenvektoren_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'lsode.txt')
            psi_imag_lsode=psi_imag_lsode*1j
            psi_t_lsode=psi_real_lsode+psi_imag_lsode
            t=t_lsode
            psi_t_lsode=psi_t_lsode.T
#            for l in tqdm(range(np.size(Anzahl_N))):

            #Eigenwerte=np.genfromtxt('build/Eigenwerte_fur_a='+str_Potential[a]+'_w='+ str_Frequenz_10000[f]+'_E='+str_Energien[e]+'_N='+str_Anzahl_N[l]+'.txt',unpack=True)
            epsilon=np.genfromtxt('build/epsilon_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'_N=' + str_Anzahl_N +'.txt')
            Eigenzustande_realteil=np.genfromtxt('build/Realpart_Eigenzustande_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'_N=' + str_Anzahl_N +'.txt')
            Eigenzustande_imagteil=np.genfromtxt('build/Imagpart_Eigenzustande_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_10000[f] +'_N=' + str_Anzahl_N +'.txt')
            Eigenzustande_realteil=Eigenzustande_realteil+1j*0
            Eigenzustande_imagteil=Eigenzustande_imagteil*1j
    #        Startzustand=H_0_eigenvektoren[:,0]
            V_phi=Eigenzustande_realteil+Eigenzustande_imagteil
            phi=phi_funktion(V_phi,0,Frequenz[f])
            # werte_1=np.zeros([1,4])+1j*0
            # werte_2=np.zeros([1,4])+1j*0
            # werte_3=np.zeros([1,4])+1j*0
            # werte_4=np.zeros([1,4])+1j*0
            # for q in range(4):
            #     werte_1[0,q]=skalarprodukt(phi[:,0],phi[:,q])
            #     werte_2[0,q]=skalarprodukt(phi[:,1],phi[:,q])
            #     werte_3[0,q]=skalarprodukt(phi[:,2],phi[:,q])
            #     werte_4[0,q]=skalarprodukt(phi[:,3],phi[:,q])
            #     # if(np.sum(werte_1) < 0.99 or np.sum(werte_1) > 1.01 or np.sum(werte_2) < 0.99 or np.sum(werte_2) > 1.01 or np.sum(werte_3) < 0.99 or np.sum(werte_3) > 1.01 or np.sum(werte_4) < 0.99 or np.sum(werte_4) > 1.01):
                #     print('nicht orthogonal genug')
                #     print( 'für:\n n=' + str(Anzahl_N[l]) + '\n w=' + str(Frequenz[f]) + '\n E=' +str(Energien[e]/10000) + '\n a=' +str(Potential[a]/100)   )
                #     print('Werte_1:',werte_1)
                #     print('Werte_1sum:',np.sum(werte_1))
                #     print('Werte_2:',werte_2)
                #     print('Werte_2sum:',np.sum(werte_2))
                #     print('Werte_3:',werte_3)
                #     print('Werte_3sum:',np.sum(werte_3))
                #     print('Werte_4:',werte_4)
                #     print('Werte_4sum:',np.sum(werte_4))
                #     quit()
                 #print('Strommittelwert für Potential='+ str(Potential[a]/10000)+'\nEnergie='+str(Energien[e]/10000)+'\nw = ' + str(Frequenz[f])  ,'Strommittel=',Stromittelwert_Foquet(Startzustand, 1, V_phi))
            I_bar_1[f]=Stromittelwert_Foquet(H_0_eigenvektoren[:,0], 1, V_phi)
            I_bar_2[f]=Stromittelwert_Foquet(H_0_eigenvektoren[:,1], 1, V_phi)
#            print(I_bar[f])
            y_lsode=Erwartungswert(Strom(1),psi_t_lsode,t_lsode)
                #print('numerischer Wert:', np.sum(y_lsode)/np.size(t_lsode))
            I_bar_lsode[f]=np.sum(y_lsode)/np.size(t_lsode)
                #Periodendauer = 2 * np.pi / Frequenz[f]
                #phi=phi_funktion(V_phi,0,Frequenz[f])
#                test1=phi_funktion(V_phi,0,Frequenz[f])
                #print('t=0',test)
#                test2=phi_funktion(V_phi,Periodendauer,Frequenz[f])
                #print('test',test1-test2)
#                t=t_lsode
#                psi_t=zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz[f],t)
                #print(Erwartungswert(Strom(1),psi_t))
            plt.figure(Figure_Zahler)
                # y=Erwartungswert(Strom(1),psi_t,t)
                #betragsquadrad(messwerte,Startzustand,V_phi,epsilon,frequenz,t):
                #plt.title('Erwartungswert des Stromes für den Startzustand ' + str(np.round(Startzustand,3))+'\n fur n=' + str(Anzahl_N[l]) + '\n w=' + str(Frequenz[f]) + ' E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/10000) )
                #plt.plot(t,y, '-r', alpha=0.5, label=r'Strom')
                #plt.plot(t,y_lsod, '-r', alpha=0.5, label=r'Strom')
                # plt.plot(t, Erwartungswert(Strom(2),psi_t), '-y', alpha=0.5, label=r'c=2')
                # plt.plot(t, Erwartungswert(Strom(3),psi_t), '-b', alpha=0.5, label=r'c=3')
                # plt.plot(t, Erwartungswert(Strom(4),psi_t), '-g', alpha=0.5, label=r'c=4')
                #T = np.linspace(np.amin(y),np.amax(y) , 2)
                #for n in range(10):
                #    plt.plot(T * 0 + Periodendauer * n, T, '--k',linewidth=0.5)
                #plt.xlabel(r'Zeit $t/ j^{-1}$')
                #plt.ylabel(r'Strom $I/c $')
                #plt.xlim(np.amin(t),np.amax(t))
                #plt.legend(loc='best')
                #plt.savefig('Plots/Potential='+ str(Potential[a]/10000)+ '/Energie='+str(Energien[e]/10000) +'/Stromerwartungswert(t)_N='+str(int(Anzahl_N[l]))+ 'w = ' + str(Frequenz[f]) + '.pdf')
                #plt.close()
            plt.figure(Figure_Zahler)
#            farbe=np.linspace(0,1,10000)
        I_bar_ges=I_bar_1+I_bar_2
        params , cov = curve_fit(gerade,Frequenz[10:50],np.real(I_bar_ges[10:50]))
        print('a=',params)
        plt.plot(Frequenz,I_bar_ges,'-'+farbe[e],alpha=0.25, label=r'$E_0=$'+str(Energien[e]/10000))
        plt.plot(Frequenz,gerade(Frequenz,*params),':'+farbe[e],label=r'$E_0=$'+str(Energien[e]/10000)+' Fit')
        #       plt.plot(Frequenz,I_bar_lsode,alpha=0.25, label=r'lsode E='+str(Energien[e]/10000))
        #plt.title('Erwartungswert des Stromes für den Startzustand ' + str(np.round(Startzustand,3))+'\n fur n=' + str(Anzahl_N) + '\n w=' + str(Frequenz[f]) + ' E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/10000) )
        #plt.plot(Energien/10000,I_bar_lsode,  alpha=0.25, label=r'Strommittelwert lsode w='+str(Frequenz[f]))
    hilfe=np.linspace(np.amin(I_bar),np.amax(I_bar),2)

        #plt.xlim(0,6)
#    plt.plot(hilfe*0+np.sqrt((Potential[a]/10000)**2+4*1**2)-Potential[a]/10000,hilfe,'k-')
#    plt.plot(hilfe*0+np.sqrt((Potential[a]/10000)**2+4*1**2)+Potential[a]/10000,hilfe,'k-')
#    plt.plot(hilfe*0+2*np.sqrt((Potential[a]/10000)**2+4*1**2),hilfe,'k-')
    plt.legend(loc='best')
    plt.xlim(0,grenze[1])
    plt.ylim(-1*10**(-8),1*10**(-8))
    plt.xlabel(r'$\omega/\frac{J}{\hbar}$ ')
    plt.ylabel(r'$\bar{\langle I \rangle}/ \frac{J\symup{e}}{\hbar}$')
    plt.tight_layout()
    plt.savefig('Plots_mittelwerte/Potential='+ str(Potential[a]/10000)+'Stromerwartungswert(t)_N='+str(int(Anzahl_N))+ '.pdf')
    Figure_Zahler=1+Figure_Zahler

#         #print(Eigenwerte)
#         #print(np.size(Eigenwerte)/2)
#
#     #print(epsilon)
#     #print(Eigenvektoren_realteil,Eigenvektoren_imagteil)
#
#     plt.title('zeitentwicklung für den Startzustand ' + str(phi[:,0]))
#     plt.plot(t,np.real(psi_t[2]),'-b',label=r'3 real')
#     plt.plot(t,np.imag(psi_t[2]),'-b',alpha=0.3,label=r'3 imag')
#     plt.plot(t,np.real(psi_t[3]),'-r',label=r'4 real')
#     plt.plot(t,np.imag(psi_t[3]),'-r',alpha=0.3,label=r'4 imag')
#     plt.plot(t,np.real(psi_t[0]),'-g',label=r'1 real')
#     plt.plot(t,np.imag(psi_t[0]),'-g',alpha=0.3,label=r'1 imag')
#     plt.plot(t,np.real(psi_t[1]),'-y',label=r'2 real')
#     plt.plot(t,np.imag(psi_t[1]),'-y',alpha=0.3,label=r'2 imag')
#     #
#     # #
#     # plt.plot(t,zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t),'-g',label=r'2')
#     #
#     # plt.plot(t,zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t),'-g',label=r'3')
#     # plt.plot(t,zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t)[0],'-g',label=r'4')
#     #
#     plt.legend(loc='best')
#     plt.savefig('Plots/zeitentwicklung phi_1 N='+str(int(Anzahl_N[l]))+ 'a='+str(Potential/100000)+'.pdf')
#
#
#
#
#
# #    print(betragsquadrad(psi_1,psi_1,phi,epsilon,t))
#     print(Startzustand)
#     t=np.linspace(0,10000,100000)
#     psi_t=zeitentwicklung(Startzustand,V_phi,epsilon,Frequenz,t)
#
#     plt.figure(Figure_Zahler )
#     Figure_Zahler=1+Figure_Zahler
#
#     #betragsquadrad(messwerte,Startzustand,V_phi,epsilon,frequenz,t):
#     plt.title('zeitentwicklung für den Startzustand ' + str(np.round(Startzustand,3))+ '\n w=' + str(Frequenz) + 'E='+str(Energien/10000) + 'a='+str(Potential/10000)  )
#     plt.plot(t,betragsquadrad(phi[:,0],psi_t),'-r',alpha=0.5,label=r'zustand 1')
#     plt.plot(t,betragsquadrad(phi[:,1],psi_t),'-y',alpha=0.5,label=r'zustand 2')
#     plt.plot(t,betragsquadrad(phi[:,2],psi_t),'-g',alpha=0.5,label=r'zustand 3')
#     plt.plot(t,betragsquadrad(phi[:,3],psi_t),'-b',alpha=0.5,label=r'zustand 4')
#     for n in range(4):
#         plt.plot(T * 0 + Periodendauer * n, T, '--k',linewidth=0.5)
#     plt.legend(loc='best')
#     plt.savefig('Plots/zeitentwicklung_der_Besetzungen_N='+str(int(Anzahl_N[l]))+ 'a='+str(Potential/100000)+'.pdf')
#
#
#     plt.figure(Figure_Zahler )
#     Figure_Zahler=1+Figure_Zahler
#     plt.title('zeitentwicklung für den Startzustand ' + str(np.round(Startzustand,3))+ '\n w=' + str(Frequenz) + 'E='+str(Energien/10000) + 'a='+str(Potential/10000)  )
#     plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 0], psi_t),'-r',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[0],3)))
#     plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 1], psi_t),'-y',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[1],3)))
#     plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 2], psi_t),'-g',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[2],3)))
#     plt.plot(t,betragsquadrad(H_0_eigenvektoren[:, 3], psi_t),'-b',alpha=0.5,label=r'eig_vek zu Eigenwert'+str(np.round(H_0_Eigenwerte[3],3)))
#     for n in range(4):
#         plt.plot(T * 0 + Periodendauer * n, T, '--k')
#     plt.legend(loc='best')
#     plt.savefig('Plots/zeitentwicklung_der_Eigenzustande_N='+str(int(Anzahl_N[l]))+ 'a='+str(Potential/100000)+'.pdf')

##################################################

#
#   print(betragsquadrad(psi_1,psi_1,phi,epsilon,t)+betragsquadrad(psi_2,psi_1,phi,epsilon,t)+betragsquadrad(psi_3,psi_1,phi,epsilon,t)+betragsquadrad(psi_4,psi_1,phi,epsilon,t))

#
#print(np.dot(phi[:,0],phi[:,1]))

#   if [(np.size(Frequenz_100000)==1)and(np.size(Potential)==1)]:
#   r i in range(np.size(phi[0,:])):
#     test(phi,i+1)


#print(Eigenvektoren_realteil,Eigenvektoren_imagteil)
#print(Eigenvektoren)






# def eig_erwartet_Matrix(Anzahl,Frequenz,Potential):
#     Eigenwerte_H_0=np.genfromtxt('Parameter/eigenwerte_von_H_0_fur_a='+str(int(Potential)) +'a.txt')
#     m=np.size(Eigenwerte_H_0)*(1+Anzahl*2)
#     Matrix=np.zeros((m,np.size(Frequenz)))
#     Frequenz=Frequenz*Potential/10000
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
# #        Matrix=np.zeros((np.size(Grosse),np.size(Frequenz_100000)))
#     return Matrix

# Eigenwerte_H_0=np.array([-2,0,0,2])
# Figure_Zahler=1
#
# #Figure_Zahler=Figure_Zahler+1
# plt.title('Eigenwerte von a='+str(Potential/10000)+' E='+str(Frequenz) )
# plt.figure(Figure_Zahler)
# for index_N , value_N in enumerate(Anzahl_N):
#         Eigenwerte = np.genfromtxt('build/Eigenwerte_fur_a='+str(int(Potential))+'_w='+str(int(Frequenz_100000))+'_E='+str(int(Energien))+'_N='+str(int(value_N))+'.txt')
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
# plt.savefig('Plots/Plot_fur'+'_a='+str(Potential/10000)+'_w='+str(Frequenz) +'_E='+str(int(Energien))+'.pdf')
# plt.close()
#
        # for j in range(n):
        #     r=j%2
        #     b=(j+1)%2
        #
        #     plt.plot(Energien/10000,Matrix_mit_Eigenwerten[j,:],'-b',color=(r,0,b),alpha=0.5)#,label='Eigenwert'+ str(j) )
            #plt.plot(Frequenz,Matrix_mit_Eigenwerten[j,:],'xr',alpha=0.5)#,label='Eigenwert'+ str(j) )



        #     Erwartete_Eigenwerte=eig_war(Eigenwerte_H_0,Anzahl,Frequenz,value_Potential)
        #     n=np.size(Eigenwerte)
        #     for i in range(n-1):
        #         plt.plot(Frequenz,Eigenwerte[i],'xb')
        #         plt.plot(Frequenz,Erwartete_Eigenwerte[i],'+r')
        #     plt.plot(Frequenz,Eigenwerte[n-1],'xb',label=r'E='+str(value_Energie/10000))
        #     plt.plot(Frequenz,Erwartete_Eigenwerte[n-1],'+r',label=r'Erwartete_Eigenwerte')
        #     print('Potential',value_Potential/10000)
        #     print('value_Energie',value_Energie/10000)
        #     print('Frequenz=', Frequenz)
        #     print('Erwartete_Eigenwerte',np.sort(Erwartete_Eigenwerte))
        #     print('Berechnete Eigenwerte',Eigenwerte)
