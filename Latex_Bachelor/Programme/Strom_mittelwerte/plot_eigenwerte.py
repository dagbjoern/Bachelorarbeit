import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from tqdm import tqdm
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

def quadrat(x,a):
    return a*x**2

Energien=np.genfromtxt('build/Durchlaufende_Energien.txt')
Potential=np.genfromtxt('build/Durchlaufende_Potentiale.txt')
Frequenz_1000=np.genfromtxt('build/Durchlaufende_Frequenzen.txt',unpack=True)
Anzahl_N=np.genfromtxt('build/Durchlaufende_N.txt',unpack=True)


str_Potential = Potential.astype(int)
str_Potential = str_Potential.astype(str)

str_Frequenz_1000 = Frequenz_1000.astype(int)
str_Frequenz_1000 = str_Frequenz_1000.astype(str)

str_Energien = Energien.astype(int)
str_Energien = str_Energien.astype(str)

str_Anzahl_N = Anzahl_N.astype(int)
str_Anzahl_N = str_Anzahl_N.astype(str)



if not os.path.exists('Plots_mittelwerte/'):
    os.makedirs('Plots_mittelwerte/')

Figure_Zahler = 1
a=np.linspace(0,np.amax(Potential)/100,100)
Frequenz = Frequenz_1000/1000

plt.figure(Figure_Zahler)
Figure_Zahler=Figure_Zahler+1
#plt.title('zeitentwicklung für den Startzustand ' + str(np.round(Startzustand,3))+'fur n=' + str(Anzahl_N[l]) + '\n w=' + str(Frequenz[f]) + ' E=' +str(Energien[e]/10000) + ' a=' +str(Potential[a]/100) )
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




#Gitterkonstante, Anzahl, Phasenverschiebung = np.genfromtxt('Parameter/Parameter_fur_a='+str_Potential[0]+'_w='+ str_Frequenz_1000[0]+'_E='+str_Energien[0]+'_N='+str_Anzahl_N[0]+'.txt',unpack=True)

print('int',Anzahl_N )
print('str', str_Anzahl_N)
print('cool' + str_Anzahl_N[0] + 'cool')
t_lsode=np.genfromtxt('build/Zeit_lsode.txt')

I_bar=Energien*0 #array der länge n E erschaffen

farbe=np.array(['r','b','g','m'])


for a in tqdm(range(np.size(Potential))):
    H_0_eigenvektoren=np.genfromtxt('Parameter/eigenvektoren_von_H_0_fur_a='+str_Potential[a]+'.txt')
    H_0_Eigenwerte=np.genfromtxt('Parameter/eigenwerte_von_H_0_fur_a='+str_Potential[a] +'.txt')
    if not os.path.exists('Plots/Potential='+ str(Potential[a]/100)):
        os.makedirs('Plots/Potential='+ str(Potential[a]/100))
    for f in tqdm(range(int(np.size(Frequenz)))):
        I_bar=Energien*0+1j*0 #array der länge von E erschaffen
        I_bar_lsode=Energien*0 #array der länge von E erschaffen
        for e in tqdm(range(np.size(Energien))):
            if not os.path.exists('Plots/Potential='+ str(Potential[a]/100) + '/Energie='+str(Energien[e]/10000)):
                os.makedirs('Plots/Potential='+ str(Potential[a]/100)+ '/Energie='+str(Energien[e]/10000))
            #               Numerche Werte Einlese
            psi_real_lsode=np.genfromtxt('build/Realpart_Eigenvektoren_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_1000[f] +'lsode.txt')
            psi_imag_lsode=np.genfromtxt('build/Imagpart_Eigenvektoren_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_1000[f] +'lsode.txt')
            psi_imag_lsode=psi_imag_lsode*1j
            psi_t_lsode=psi_real_lsode+psi_imag_lsode
            t=t_lsode
            psi_t_lsode=psi_t_lsode.T
            for l in tqdm(range(np.size(Anzahl_N))):
                #Eigenwerte=np.genfromtxt('build/Eigenwerte_fur_a='+str_Potential[a]+'_w='+ str_Frequenz_1000[f]+'_E='+str_Energien[e]+'_N='+str_Anzahl_N[l]+'.txt',unpack=True)
                epsilon=np.genfromtxt('build/epsilon_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_1000[f] +'_N=' + str_Anzahl_N[l] +'.txt')
                Eigenzustande_realteil=np.genfromtxt('build/Realpart_Eigenzustande_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_1000[f] +'_N=' + str_Anzahl_N[l] +'.txt')
                Eigenzustande_imagteil=np.genfromtxt('build/Imagpart_Eigenzustande_fur_a='+str_Potential[a]+'_E='+str_Energien[e]+'_w=' + str_Frequenz_1000[f] +'_N=' + str_Anzahl_N[l] +'.txt')
                Eigenzustande_realteil=Eigenzustande_realteil+1j*0
                Eigenzustande_imagteil=Eigenzustande_imagteil*1j
                Startzustand=H_0_eigenvektoren[:,0]
                V_phi=Eigenzustande_realteil+Eigenzustande_imagteil
                #print('Strommittelwert für Potential='+ str(Potential[a]/100)+'\nEnergie='+str(Energien[e]/10000)+'\nw = ' + str(Frequenz[f])  ,'Strommittel=',Stromittelwert_Foquet(Startzustand, 1, V_phi))
                I_bar[e]=Stromittelwert_Foquet(Startzustand, 1, V_phi)
                y_lsode=Erwartungswert(Strom(1),psi_t_lsode,t_lsode)
                #print('numerischer Wert:', np.sum(y_lsode)/np.size(t_lsode))
                I_bar_lsode[e]=np.sum(y_lsode)/np.size(t_lsode)
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
        params , cov = curve_fit(quadrat,Energien/10000,np.real(I_bar))
        plt.plot(Energien/10000,I_bar, '-'+farbe[f], alpha=0.25, label=r'$ \omega_F=$'+str(Frequenz[f]))
#        plt.plot(Energien/10000,quadrat(Energien/10000,*params),label=r'fit_w='+str(Frequenz[f]) )
        plt.plot(Energien/10000,I_bar_lsode, ':'+farbe[f], alpha=1, label=r'$\omega_M=$'+str(Frequenz[f]))
        plt.xlabel(r'$E_0/ \frac{J}{d\symup{e}}$')
        plt.ylabel(r'$\bar{\langle I \rangle} / \frac{J\symup{e}}{\hbar} $')
        plt.legend(loc='best')
    plt.xlim(0,0.1)
    plt.tight_layout()
    plt.savefig('Plots_mittelwerte/Potential='+ str(Potential[a]/100)+'Stromerwartungswert(t)_N='+str(int(Anzahl_N[l]))+ '.pdf')
    Figure_Zahler=1+Figure_Zahler
