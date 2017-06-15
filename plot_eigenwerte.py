import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
# from uncertainties.unumpy import (nominal_values as noms,
#                                   std_devs as stds)
# from scipy import stats
# from scipy.optimize import curve_fit
# import scipy.constants as const
# from uncertainties import ufloat
# from scipy.integrate import trapz
# from scipy.integrate import simps
# from tabulate import tabulate


Energien=np.genfromtxt('build/Durchlaufende_Energien.txt')
Potential=np.genfromtxt('build/Durchlaufende_Potentiale.txt')
frequenzzahler , frequenznenner =np.genfromtxt('build/Durchlaufende_Frequenzen.txt')

erwartete_Eigenwerte=[-2,0,0,2,3,5,5,7,-3,-5,-5,-7]

Figure_Zahler=1

for index_Energie , value_Energie in enumerate(Energien):
    for index_Potential , value_Potenial in enumerate(Potential):
        Eigenwerte=np.genfromtxt('build/Eigenwerte_fur_a='+str(int(value_Potenial))+'E='+str(int(value_Energie))+'.txt',unpack=True)
        Gitterkonstante,Frequenz,Anzahl,Phasenverschiebung=np.genfromtxt('build/Parameter_fur_a='+str(int(value_Potenial))+'E='+str(int(value_Energie))+'.txt',unpack=True)

        x=np.linspace(1,len(Eigenwerte),len(Eigenwerte))
        plt.figure(Figure_Zahler)
        Figure_Zahler=Figure_Zahler+1
        plt.plot(x,Eigenwerte,'xb',label=r'E='+str(value_Energie/100))
        for i , value in enumerate(erwartete_Eigenwerte):
            plt.plot(x,value+x*0,'--k')
        plt.plot(x,Eigenwerte,'xr',label=r'E=Test')
        plt.xlabel('Eigenwerte')
        plt.ylabel('E/J')
        plt.savefig('Plots/Plotfur_a='+str(value_Potenial/100)+'_E='+str(value_Energie/100) +'.pdf')
        plt.close()

#     Energien =np.genfromtxt('build/Energien.txt',unpack=True)
#
#
#
# E_all=np.array([e1,e2,e3,e4,e5])
# print(E_all[1,:])
#
# for i in range(0,len(E_all)):
#     #name=np.array([i])
#     #name=name.astype(str)
#     #test=np.array(i)
#     e_test=np.array(E_all[i,:])
#     for index, value in enumerate(e_test):
#         e_test[index]= e_test[0]+index*Frequenz
#          plt.plot(x,value+x*0,'--k')
# #         plt.plot(x,e3,'-c',label=r'E=1')
# #         #plt.plot(x,e2,'-r',label=r'E=2')
# #         #plt.plot(x,e4,'-g',label=r'E=3')
# #         # plt.plot(x,e5,'-m',label=r'E=4')
#     plt.legend(loc='best')
#     plt.savefig('build/Eigenwerte_Energie'+str(Energien[i]) +'.pdf')
# cool=np.array([5])
# a=cool.astype(str)
# print(a)
# print('cool'+ a[0])
#
# plt.figure(2)
# plt.plot(e1,e1*0,'xb',label=r'E=0')
# plt.plot(e3,e1*0,'xc',label=r'E=1')
# plt.plot(e2,e1*0,'xr',label=r'E=2')
# plt.plot(e4,e1*0,'xg',label=r'E=3')
# plt.plot(e5,e1*0,'xm',label=r'E=4')
# plt.legend(loc='best')
# plt.savefig('build/Eigenwerte_Energien1.pdf')
#
#
# e1 , e2 ,e3, e4, e5 =np.genfromtxt('Eigenwerte_Frequenz.txt',unpack=True)
#
# x=np.linspace(0,len(e1)-1,len(e1))
#
# plt.figure(3)
# plt.plot(x,e1,'-b',label=r'w=1')
# plt.plot(x,e3,'-c',label=r'w=2')
# plt.plot(x,e2,'-r',label=r'w=3')
# plt.plot(x,e4,'-g',label=r'w=4')
# plt.plot(x,e5,'-m',label=r'w=5')
# plt.legend(loc='best')
# plt.savefig('build/Eigenwerte_Frequenz.pdf')
#
#
# plt.figure(4)
# plt.plot(e1,e1*0,'xb',label=r'E=0')
# plt.plot(e3,e1*0,'xc',label=r'E=1')
# plt.plot(e2,e1*0,'xr',label=r'E=2')
# plt.plot(e4,e1*0,'xg',label=r'E=3')
# plt.plot(e5,e1*0,'xm',label=r'E=4')
# plt.legend(loc='best')
# plt.savefig('build/Eigenwerte_Frequenz1.pdf')
#
#
# e1 , e2 ,e3, e4, e5 =np.genfromtxt('Eigenwerte_Gitter.txt',unpack=True)
#
# x=np.linspace(0,len(e1)-1,len(e1))
#
# plt.figure(5)
# plt.plot(x,e1,'-b',label=r'a=1')
# plt.plot(x,e3,'-c',label=r'a=2')
# plt.plot(x,e2,'-r',label=r'a=3')
# plt.plot(x,e4,'-g',label=r'a=4')
# plt.plot(x,e5,'-m',label=r'a=5')
# plt.legend(loc='best')
# plt.savefig('build/Eigenwerte_Gitter.pdf')
#
# plt.figure(6)
# plt.plot(e1,e1*0,'xb',label=r'E=0')
# plt.plot(e3,e1*0,'xc',label=r'E=1')
# plt.plot(e2,e1*0,'xr',label=r'E=2')
# plt.plot(e4,e1*0,'xg',label=r'E=3')
# plt.plot(e5,e1*0,'xm',label=r'E=4')
# plt.legend(loc='best')
# plt.savefig('build/Eigenwerte_Gittter1.pdf')
