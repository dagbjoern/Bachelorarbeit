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

e1 , e2 ,e3, e4,e5 =np.genfromtxt('eigenwerte.txt',unpack=True)

x=np.linspace(0,len(e1)-1,len(e1))
print('cool,cool,cool')
plt.figure(1)
plt.plot(x,e1,'-b',label=r'1')
plt.plot(x,e3,'-c',label=r'2')
plt.plot(x,e2,'-r',label=r'3')
plt.plot(x,e4,'-g',label=r'4')
plt.plot(x,e5,'-m',label=r'5')

plt.legend(loc='best')
plt.savefig('build/eigenwerte.pdf')
