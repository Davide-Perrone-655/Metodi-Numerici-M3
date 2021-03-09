import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from Veryfunlib import veryfunlib as vfl

#avvio il generatore random, usare MT19937 al posto di PCG64 per il mersenne twister

vfl.setseed(15)

#mi dice quanta roba c'è
camp, delta, eta = np.loadtxt('shape_10.dat', unpack=True, max_rows=1)
camp=int(camp)


#carico i dati
path_data=np.zeros(camp)

path_data = np.loadtxt('misure_prova.dat', unpack=True)



#faccio le medie e le varianze

path_mean=np.mean(path_data)
path_var=np.sqrt(np.var(path_data, ddof=1))


print('medie')
print(path_mean)
print('-')
print(path_var)


bootlen=10
print(camp)
print(2**bootlen)

#bootstrap
#path_boot=vfl.bootstrap(path_data, bootlen, est='mean', rep=30, R_print=True)
path_block=vfl.blocking(path_data, bootlen, R_print=True)

print(path_block)
	







