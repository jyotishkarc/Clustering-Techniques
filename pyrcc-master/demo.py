import pyrcc
import numpy as np
from sklearn.metrics import adjusted_rand_score

X = []
Y = []

pendata = np.genfromtxt('pendigits/pendigits.txt', delimiter=',')

X = pendata[:,:-1]
Y = pendata[:,-1]

clusterer = pyrcc.RccCluster(measure='cosine')

P = clusterer.fit(X)

print('ARI: {}'.format(adjusted_rand_score(Y, P)))