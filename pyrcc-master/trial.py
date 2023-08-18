import pyrcc
import numpy as np
from sklearn import datasets
from sklearn.metrics import adjusted_rand_score

iris = datasets.load_iris()
X = iris.data[:, :2]
Y = iris.target

# pendata = np.genfromtxt('pendigits/pendigits.txt', delimiter=',')

# X = pendata[:,:-1]
# Y = pendata[:,-1]

clusterer = pyrcc.RccCluster(measure='cosine')

P = clusterer.fit(X)

print('ARI: {}'.format(adjusted_rand_score(Y, P)))