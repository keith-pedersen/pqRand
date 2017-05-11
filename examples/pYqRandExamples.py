import sys
sys.path.append("../")
from pYqRand import *
import matplotlib.pyplot as plt

gen = engine()
dist = log_normal(0., 1.)

hist = plt.hist([dist(gen) for _ in xrange(10**6)], 100, normed=1, histtype="step")
plt.show()
