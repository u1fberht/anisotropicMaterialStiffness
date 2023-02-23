# Importing appropriate modules

import math
import numpy as np
import matplotlib.pyplot as plt

# Material propertires

theta = 19.9 * (2*(math.pi)/360)

E1 = 180 * 1000000000
E2 = 9 * 1000000000
G12 = 4 * 1000000000
nu12 = 0.25
nu21 = nu12*(E2/E1)

m = math.cos(theta)
n = math.sin(theta)

# Transpose Matrices

T1 = [[pow(m,2), pow(n,2), 2*m*n],
	  [pow(n,2), pow(m,2), -2*m*n],
	  [-m*n, m*n, pow(m,2)-pow(n,2)]]

T2 = [[pow(m,2), pow(n,2), m*n],
	  [pow(n,2), pow(m,2), -m*n],
	  [-2*m*n, 2*m*n, pow(m,2)-pow(n,2)]]

print(theta)