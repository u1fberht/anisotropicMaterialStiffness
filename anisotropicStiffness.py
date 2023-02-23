# Importing appropriate modules

import math
import numpy as np
import matplotlib.pyplot as plt

# Material propertires - All values are given

theta = 19.9 * (2*(math.pi)/360)

m = math.cos(theta)
n = math.sin(theta)

E_1 = 180000 * 1000000
E_2 = 9000 * 1000000
G_12 = 4000 * 1000000
nu_12 = 0.25
nu_21 = nu_12*(E_2/E_1)

sigma_xx = 1534 * 1000000
sigma_yy = 831 * 1000000
tau_xy = 409 * 1000000

# Stress Transformation

T1 = [
		[pow(m,2), pow(n,2), 2*m*n],
		[pow(n,2), pow(m,2), -2*m*n],
		[-m*n, m*n, pow(m,2)-pow(n,2)]
	 ]

T2 = [
		[pow(m,2), pow(n,2), m*n],
		[pow(n,2), pow(m,2), -m*n],
		[-2*m*n, 2*m*n, pow(m,2)-pow(n,2)]
	 ]

# Complaince Matrix, S

sMatrix = [
			[(1/E_1), (-nu_21/E_2), 0],
			[(-nu_12/E_1), (1/E_2), 0],
			[(0), (0), (1/G_12)]
		  ]

# C-Matrix is simply just the inverse of our S comliance matrix

cMatrix = np.linalg.inv(sMatrix)

s_1 = [sMatrix[0]]
s_2 = [sMatrix[1]]
s_3 = [sMatrix[2]]


# Stress Matrix

sigma_x = [
			[sigma_xx], 
			[sigma_yy], 
			[tau_xy]
		  ]

sigma_1 = np.matmul(T1, sigma_x)

# Epsilon Matrix

epsilon_x = np.matmul(sMatrix, sigma_x)
epsilon_1 = np.matmul(T2, epsilon_x)

# Constructing shear strain variable

gamma_12 = (sigma_1[2][0])/(G_12)
print(gamma_12)