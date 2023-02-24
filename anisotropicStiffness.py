# Importing appropriate modules
import math
import numpy as np
import matplotlib.pyplot as plt

# Material propertires - All values are given
theta = math.radians(19.9)

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

# Compliance Matrix, S
sMatrix = [
			[(1/E_1), (-nu_21/E_2), 0],
			[(-nu_12/E_1), (1/E_2), 0],
			[(0), (0), (1/G_12)]
		  ]

# C-Matrix is simply just the inverse of our S comliance matrix
cMatrix = np.linalg.inv(sMatrix)

# Stress Matrix
sigma_x = [
			[sigma_xx], 
			[sigma_yy], 
			[tau_xy]
		  ]

# Variant 1
sigma_1 = np.matmul(T1, sigma_x)

epsilon_1_var1 = np.matmul(sMatrix, sigma_1)
inverseT2 = np.linalg.inv(T2)
epsilon_x_var1 = np.matmul(inverseT2, epsilon_1_var1)

gamma_12 = (sigma_1[2][0])/(G_12)

# Variant 2
s_x = sMatrix[0]

epsilon_x_var2 = np.matmul(sMatrix, sigma_x)
epsilon_1_var2 = np.matmul(sMatrix, sigma_1)

# Variant Comparison
print(sigma_1)
print(sigma_x)
print("epsilon_1_var1: ")
print(epsilon_1_var1)
print("epsilon_x_var1: ")
print(epsilon_x_var1)
print("epsilon_x_var2: ")
print(epsilon_x_var2)
print("epsilon_x_var2: ")
print(epsilon_1_var2)

# Defining epsilon and sigma arrays

angles = []
epsilon_11 = []
epsilon_22 = []
gamma_12 = []
sigma_11 = []
sigma_22 = []
tau_12 =[]
epsilon_xx = []
epsilon_yy = []
gamma_xy = []

# Variant 1 for-loop

for x in range (91):
	angles.append(x)
	thetaPrime = math.radians(x)
	theta = x

	m = math.cos(thetaPrime)
	n = math.sin(thetaPrime)

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

	sigma_1 = np.matmul(T1, sigma_x)
	epsilon_1_var1 = np.matmul(sMatrix, sigma_1)
	inverseT2 = np.linalg.inv(T2)
	epsilon_x_var1 = np.matmul(inverseT2, epsilon_1_var1)

	epsilon_11.append(epsilon_1_var1[0]*100)
	epsilon_22.append(epsilon_1_var1[1]*100)
	gamma_12.append(epsilon_1_var1[2]*100)

	sigma_11.append(sigma_1[0])
	sigma_22.append(sigma_1[1])
	tau_12.append(sigma_1[2])

	epsilon_xx.append(epsilon_x_var1[0]*100)
	epsilon_yy.append(epsilon_x_var1[1]*100)
	gamma_xy.append(epsilon_x_var1[2]*100)

# Variant 1 - Epsilon_1 Graph
plt.figure()
plt.plot(angles, epsilon_11, label='Epsilon_11')
plt.plot(angles, epsilon_22, label='Epsilon_22')
plt.plot(angles, gamma_12, label='Gamma_12')

plt.xlabel('Theta (deg)')
plt.ylabel('% Strain')
plt.legend()

plt.title('Variant 1 - % Strain vs Angle | Non-Material Orientation')

# Variant 1 - Stress_1 Graph
plt.figure()
plt.plot(angles, sigma_11, label='Sigma_11')
plt.plot(angles, sigma_22, label='Sigma_22')
plt.plot(angles, tau_12, label='Tau_12')

plt.xlabel('Theta (deg)')
plt.ylabel('Stress (Pa)')
plt.legend()

plt.title('Variant 1 - Stress (Pa) vs Angle | Material Orientation')

# Variant 1 - Epsilon_x Graph
plt.figure()
plt.plot(angles, epsilon_xx, label='Epsilon_xx')
plt.plot(angles, epsilon_yy, label='Epsilon_yy')
plt.plot(angles, gamma_xy, label='Gamma_xy')

plt.xlabel('Theta (deg)')
plt.ylabel('% Strain')
plt.legend()

plt.title(' Variant 1 - % Strain vs Angle | Non-Material Orientation')

# Variant 2
epsilon_11_var2 = []
epsilon_22_var2 = []
gamma_12_var2 = []
sigma_11_var2 = []
sigma_22_var2 = []
tau_12_var2 =[]
epsilon_xx_var2 = []
epsilon_yy_var2 = []
gamma_xy_var2 = []

# Variant 2 for-loop

for x in range (91):
	thetaPrime = math.radians(x)
	theta = x

	m = math.cos(thetaPrime)
	n = math.sin(thetaPrime)

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

	sigma_1_var2 = np.matmul(T1, sigma_x)
	s_x = sMatrix[0]
	epsilon_x_var2 = np.matmul(sMatrix, sigma_x)
	epsilon_1_var2 = np.matmul(sMatrix, sigma_1_var2)

	epsilon_11_var2.append(epsilon_1_var2[0]*100)
	epsilon_22_var2.append(epsilon_1_var2[1]*100)
	gamma_12_var2.append(epsilon_1_var2[2]*100)

	sigma_11_var2.append(sigma_1_var2[0])
	sigma_22_var2.append(sigma_1_var2[1])
	tau_12_var2.append(sigma_1_var2[2])

	epsilon_xx_var2.append(epsilon_x_var2[0]*100)
	epsilon_yy_var2.append(epsilon_x_var2[1]*100)
	gamma_xy_var2.append(epsilon_x_var2[2]*100)

# Variant 2 - Epsilon_1 Graph
plt.figure()
plt.plot(angles, epsilon_11_var2, label='Epsilon_11')
plt.plot(angles, epsilon_22_var2, label='Epsilon_22')
plt.plot(angles, gamma_12_var2, label='Gamma_12')

plt.xlabel('Theta (deg)')
plt.ylabel('% Strain')
plt.legend()

plt.title('Variant 2 - % Strain vs Angle | Non-Material Orientation')

# Variant 2 - Stress_1 Graph
plt.figure()
plt.plot(angles, sigma_11_var2, label='Sigma_11')
plt.plot(angles, sigma_22_var2, label='Sigma_22')
plt.plot(angles, tau_12_var2, label='Tau_12')

plt.xlabel('Theta (deg)')
plt.ylabel('Stress (Pa)')
plt.legend()

plt.title('Variant 2 - Stress (Pa) vs Angle | Material Orientation')

# Variant 3 - Epsilon_x Graph
plt.figure()
plt.plot(angles, epsilon_xx_var2, label='Epsilon_xx')
plt.plot(angles, epsilon_yy_var2, label='Epsilon_yy')
plt.plot(angles, gamma_xy_var2, label='Gamma_xy')

plt.xlabel('Theta (deg)')
plt.ylabel('% Strain')
plt.legend()

plt.title('Variant 2 - % Strain vs Angle | Non-Material Orientation')
plt.show()