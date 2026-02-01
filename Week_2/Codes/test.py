import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


### Choosing the initial Parameters
# Flight Params
rho     = 1.225
g       = 9.81
b       = 2
AR      = 9.686999647
S       = b**2 / AR
payload_mass_drop = 0.5
payload_mass = 2.1

V_climb = 20
V_cruise = 22
V_des = 20

cruise_dist_1 = 30000
cruise_dist_2 = 30000

height_cruise = 200

gamma_climb = 8
gamma_des = 8

C_d_0 = 0.03
e = 0.8

energy_density = 237600

Initial_weight = 8
g = 9.8
def Power_climb(V_cl,gamma_climb,Weight):
    C_l = Weight*np.cos(np.deg2rad(gamma_climb))/(0.5*rho*V_cl**2*S)
    C_d = C_d_0 + C_l**2/(np.pi*AR*e)
    D_cl = 0.5*rho*V_cl**2*S*C_d
    T_cl = D_cl + Weight*np.sin(np.deg2rad(gamma_climb))
    P_cl = T_cl*V_cl
    return P_cl

def Power_cruise(V_cr,Weight):
    P_cr = Power_climb(V_cr,0,Weight)
    return P_cr

def Power_descent(V_des,gamma_des,Weight):
    P_des = max(Power_climb(V_des,gamma_des,Weight),0)
    return P_des

def Power_loiter(V_loi,Weight):
    P_loi = Power_cruise(V_loi,Weight)
    return P_loi
t_cr_1 = cruise_dist_1/V_cruise
t_cr_2 = cruise_dist_2/V_cruise

t_cl_1 = height_cruise/(V_climb*np.sin(np.deg2rad(gamma_climb)))
t_des_1 = height_cruise/(V_des*np.sin(np.deg2rad(gamma_des)))

t_cr_1,t_cr_2,t_cl_1,t_des_1
def Total_energy(Weight):
    energy_climb = Power_climb(V_climb,gamma_climb,Weight)*t_cl_1
    energy_cruise_1 = Power_cruise(V_cruise,Weight)*t_cr_1
    energy_cruise_2 = Power_cruise(V_cruise,Weight-payload_mass_drop*g)*t_cr_2
    energy_descent = Power_descent(V_des,gamma_des,Weight-payload_mass_drop*g)*t_des_1
    Total_energy = energy_climb + energy_cruise_1 + energy_cruise_2 + energy_descent
    return Total_energy
Total_energy(10*g)/energy_density
### Empty Weight Graph
Data = pd.read_excel("Drone_dataset.xlsx")
Data
Mtow_array = Data["MTOW (kg)"].to_numpy()
Empty_weight_array = Data["Empty Weight (kg)"].to_numpy()

Empty_weight_ratios = Empty_weight_array/Mtow_array
logx = np.log10(Mtow_array)
logy = np.log10(Empty_weight_ratios)
A = np.column_stack((logx,np.ones_like(logx)))
B = logy

Soln, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

L = Soln[0]
loga = Soln[1]
a = 10**loga

print("L =", L)
print("loga =", loga)
print("a =", a)
a*12**L
W_0 = Initial_weight
Weight_array = []

for i in range(100):
    Weight_array.append(W_0)
    W_b = Total_energy(W_0*g)/energy_density
    W_0 = (payload_mass/(1 - (W_b/W_0) - a*W_0**L)).item()

print(Weight_array)

plt.figure(1)
plt.plot(Weight_array)
plt.show()