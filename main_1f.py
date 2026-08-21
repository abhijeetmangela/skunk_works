#
import numpy as np
import matplotlib.pyplot as plt

#
'''
Keywords:

t: total
_c: centroidal coordinates
s: support
r: root
2: x direction (horizontal)
3: y direction (vertical)
yc: camber
u: upper
l: lower
L: total length
d: differential
_sum: integral (riemann sum)

'''

#
'''
Mistakes
Copy Pasting Code and Modifying It (Might lead to Errors if anything is Missed)
Minimize the Amount of Function Imports & No of Cells
Label Every Plot as it is Made
Double Counting Check
Distinct Variable Names

Concept
Cummulative Summation

'''

#
'''
New Functions

np.concatenate
np.vstack
np.diff
np.cumsum
np.squeeze
plt.colorbar
plt.quiver

'''

#
#1 Shear Force and Bending Moment Diagram

Rs = -800
Ms = -920

xi = np.linspace(0, 0.4, 4)
xj = np.linspace(0.4, 1, 6)

V3i = 800
V3j = 1000

M2i = -920 + 800*xi
M2j = -1000 + 1000*xj


# Plot

# plt.plot(xi, V3i*np.ones_like(xi), label='x = (0, 0.4)')
# plt.plot(xj, V3j*np.ones_like(xj), label='x = (0.4, 1)')
# plt.title('Shear Force Diagram')
# plt.ylim((0, 1200))
# plt.ylabel('Shear Force: V (N)')
# plt.xlabel('Position: x (m)')
# plt.legend()
# plt.grid()
# plt.show()

# plt.plot(xi, M2i, label='x = (0, 0.4)')
# plt.plot(xj, M2j, label='x = (0.4, 1)')
# plt.title('Bending Moment Diagram')
# plt.legend()
# plt.ylabel('Bending Moment: M (N)')
# plt.xlabel('Position: x (m)')
# plt.ylim((-1000, 100))
# plt.grid()
# plt.show()

#
# import math as mt

# # NACA 4 Digit Airfoil Generation

# m = 2   # maximum camber
# p = 3   # maximum camber position
# t = 19  # thickness

# c = 0.15   # chord length

# m = 0.01*m  # maximum camber in % of chord
# p = 0.10*p  # maximum camber position in tenths of chord
# t = 0.01*t  # thickness in % of chord


# # Coefficients for 4 digit series
# a0 =  1.4845
# a1 = -0.6300
# a2 = -1.7580
# a3 =  1.4215
# a4 = -0.5075

# n = 199 # number of points along the chord
# x = np.linspace(0,c,n) # x coordinate of points along the chord
# y   = np.zeros(n) # x coordinate of points along the chord
# yc  = np.zeros(n) # y coordinate of the camber line
# dyc = np.zeros(n) # gradient of the camber line
# yt  = np.zeros(n) # thickness distribution
# xu  = np.zeros(n) # x coordinate of the upper surface
# yu  = np.zeros(n) # y coordinate of the upper surface
# xl  = np.zeros(n) # x coordinate of the lower surface
# yl  = np.zeros(n) # y coordinate of the lower surface
# for i in range(n):
#     if  (x[i]/c < p):
#         yc[i]  = (c*m/p**2)*(2*p*(x[i]/c)-(x[i]/c)**2)
#         dyc[i] = ((2*m)/p**2)*(p-(x[i]/c))
#     else:
#         yc[i]  = (c*m/(1-p)**2)*(1-2*p+2*p*(x[i]/c)-(x[i]/c)**2)
#         dyc[i] = ((2*m)/(1-p)**2)*(p-(x[i]/c))
        
# for i in range(n):
#     yt[i] = (t*c)*(a0*mt.sqrt(x[i]/c)+a1*(x[i]/c)+a2*(x[i]/c)**2+a3*(x[i]/c)**3+a4*(x[i]/c)**4)
#     teta  = mt.atan(dyc[i])
#     xu[i] = x[i]  - yt[i]*mt.sin(teta)
#     xl[i] = x[i]  + yt[i]*mt.sin(teta)
#     yu[i] = yc[i] + yt[i]*mt.cos(teta)
#     yl[i] = yc[i] - yt[i]*mt.cos(teta)


# plt.xlim(-0.2,c+0.2)
# plt.ylim(-c/3,c/3)
# plt.scatter(xu,yu,color='deepskyblue', s=1)   
# plt.scatter(xl,yl,color='deepskyblue', s=1)
# plt.scatter(x,yc, s=1) 
# plt.axis('equal')
# # plt.yticks([])
# # plt.xticks([])
# plt.show()

#
# NACA XFLR

naca2319 = np.loadtxt("s7055.dat")

xul = naca2319[:, 0].T*0.15
yul = naca2319[:, 1].T*0.15

n = xul.shape[0]

# print(xul.shape)
# plt.scatter(xul, yul, s=2)
# plt.axis('equal')
# plt.show()

nh = int(n/2)

xu = xul[1:nh]
xu = xu[::-1]
xl = xul[nh+1:]
yu = yul[1:nh]
yl = yul[nh+1:]
yu = yu[::-1]

xle, yle = xul[nh], yul[nh]
xte, yte = xul[0], yul[0]


#
#2 NACA 2319

# Length Summation

dlu = np.sqrt(np.diff(xu)**2 + np.diff(yu)**2)
dll = np.sqrt(np.diff(xl)**2 + np.diff(yl)**2)

Lu = np.sum(dlu)
Ll = np.sum(dll)

Lle = np.sqrt((xle - xu[0])**2 + (yle - yu[0])**2) + np.sqrt((xle - xl[0])**2 + (yle - yl[0])**2)
Lte = np.sqrt((xte - xu[-1])**2 + (yte - yu[-1])**2) + np.sqrt((xte - xl[-1])**2 + (yte - yl[-1])**2)

L = Lu + Ll + Lte + Lle
# print(f'Total Length (Perimeter) of Airfoil: {L:.4g} m')

ta = 0.3e-3
As = L*ta
# print(f'Sectional Area of Airfoil: {As:.4g} m^2')

# 
#2 Finding Centroid

xmu = (xu[1:] + xu[:-1]) / 2
ymu = (yu[1:] + yu[:-1]) / 2

xml = (xl[1:] + xl[:-1]) / 2
yml = (yl[1:] + yl[:-1]) / 2

dAu = ta * dlu
dAl = ta * dll

xdA_le = (xu[0] + xle)/2*np.sqrt((xu[0] - xle)**2 + (yu[0] - yle)**2)*ta + (xl[0] + xle)/2*np.sqrt((xl[0] - xle)**2 + (yl[0] - yle)**2)*ta
xdA_te = (xte + xu[-1])/2*np.sqrt((xte - xu[-1])**2 + (yu[-1] - yte)**2)*ta + (xte + xl[-1])/2*np.sqrt((xte - xl[-1])**2 + (yl[-1] - yte)**2)*ta

ydA_le = (yu[0] + yle)/2*np.sqrt((xu[0] - xle)**2 + (yu[0] - yle)**2)*ta + (yl[0] + yle)/2*np.sqrt((xl[0] - xle)**2 + (yl[0] - yle)**2)*ta
ydA_te = (yu[-1] + yte)/2*np.sqrt((xte - xu[-1])**2 + (yu[-1] - yte)**2)*ta + (yte + yl[-1])/2*np.sqrt((xte - xl[-1])**2 + (yl[-1] - yte)**2)*ta


xdA_sum = np.sum(xmu*dAu) + np.sum(xml*dAl) + xdA_le + xdA_te
ydA_sum = np.sum(ymu*dAu) + np.sum(yml*dAl) + ydA_le + ydA_te

# print(xdA_sum)
# print(ydA_sum)

x_c, y_c = xdA_sum/As, ydA_sum/As
# print(f"Centroid Location {x_c:.5f}, {y_c:.5f} m")

#
#2 Labeled Plot of Airfoil Section

# plt.plot(xu, yu, 'b')   
# plt.plot(xl, yl, 'b')
# plt.plot(xle, yle, 'b')
# plt.plot(xte, yte, 'b')
# # plt.plot(x, yc, 'r--', label='Camber')
# plt.scatter(x_c, y_c, c='g', s=7, label=f'Centroid: {x_c:.4f}, {y_c:.4f}')
# plt.title('NACA 2319')
# plt.grid()
# plt.legend()
# plt.axis('equal')
# plt.show()

#
#3 Second Area Moments

# Centroidal Coordinates

xmu_c = xmu - x_c
ymu_c = ymu - y_c

xml_c = xml - x_c
yml_c = yml - y_c

x2dA_le = ((xu[0] + xle)/2 - x_c)**2*np.sqrt((xu[0] - xle)**2 + (yu[0] - yle)**2)*ta + ((xl[0] + xle)/2 - x_c)**2*np.sqrt((xl[0] - xle)**2 + (yl[0] - yle)**2)*ta
x2dA_te = ((xte + xu[-1])/2 - x_c)**2*np.sqrt((xte - xu[-1])**2 + (yu[-1] - yte)**2)*ta + ((xte + xl[-1])/2 - x_c)**2*np.sqrt((xte - xl[-1])**2 + (yl[-1] - yte)**2)*ta

y2dA_le = ((yu[0] + yle)/2 - y_c)**2*np.sqrt((xu[0] - xle)**2 + (yu[0] - yle)**2)*ta + ((yl[0] + yle)/2 - y_c)**2*np.sqrt((xl[0] - xle)**2 + (yl[0] - yle)**2)*ta
y2dA_te = ((yu[-1] + yte)/2 - y_c)**2*np.sqrt((xte - xu[-1])**2 + (yu[-1] - yte)**2)*ta + ((yte + yl[-1])/2 - y_c)**2*np.sqrt((xte - xl[-1])**2 + (yl[-1] - yte)**2)*ta

xydA_le = ((xu[0] + xle)/2 - x_c)*((yu[0] + yle)/2 - y_c)*np.sqrt((xu[0] - xle)**2 + (yu[0] - yle)**2)*ta + ((xl[0] + xle)/2 - x_c)*((yl[0] + yle)/2 - y_c)*np.sqrt((xl[0] - xle)**2 + (yl[0] - yle)**2)*ta
xydA_te = ((yu[-1] + yte)/2 - y_c)*((xu[-1] + xte)/2 - x_c)*np.sqrt((xte - xu[-1])**2 + (yu[-1] - yte)**2)*ta + ((yte + yl[-1])/2 - y_c)*((xte + xl[-1])/2 - x_c)*np.sqrt((xte - xl[-1])**2 + (yl[-1] - yte)**2)*ta


I22 = np.sum(ymu_c**2 * dAu) + np.sum(yml_c**2 * dAl) + y2dA_le + y2dA_te
I33 = np.sum(xmu_c**2 * dAu) + np.sum(xml_c**2 * dAl) + x2dA_le + x2dA_te
I23 = np.sum(xmu_c*ymu_c*dAu) + np.sum(xml_c*yml_c*dAl) + xydA_le + xydA_te

# print(f"I22 = {I22:.5g} m^4")
# print(f"I33 = {I33:.5g} m^4")
# print(f"I23 = {I23:.5g} m^4")

#
xmu_le_c = np.array([(xle+xu[0])/2 - x_c, (xle+xl[0])/2 - x_c])
xmu_te_c = np.array([(xte+xu[-1])/2 - x_c, (xte+xl[-1])/2 - x_c])

xmt_c = np.concatenate((np.array([xmu_te_c[0]]), xmu_c[::-1], xmu_le_c, xml_c, np.array([xmu_te_c[-1]])), axis=0)
# print(xmt_c.shape)

ymu_le_c = np.array([(yle+yu[0])/2 - y_c, (yle+yl[0])/2 - y_c])
ymu_te_c = np.array([(yte+yu[-1])/2 - y_c, (yte+yl[-1])/2 - y_c])

ymt_c = np.concatenate((np.array([ymu_te_c[0]]), ymu_c[::-1], ymu_le_c, yml_c, np.array([ymu_te_c[-1]])), axis=0)

# [markdown]
# ![alt text](image.png)

#
#4 Axial Stress Distribution (Root)

Ey = 70e9
nu = 0.36
M2r = -920
M3r = 0

Ieqv = 1/(I22*I33 - I23**2)
Imat = np.array([[I33, I23],
                 [I23, I22]])

def axialstress_calc(M2_, M3_):
    Mvec = np.array([[M2_],
                    [M3_]])

    X_stack = np.vstack((ymt_c, -xmt_c))
    # print(X_stack.shape)
    X_c = X_stack.T
    # print(X_c.shape)

    S11 = (X_c @ Imat @ Mvec) * Ieqv
    # print(S11.shape)
    return S11

#
#4 Visualize Axial Stress Distribution

# plt.scatter(0, 0, c='r', marker='x', label='Centroid')
# plt.scatter(xmt_c, ymt_c, c=np.squeeze(S11)*1e-9, cmap='jet', s=10)
# plt.colorbar(label=r'$\sigma_{11}$ Axial Stress (GPa)')
# plt.axis('equal')
# plt.title(r'Axial Stress $\sigma_{11}$ Distribution on Airfoil Section')
# plt.xlabel('X (Axis-2) Coordinate (m)')
# plt.ylabel('Y (Axis-3) Coordinate (m)')
# plt.legend()
# plt.grid()
# plt.show()

# [markdown]
# ![image.png](attachment:image.png)

#
#5 Shear Flow Distribution (Tip)

V2t = 0
V3t = 1000

Imat2 = np.array([[I22, -I23],
                  [-I23, I33]])

def shearflow_calc(V2_, V3_):
    Vvec = np.array([[V2_],
                    [V3_]])

    dxu = np.diff(xu[::-1])
    dxl = np.diff(xl)

    dyu = np.diff(yu[::-1])
    dyl = np.diff(yl)

    dx_te = np.array([xu[-1] - xte, xte - xl[-1]])
    dx_le = np.array([xle - xu[0], xl[0] - xle])

    dy_te = np.array([yu[-1] - yte, yte - yl[-1]])
    dy_le = np.array([yle - yu[0], yl[0] - yle])

    dxt = np.concatenate((np.array([dx_te[0]]), dxu, dx_le, dxl, np.array([dx_te[-1]])), axis=0)
    dyt = np.concatenate((np.array([dy_te[0]]), dyu, dy_le, dyl, np.array([dy_te[-1]])), axis=0)

    dst = np.sqrt(dxt**2 + dyt**2)

    txtds = xmt_c*ta*dst
    tytds = ymt_c*ta*dst


    tXds = np.vstack((txtds, tytds)).T
    dqb = -(tXds @ Imat2 @ Vvec) * Ieqv

    qb = np.cumsum(dqb, axis=0)
    # print(qb.shape)

    # [markdown]
    # ![image.png](attachment:image.png)

    #
    #5 Determining q0 & Plotting qs

    q0 = -np.sum(np.squeeze(qb)*dst) / np.sum(dst)
    # print(q0)

    qs = np.squeeze(qb) + q0
    # print(qs.shape)
    return qs, dxt, dyt

# Plot
# plt.scatter(0, 0, c='r', marker='x', label='Centroid')
# plt.scatter(xmt_c, ymt_c, c=qs, cmap='jet', s=10)
# plt.colorbar(label=r'$q(s)$ Shear Flow (N/m)')
# plt.axis('equal')
# plt.title(r'Shear Flow $q(s)$ Distribution on Airfoil Section')
# plt.xlabel('X (Axis-2) Coordinate (m)')
# plt.ylabel('Y (Axis-3) Coordinate (m)')
# plt.legend()
# plt.grid()
# plt.show()