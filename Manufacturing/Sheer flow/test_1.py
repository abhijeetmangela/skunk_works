import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def find_shear(data2, a_boom2, Vy, L):

    # --------------------------------------------------------
    # Calculate centroid of boom idealisation
    # --------------------------------------------------------

    cg_new = np.zeros(2)

    cg_new[0] = (
        np.sum(
            data2[0, :L - 1] * a_boom2
        )
        / np.sum(a_boom2)
    )

    cg_new[1] = (
        np.sum(
            data2[1, :L - 1] * a_boom2
        )
        / np.sum(a_boom2)
    )


    # --------------------------------------------------------
    # Shift coordinates to boom centroid
    # --------------------------------------------------------

    data3 = np.zeros_like(data2)

    data3[0, :] = data2[0, :] - cg_new[0]
    data3[1, :] = data2[1, :] - cg_new[1]


    # --------------------------------------------------------
    # Calculate section properties
    # --------------------------------------------------------

    Ixy2 = np.sum(
        data3[0, :L - 1] *
        data3[1, :L - 1] *
        a_boom2
    )

    Ix2 = np.sum(
        data3[1, :L - 1]**2 *
        a_boom2
    )

    Iy2 = np.sum(
        data3[0, :L - 1]**2 *
        a_boom2
    )


    den2 = (
        Ix2 * Iy2 -
        Ixy2**2
    )


    Kxy2 = Ixy2 / den2
    Kx2 = Ix2 / den2
    Ky2 = Iy2 / den2


    # --------------------------------------------------------
    # Calculate first moments
    # --------------------------------------------------------

    Qy2 = np.zeros(L - 1)
    Qx2 = np.zeros(L - 1)

    qs2 = np.zeros(L - 1)


    # MATLAB:
    # for i = 1:L-2

    for i in range(L - 2):

        Qy2[i + 1] = (
            Qy2[i] +
            data3[0, i + 1] *
            a_boom2[i + 1]
        )

        Qx2[i + 1] = (
            Qx2[i] +
            data3[1, i + 1] *
            a_boom2[i + 1]
        )

        qs2[i + 1] = Vy * (
            Kxy2 * Qy2[i + 1]
            -
            Ky2 * Qx2[i + 1]
        )


    return qs2

# ============================================================
# Skin without spar
# ============================================================

# Clear equivalent
plt.close('all')


# ------------------------------------------------------------
# Read airfoil data
# ------------------------------------------------------------

T = pd.read_csv(
    'fx63137.dat.txt',
    sep=r'\s+',
    header=None
)

data1 = T.to_numpy(dtype=float)

c = 0.234
data = data1 * c

L = len(data)


# ------------------------------------------------------------
# Calculate panel lengths and midpoint coordinates
# ------------------------------------------------------------

dis = np.zeros(L - 1)
mid_data = np.zeros((L - 1, 2))

for i in range(L - 1):
    dis[i] = abs(
        np.sqrt(
            (data[i + 1, 0] - data[i, 0])**2 +
            (data[i + 1, 1] - data[i, 1])**2
        )
    )

    mid_data[i, :] = 0.5 * (data[i, :] + data[i + 1, :])


# ------------------------------------------------------------
# Calculate centroid
# ------------------------------------------------------------

cg = np.zeros(2)

cg[0] = np.sum(mid_data[:, 0] * dis) / np.sum(dis)
cg[1] = np.sum(mid_data[:, 1] * dis) / np.sum(dis)


# ------------------------------------------------------------
# Shift coordinates to centroid
# ------------------------------------------------------------

data[:, 0] = data[:, 0] - cg[0]
data[:, 1] = data[:, 1] - cg[1]


# ------------------------------------------------------------
# Calculate neutral axis angle
# ------------------------------------------------------------

I_num = 0.0
I_den = 0.0

for i in range(L - 1):

    I_num += (
        (data[i, 0] + data[i + 1, 0]) *
        (data[i, 1] + data[i + 1, 1]) *
        0.25 * dis[i]
    )

    I_den += (
        (data[i, 1] + data[i + 1, 1]) *
        (data[i, 1] + data[i + 1, 1]) *
        0.25 * dis[i]
    )

tan_alpha = I_num / I_den


# ------------------------------------------------------------
# Distance from neutral axis
# ------------------------------------------------------------

p_d = np.zeros(L - 1)

for i in range(L - 1):
    p_d[i] = abs(
        data[i, 1] - tan_alpha * data[i, 0]
    ) / np.sqrt(1 + tan_alpha**2)


# ------------------------------------------------------------
# Calculate boom areas
# ------------------------------------------------------------

a_boom = np.zeros(L - 1)

for i in range(L - 1):

    if i < L - 2:

        a_boom[i] += (
            (1 / 6) * dis[i] * 0.0005 *
            (2 + p_d[i + 1] / p_d[i])
        )

        a_boom[i + 1] += (
            (1 / 6) * dis[i] * 0.0005 *
            (2 + p_d[i] / p_d[i + 1])
        )

    else:

        a_boom[i] += (
            (1 / 6) * dis[i] * 0.0005 *
            (2 + p_d[0] / p_d[i])
        )

        a_boom[0] += (
            (1 / 6) * dis[i] * 0.0005 *
            (2 + p_d[i] / p_d[0])
        )


# ============================================================
# Applied shear force
# ============================================================

Vy = float(
    input('The shear force to be considered [in N] is: ')
)

data2 = data.T

qs = find_shear(data2, a_boom, Vy, L)


# ------------------------------------------------------------
# Plot shear flow without spar
# ------------------------------------------------------------

plt.figure(figsize=(14, 8))

half = int(0.5 * (L - 1))

plt.plot(
    data2[0, :half] + cg[0],
    -qs[:half],
    linewidth=2,
    color='blue',
    label='Upper Surface'
)

plt.plot(
    data2[0, half:L - 1] + cg[0],
    -qs[half:L - 1],
    linewidth=2,
    color='red',
    label='Lower Surface'
)

plt.grid(True)

plt.xlabel('x [m]')
plt.ylabel(r'$q_s$ [N/m]')

plt.legend(
    fontsize=18,
    prop={'family': 'Helvetica', 'weight': 'bold'}
)

plt.title(
    'Shear Flow Distribution without spar',
    fontsize=24,
    fontweight='bold'
)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()


# ============================================================
# Zero net shear flow correction
# ============================================================

q0 = -np.sum(qs * dis) / np.sum(dis)

qs = qs + q0

shear = qs


# ============================================================
# Stringer requirement
# ============================================================

a = float(
    input('The ribs spacing in m: ')
)

D = (
    np.pi**2 *
    70 *
    (10**9) *
    (0.0005**3) *
    1.62
) / (
    12 *
    2.5 *
    1.5 *
    1.5 *
    1.5 *
    (1 - 0.33**2)
)

Kss = (
    5.34 +
    (4 * 0.234 * 0.234) / (a**2)
)

N = (
    Kss * D
) / (
    0.234 * 0.234
)


print('Maximum shear flow =', np.max(np.abs(shear)))
print('Critical shear flow N =', N)

if N > np.max(np.abs(shear)):
    print('No stringer needed without spar')
else:
    print('Stringer needed without spar')


# ============================================================
# Add stringers
# ============================================================

str_area = 10 * 0.001 * 0.001

a_boom2 = a_boom.copy()


iu = int(
    input('# stringers you want on upper surface: ')
)

iu_ind = list(map(
    int,
    input(
        'iu stringer indices for 1 to 48 points: '
    ).split()
))


il = int(
    input('# stringers you want on lower surface: ')
)

il_ind = list(map(
    int,
    input(
        'il stringer indices for 1 to 48 points: '
    ).split()
))


# ------------------------------------------------------------
# Upper surface stringers
# ------------------------------------------------------------

for i in range(iu):

    # MATLAB:
    # a_boom2(iu_ind(i+1)) = a_boom2(iu_ind(i+1)) + str_area;

    # Python requires -1 because MATLAB is 1-based
    index = iu_ind[i + 1] - 1

    a_boom2[index] += str_area


# ------------------------------------------------------------
# Lower surface stringers
# ------------------------------------------------------------

for i in range(il):

    # MATLAB:
    # il_ind(i+1) + 0.5*(L+1)

    index = (
        il_ind[i + 1]
        + int(0.5 * (L + 1))
        - 1
    )

    a_boom2[index] += str_area


# ============================================================
# Shear flow with stringers
# ============================================================

qs2 = find_shear(
    data2,
    a_boom2,
    Vy,
    L
)

q02 = -np.sum(qs2 * dis) / np.sum(dis)


# ------------------------------------------------------------
# Determine maximum shear flow in each panel
# ------------------------------------------------------------

q_max2 = np.zeros(iu + il + 2)
dis2 = np.zeros(iu + il + 2)


print('\nUpper stringer location')

for i in range(iu + 1):

    start = iu_ind[i] - 1
    end = iu_ind[i + 1] - 1

    q_max2[i] = np.max(
        np.abs(qs2[start:end] + q02)
    )

    dis2[i] = np.sum(dis[start:end])

    print(
        data[start, 0] + cg[0]
    )


print('\nLower stringer location')

# Shift lower indices as in MATLAB
il_ind_shifted = [
    x + int(0.5 * (L - 1))
    for x in il_ind
]

for i in range(il + 1):

    start = il_ind_shifted[i]
    end = il_ind_shifted[i + 1]

    q_max2[iu + 1 + i] = np.max(
        np.abs(qs2[start:end] + q02)
    )

    dis2[iu + 1 + i] = np.sum(
        dis[start:end]
    )

    print(
        data[start, 0] + cg[0]
    )


# ------------------------------------------------------------
# Stringer critical shear flow
# ------------------------------------------------------------

Kss2 = (
    5.34 +
    (4 * dis2**2) / (a**2)
)

N2 = (
    Kss2 * D
) / (
    dis2**2
)

print('\nN2 / q_max2 =')
print(N2 / q_max2)


# ============================================================
# Plot shear flow with stringers
# ============================================================

plt.figure(figsize=(14, 8))

plt.plot(
    data2[0, :half] + cg[0],
    -qs2[:half],
    linewidth=2,
    color='blue',
    label='Upper Surface'
)

plt.plot(
    data2[0, half:L - 1] + cg[0],
    -qs2[half:L - 1],
    linewidth=2,
    color='red',
    label='Lower Surface'
)

plt.grid(True)

plt.xlabel('x [m]')
plt.ylabel(r'$q_s$ [N/m]')

plt.legend(
    fontsize=18,
    prop={'family': 'Helvetica', 'weight': 'bold'}
)

plt.title(
    'Shear Flow Distribution without spar with stringers',
    fontsize=24,
    fontweight='bold'
)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()


# ============================================================
# Plot airfoil with panels
# ============================================================

plt.figure(figsize=(14, 8))

# Upper panels
for i in range(iu + 1):

    start = iu_ind[i] - 1
    end = iu_ind[i + 1]

    plt.plot(
        data[start:end, 0] + cg[0],
        data[start:end, 1] + cg[1],
        linewidth=2
    )

    mid = round(
        0.5 * (start + end - 1)
    )

    plt.text(
        data[mid, 0] + cg[0],
        data[mid, 1] + cg[1] + 0.00234,
        f'Panel {i + 1}',
        fontsize=18,
        fontweight='bold'
    )


# Lower panels
for i in range(il + 1):

    start = il_ind_shifted[i]
    end = il_ind_shifted[i + 1] + 1

    plt.plot(
        data[start:end, 0] + cg[0],
        data[start:end, 1] + cg[1],
        linewidth=2
    )

    mid = round(
        0.5 * (start + end - 1)
    )

    plt.text(
        data[mid, 0] + cg[0],
        data[mid, 1] + cg[1] + 0.00234,
        f'Panel {iu + 1 + i}',
        fontsize=18,
        fontweight='bold'
    )


plt.axis('equal')
plt.grid(True)

plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.title(
    'Section 3 with Stringers',
    fontsize=24,
    fontweight='bold'
)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()


# ============================================================
# Plot idealised airfoil
# ============================================================

plt.figure(figsize=(14, 8))

plt.plot(
    data[:, 0] + cg[0],
    data[:, 1] + cg[1],
    linewidth=2,
    color='blue'
)

plt.scatter(
    data[:, 0] + cg[0],
    data[:, 1] + cg[1],
    s=50,
    marker='o',
    facecolors='red'
)

plt.axis('equal')
plt.grid(True)

plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.title(
    'FX63137 Airfoil Idealised',
    fontsize=24,
    fontweight='bold'
)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlim([
    np.min(data[:, 0] + cg[0]) - 0.01,
    np.max(data[:, 0] + cg[0]) + 0.01
])

plt.ylim([
    np.min(data[:, 1] + cg[1]) - 0.011,
    np.max(data[:, 1] + cg[1]) + 0.011
])

plt.savefig(
    'airfoil + booms.png',
    dpi=300,
    bbox_inches='tight'
)

plt.show()


# ============================================================
# FUNCTION: FIND SHEAR
# ============================================================


