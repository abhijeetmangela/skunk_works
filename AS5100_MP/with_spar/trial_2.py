# ---------------------------------------------------------------------
# SPAR IDEALISATION
# ---------------------------------------------------------------------

spar_height = 0.028
flange_length = 0.03
flange_thickness = 1e-3
web_thickness = 1e-3

# Spar indices on the existing skin discretisation
spar_upper_index = 25
spar_lower_index = 57

# Spar upper/lower coordinates
spar_upper_x = x[spar_upper_index]
spar_upper_y = y[spar_upper_index]

spar_lower_x = x[spar_lower_index]
spar_lower_y = y[spar_lower_index]


# ---------------------------------------------------------------------
# 1. ADD SPAR FLANGES TO EXISTING SKIN BOOMS
# ---------------------------------------------------------------------

# Area of one flange
A_flange = flange_length * flange_thickness

# Copy the skin boom array
B_total = B_skin.copy()

# Add upper and lower flange areas to the existing skin booms
B_total[spar_upper_index] += A_flange
B_total[spar_lower_index] += A_flange


# ---------------------------------------------------------------------
# 2. WEB BOOMS
# ---------------------------------------------------------------------

num_web_booms = 10

# Web area
A_web = web_thickness * spar_height

# Area assigned to each web boom
A_web_boom = A_web / num_web_booms

# Coordinates of the 10 web booms
# Exclude the flange locations themselves
web_y = np.linspace(
    spar_upper_y,
    spar_lower_y,
    num_web_booms + 2
)[1:-1]

# Web is vertical, so x is constant
web_x = np.full(num_web_booms, spar_upper_x)

# Areas of the web booms
B_web = np.full(num_web_booms, A_web_boom)


# ---------------------------------------------------------------------
# 3. ADD WEB BOOMS TO THE COMPLETE BOOM ARRAY
# ---------------------------------------------------------------------

x_total = np.concatenate([
    x,
    web_x
])

y_total = np.concatenate([
    y,
    web_y
])

B_total = np.concatenate([
    B_total,
    B_web
])


# ---------------------------------------------------------------------
# CHECK
# ---------------------------------------------------------------------

print(f"Flange area       = {A_flange:.6e} m^2")
print(f"Total web area    = {A_web:.6e} m^2")
print(f"Web boom area     = {A_web_boom:.6e} m^2")
print(f"Number of booms   = {len(B_total)}")