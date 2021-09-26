###
# Elspeth Lee - Sep 2021
# Python limb darkening testing code for gCMCRT
# Samples a random transit chord on the transit annulus - projects it back onto the stellar disk,
# then calculates directly the stellar zenith angle of the chord,
# then calculates the "exact" limb darkening coefficent for that chord given a limb darkening law.
# Then makes some neato plots for the number of samples
###

import numpy as np
import matplotlib.pylab as plt
import cartopy.crs as ccrs

# Constants
GM = 1.2668653e17
G =  6.67430e-11
au = 1.495978707e11
Rj = 7.1492e7
Mj = GM/G
Rsun = 6.95700e8
Msun = 1.98847e30


# Number of random samples wanted
nsamp = 250

# limb darkening law integer picker
ilimb = 2

# Coefficents for the limb darkening laws
c1 = 0.077
c2 = 0.311
c3 = 0.0
c4 = 0.0

# System properties - zp = planetary TOA altitude, Rz = total height of atmosphere + planet
Rp = 1.138 * Rj
zp = 0.1 * Rp
Rz = Rp + zp
Rs = 0.805 * Rsun
a = 0.03142 * au

# Inclination in degrees
inc = 85.51

# Phase of planet (0-1) - 0.0 is mid transit
phase = 0.005

print('---------------------------')
print('number of samples: ',nsamp)
print('Limb darkening law integer: ', ilimb)
print('Limb darkening coefficents:',c1,c2,c3,c4)
print('R planet (Rj): ', Rp/Rj)
print('TOA altitude (m) : ', zp)
print('R + TOA alt. (Rj): ', Rz/Rj)
print('R star (Rsun): ', Rs/Rsun)
print('semi-major axis (au): ', a/au)
print('inclination (deg): ', inc)
print('orbital phase: ', phase)
print('---------------------------')

#### Begin calculations ####

# Angle array for limb darkening law
nang = 1800
thetas = np.linspace(0,90,nang)
mus = np.cos(thetas * np.pi/180.0)

# Could use switches here but i'm lazy
if (ilimb == 1):
  # Schwartzchild
  Imus = 1.0 - c1*(1.0 - mus)
elif (ilimb == 2):
  # Quadratic
  Imus = 1.0 - c1*(1.0 - mus) - c2*(1.0 - mus)**2
elif (ilimb == 3):
  # Square root law
  Imus = 1.0 - c1*(1.0 - mus) - c2*(1.0 - np.sqrt(mus))
elif (ilimb == 4):
  # Logarithmic
  Imus = 1.0 - c1*(1.0 - mus) - c2*mus*np.log(mus)
elif (ilimb == 5):
  # Exponental law
  Imus = 1.0 - c1*(1.0 - mus) - c2/(1.0 - np.exp(mus))
elif (ilimb == 6):
  # Sing (2009) three paramater
  Imus = 1.0 - c1*(1.0 - mus) - c2*(1.0 - mus**(3.0/2.0)) - c3*(1.0 - mus**2)
elif (ilimb == 7):
  # Claret (2000) four parameter
  Imus = 1.0 - c1*(1.0 - np.sqrt(mus)) - c2*(1.0 - mus) - c3*(1.0 - mus**(3.0/2.0)) - c4*(1.0 - mus**2)

# Calculate the phase angle in degrees from the 0 longitude, where 90 would be direct facing.
# This makes it analogous to the impact parameter, i call this the phasepact parameter.
# We define negative impact parameter (negative inclination) as going southward, and positive as going northward.
# We define negative phasepact parameter (>0.5 phase) as going westward, and positive as going eastward.
if (phase <= 0.5):
  phase = 90.0 - phase*360.0
else:
  phase = -90.0 - (360.0 - phase*360.0)

inc = inc * np.pi/180.0
phase = phase * np.pi/180.0

# Impact parameter and phasepact paramater
b = -(a * np.cos(inc))/Rs
if (inc < 0.0):
  b = abs(b)

ph = (a * np.cos(phase))/Rs

print('impact parameter: ', -b)
print('phasepact parameter: ', ph)

if (abs(b) > 1.0 or abs(ph) > 1.0):
  print('!!! Warning !!! - Planetary disk outside of stellar disk')

# Find central lon,lat
zs_cent = b * Rs
ys_cent = ph * Rs
xs_cent = np.sqrt(Rs**2 - zs_cent**2 - ys_cent**2)
theta_cent = np.arccos(zs_cent/Rs) - np.pi/2.0
phi_cent = np.arctan2(ys_cent, xs_cent)
mu_cent = np.cos(theta_cent) * np.cos(phi_cent)
Imu_cent = 1.0 - c1*(1.0 - mu_cent) - c2*(1.0 - mu_cent)**2
theta_cent = theta_cent * 180.0/np.pi
phi_cent = phi_cent * 180.0/np.pi

print('Central (phi, theta): ', phi_cent, theta_cent)


# Start sampling random annulli
mu_z = np.zeros(nsamp)
Imu_z = np.zeros(nsamp)
theta_arr = np.zeros(nsamp)
phi_arr = np.zeros(nsamp)

for n in range(nsamp):

    # Sample random position on the transmission annulus
    rr2 = np.sqrt(Rp**2 + (Rz**2 - Rp**2)*np.random.random_sample())

    # Sample uniform disk
    #rr2 = np.sqrt(Rz**2*np.random.random_sample())

    #print(Rp,rr2,Rz)

    # Sample ann theta
    ann_theta = np.random.random_sample() * np.pi * 2.0

    # z,y,x position on planetary sphere
    zp = rr2 * np.cos(ann_theta)
    yp = rr2 * np.sin(ann_theta)
    xp = np.sqrt(Rz**2 - zp**2 - yp**2)

    #print(zp,yp,xp)

    # z,y,z posotion on stellar sphere (assumed parallel beams)
    zs = zp + b * Rs
    ys = yp + ph * Rs
    xs = Rs**2 - zs**2 - ys**2
    if (xs > 0.0):
      xs = np.sqrt(xs)
    else:
      print('Packet outside of stellar disk, calculate anyway cause why not')
      xs = np.sqrt(xs)

    #print(zs,ys,xs)

    # Find longitude and latitude coordinate on star (radians)
    theta = np.arccos(zs/Rs) - np.pi/2.0
    phi = np.arctan2(ys, xs)

    # Convert to degrees
    theta_arr[n] = theta * 180.0/np.pi
    phi_arr[n] = phi * 180.0/np.pi

    #print(n, theta_arr[n], phi_arr[n])

    # Calculate zenith angle of star and the limb darkening law
    mu_z[n] = np.cos(theta) * np.cos(phi)
    # Could use switches here but i'm lazy
    if (ilimb == 1):
      # Schwartzchild
      Imu_z[n] = 1.0 - c1*(1.0 - mu_z[n])
    elif (ilimb == 2):
      # Quadratic
      Imu_z[n] = 1.0 - c1*(1.0 - mu_z[n]) - c2*(1.0 - mu_z[n])**2
    elif (ilimb == 3):
      # Square root law
      Imu_z[n] = 1.0 - c1*(1.0 - mu_z[n]) - c2*(1.0 - np.sqrt(mu_z[n]))
    elif (ilimb == 4):
      # Logarithmic
      Imu_z[n] = 1.0 - c1*(1.0 - mu_z[n]) - c2*mu_z[n]*np.log(mu_z[n])
    elif (ilimb == 5):
      # Exponental law
      Imu_z[n] = 1.0 - c1*(1.0 - mu_z[n]) - c2/(1.0 - np.exp(mu_z[n]))
    elif (ilimb == 6):
      # Sing (2009) three paramater
      Imu_z[n] = 1.0 - c1*(1.0 - mu_z[n]) - c2*(1.0 - mu_z[n]**(3.0/2.0)) - c3*(1.0 - mu_z[n]**2)
    elif (ilimb == 7):
      # Claret (2000) four parameter
      Imu_z[n] = 1.0 - c1*(1.0 - np.sqrt(mu_z[n])) - c2*(1.0 - mu_z[n]) - c3*(1.0 - mu_z[n]**(3.0/2.0)) - c4*(1.0 - mu_z[n]**2)

    #print(mu_z[n], Imu_z[n])

# Chose to plot a specific limb darkening law
fig = plt.figure()

plt.plot(mus,Imus)
plt.scatter(mu_z,Imu_z,s=5)
plt.plot(mu_cent,Imu_cent,'x',markersize=5)
plt.xlabel('mu')
plt.ylabel('I(mu)')
plt.title('LD Coefficent')


fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global()
ax.plot(phi_arr,theta_arr,'o',transform=ccrs.PlateCarree(),markersize=1)
ax.plot(phi_cent,theta_cent,'x',transform=ccrs.PlateCarree(),markersize=5)
gl = ax.gridlines(draw_labels=True)
plt.title('Flat projection')

# fig = plt.figure(figsize=(10, 5))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.EqualEarth())
# ax.set_global()
# ax.plot(phi_arr,theta_arr,'o',transform=ccrs.PlateCarree(),markersize=1)
# ax.plot(phi_cent,theta_cent,'x',transform=ccrs.PlateCarree(),markersize=5)
# gl = ax.gridlines(draw_labels=True)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(0,0))
ax.set_global()
ax.plot(phi_arr,theta_arr,'o',transform=ccrs.PlateCarree(),markersize=1)
ax.plot(phi_cent,theta_cent,'x',transform=ccrs.PlateCarree(),markersize=5)
gl = ax.gridlines(draw_labels=True)
plt.title('Orthographic projection (0,0)')

#plt.savefig('Ortho_0_0.png', dpi=300)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(phi_cent,theta_cent))
ax.set_global()
ax.plot(phi_arr,theta_arr,'o',transform=ccrs.PlateCarree(),markersize=1)
ax.plot(phi_cent,theta_cent,'x',transform=ccrs.PlateCarree(),markersize=5)
gl = ax.gridlines(draw_labels=True)
plt.title('Orthographic projection - (phi_cent,theta_cent)')

plt.show()
