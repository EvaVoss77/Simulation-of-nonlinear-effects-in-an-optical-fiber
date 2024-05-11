import functions
# Pulse parameters
T0 = 20         # Pulse width [ps]
ts = 0.2        # Sampling time
n = 15
N = np.power(2, n)
T = N * ts       # Window width [ps]
t = np.arange(-N/2, N/2) * ts  # Temporal vector [ps]
A = 10                                # Power pulse

# Transmission parameters
beta_2 = -21 # ps^2/km
alph = 0.04
L_D = T0 ** 2 / abs(beta_2)
z = np.arange(0, 40000)
dz = 10
z_tot = 50

# Pulse on temporal domain
U0 = np.sqrt(A) * np.exp(-t * t / (2 * T0 * T0))
# Frequency domain
FT_U0 = np.fft.fft(U0)
frec = np.arange(-N/2, N/2) / (ts * N)

# Plot the pulse in the time and frequency domain at different propagation points for a case in which there is not attenuation and a nonlinear component $\beta_2$.
fig, ax = plt.subplots(int(z_tot/dz) + 1, 2)
fig.set_figheight(12)
fig.set_figwidth(16)

# Frequency domine
j = 0
ax[j, 0].set_title("Effect of $\\beta_2$ on pulse propagation in frequency domain for $\\alpha = 0$")
for i in range(0, z_tot + 1, dz):
  U = np.fft.fftshift(F_disp2(FT_U0, z[i], beta_2, ts, N))
  ax[j, 0].plot(frec, np.abs(U)**2, label='z =' + str(i) + 'km')
  ax[j, 0].set_xlim([-0.05, 0.05])
  ax[j, 0].set_ylabel("Amplitude")
  ax[j, 0].legend()
  j += 1
ax[j - 1, 0].set_xlabel("Frequency [ps]")

# Time domine
j = 0
ax[j, 1].set_title("Effect of $\\beta_2$ on pulse propagation in the time domain for $\\alpha = 0$")
for i in range(0, z_tot + 1, dz):
  U = np.fft.ifft(F_disp2(FT_U0, z[i], beta_2, ts, N))
  ax[j, 1].plot(t, np.abs(U)**2, label='z =' + str(i) + ' km')
  ax[j, 1].legend()
  ax[j, 1].set_xlim([-300, 300])
  ax[j, 1].set_ylabel("Amplitude")
  j += 1
ax[j - 1, 1].set_xlabel("Time [ps]")
plt.show()

# Plot the pulse in the time and frequency domain at different propagation points for a case in which there is attenuation and a nonlinear component $\beta_2$.
fig, ax2 = plt.subplots(int(z_tot/dz) + 1, 2)
fig.set_figheight(12)
fig.set_figwidth(16)

# Frequency domine
j = 0
ax2[j, 0].set_title("Effect of $\\beta_2$ on pulse propagation in frequency domain for $\\alpha \\neq 0$")
for i in range(0, z_tot + 1, dz):
  U = np.fft.fftshift(F_disp2_at(FT_U0,i,beta_2,alph,ts,N))
  ax2[j, 0].plot(frec, np.abs(U)**2, label='z =' + str(i) + 'km')
  ax2[j, 0].set_xlim([-0.05, 0.05])
  ax2[j, 0].set_ylim([-1E5,7E5])
  ax2[j, 0].set_ylabel("Amplitude")
  ax2[j, 0].legend()
  j += 1
ax2[j - 1, 0].set_xlabel("Frequency [ps]")

# Time domine
j = 0
ax2[j, 1].set_title("Effect of $\\beta_2$ on pulse propagation in the time domain for $\\alpha \\neq 0$")
for i in range(0, z_tot + 1, dz):
  U = np.fft.ifft(F_disp2_at(FT_U0, i, beta_2, alph, ts, N))
  ax2[j, 1].plot(t, np.abs(U)**2, label='z =' + str(i) + ' km')
  ax2[j, 1].legend()
  ax2[j, 1].set_xlim([-300, 300])
  ax2[j, 1].set_ylim([-1,11])
  ax2[j, 1].set_ylabel("Amplitude")
  j += 1

ax2[j - 1, 1].set_xlabel("Time [ps]")
plt.show()



#################### Verify that the energy of the pulse decreases exponentially ###############
alph = 0.04
z = np.arange(0, 80)
Energy = []
for c in z:
  U = np.fft.ifft(F_disp2_at(FT_U0, c, beta_2, alph, ts, N))
  e = np.trapz(np.abs(U) ** 2) * 1E-3
  Energy.append(e)

fig, ax = plt.subplots(1)
fig.set_figheight(5)
fig.set_figwidth(8)

plt.plot(z, Energy)
plt.xlabel('Distance Traveled [km]')
plt.ylabel('Energy [pJ]')
plt.title('Exponential Decay of Pulse Energy over Distance')



################# Verify that the pulse widens as indicated by the analytical function ###########
L_D = T0**2 / abs(beta_2)  # km
z = np.arange(0, 100)
Width_sim = []
Width_teo = T0 * np.sqrt(1 + (z / L_D) ** 2)

print(GaussianWidth(U0 * U0, t))
i = 0
for c in z:
    U = np.fft.ifft(F_disp2_at(FT_U0, c, beta_2, alph, ts, N))
    Width_sim.append(GaussianWidth(np.abs(U) ** 2, t))

fig, ax = plt.subplots()
fig.set_figheight(5)
fig.set_figwidth(8)
plt.title('Width of the pulse with attenuation')
ax.plot(z, Width_sim, label='Simulation')
ax.plot(z, Width_teo, label='Analytical')
plt.legend()
plt.xlabel('Distance traveled [km]')
plt.ylabel('Pulse Width [ps]')