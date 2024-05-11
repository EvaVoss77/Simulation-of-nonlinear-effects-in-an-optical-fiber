from functions import *

##################################### Observation of pulses with attenuation, dispersion of 2nd and 3rd order #################
beta_2 = -20
beta_3 = 10
z_tot = 100
dz = 20

fig, ax3 = plt.subplots(int(z_tot/dz) + 1, 2)
fig.set_figheight(12)
fig.set_figwidth(16)

# Frequency domain
j = 0
ax3[j, 0].set_title("Effect of $\\beta_2 + \\beta_3$ on pulse propagation in frequency domain for $\\alpha = 0$")
for i in range(0, z_tot + 1, dz):
  U = np.fft.fftshift(F_disp23(FT_U0, i, beta_2, beta_3, ts, N))
  ax3[j, 0].plot(frec, np.abs(U)**2, label='z =' + str(i) + 'km')
  ax3[j, 0].set_xlim([-0.05, 0.05])
  ax3[j, 0].set_ylim([-1E5, 7E5])
  ax3[j, 0].set_ylabel("Amplitude")
  ax3[j, 0].legend()
  j += 1
ax3[j - 1, 0].set_xlabel("Frequency [ps]")

# Time domain
j = 0
ax3[j, 1].set_title("Effect of $\\beta_2 + \\beta_3$ on pulse propagation in the time domain for $\\alpha = 0$")
for i in range(0, z_tot + 1, dz):
  U = np.fft.ifft(F_disp23(FT_U0, i, beta_2, beta_3, ts, N))
  ax3[j, 1].plot(t, np.abs(U)**2, label='z =' + str(i) + ' km')
  ax3[j, 1].legend()
  ax3[j, 1].set_xlim([-300, 300])
  ax3[j, 1].set_ylim([-1, 11])
  ax3[j, 1].set_ylabel("Amplitude")
  j += 1

ax3[j - 1, 1].set_xlabel("Time [ps]")
plt.show()

########################  Variation of pulse width under the effect of $\beta_2$ + $\beta_3$ ###############################################
z = np.arange(0,200)
Width_sim = []
Width_teo = T0*np.sqrt(1+(z/L_D)**2)


for c in z:
  U = np.fft.ifft(F_disp23(FT_U0,c,beta_2,beta_3,ts,N))
  Width_sim.append(GaussianWidth(np.abs(U)**2,t))

fig, ax = plt.subplots(1)
fig.set_figheight(5)
fig.set_figwidth(8)
plt.title('Width of the pulse over the effect of $\\beta_2 + \\beta_3$')
plt.xlabel('Distance traveled [km]')
plt.ylabel('Pulse Width [ps]')
plt.plot(z, Width_sim)
plt.plot(z,Width_teo, label='Analytical')