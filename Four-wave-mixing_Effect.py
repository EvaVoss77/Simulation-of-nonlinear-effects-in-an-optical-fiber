'''Propagation of two pulse sequences is studied at wavelengths λ1 and λ2 separated by 50 GHz, for z_tot km of SMF and TW-RS fibers, neglecting the effect of attenuation. Four-wave mixing generation is observed for cases where the input power is P01=1, P02=10, and P03=100 mW for each wavelength.'''
ts = 0.5       # Sampling time      V=10000G
n = 14
N = np.power(2, n)
T = N * ts       # Window width [ps]  fs~3G
t = np.arange(-N/2, N/2) * ts  # Temporal vector [ps]
frec = np.arange(-N/2, N/2) / (ts * N)

# Pulse parameters
z_tot = 10
h = 0.05

gamma = 1.2      # Non-linear coefficient
alph = 0         # Attenuation coefficient
Ld = 20          # Characteristic length
beta_2 = -21

T0 = 25         # Pulse width [ps]
df = 0.1       # Pulse separation

P01 = 1         # 1 mW
P02 = 10        # 10 mW
P03 = 100       # 100 mW

# Superposition of both gaussian pulses for each power value
U01 = superposition_pulse_generator(P01, T0, df)
U02 = superposition_pulse_generator(P02, T0, df)
U03 = superposition_pulse_generator(P03, T0, df)

U1 = StepFourier_NL_L(U01, h, z_tot, alph, beta_2, gamma, ts, N, t)
U2 = StepFourier_NL_L(U02, h, z_tot, alph, beta_2, gamma, ts, N, t)
U3 = StepFourier_NL_L(U03, h, z_tot, alph, beta_2, gamma, ts, N, t)

# Frequency responses when pulses have traveled 0 and z_tot km
F01 = np.fft.fftshift(np.fft.fft(U01))
F02 = np.fft.fftshift(np.fft.fft(U02))
F03 = np.fft.fftshift(np.fft.fft(U03))

F1 = np.fft.fftshift(np.fft.fft(U1))
F2 = np.fft.fftshift(np.fft.fft(U2))
F3 = np.fft.fftshift(np.fft.fft(U3))

# Plot for each pulse at z=0 km and z=z_tot km
fig, ax = plt.subplots(3, 2, figsize=(15, 15))

# Plot for P01
ax[0, 0].plot(t, abs(U01) ** 2, label='z=0')
ax[0, 0].plot(t, abs(U1) ** 2, label='z=' + str(z_tot) + ' km')
ax[0, 0].set(xlim=[-100, 100], xlabel='Time [ps]', ylabel='Intensity', title='Temporal Response P01')
ax[0, 0].legend()

ax[0, 1].semilogy(frec, abs(F01) ** 2, label='z=0')
ax[0, 1].semilogy(frec, abs(F1) ** 2, label='z=' + str(z_tot) + ' km')
ax[0, 1].set(xlabel='Frequency [GHz]', ylabel='Intensity', title='Frequency Response - P01')
ax[0, 1].legend()

# Plot for P02
ax[1, 0].plot(t, abs(U02) ** 2, label='z=0')
ax[1, 0].plot(t, abs(U1) ** 2, label='z=' + str(z_tot) + ' km')
ax[1, 0].set(xlim=[-100, 100], xlabel='Time [ps]', ylabel='Intensity', title='Temporal Response P02')
ax[1, 0].legend()

ax[1, 1].semilogy(frec, abs(F02) ** 2, label='z=0')
ax[1, 1].semilogy(frec, abs(F1) ** 2, label='z=' + str(z_tot) + ' km')
ax[1, 1].set(xlabel='Frequency [GHz]', ylabel='Intensity', title='Frequency Response - P02')
ax[1, 1].legend()

# Plot for P03
ax[2, 0].plot(t, abs(U03) ** 2, label='z=0')
ax[2, 0].plot(t, abs(U1) ** 2, label='z=' + str(z_tot) + ' km')
ax[2, 0].set(xlim=[-100, 100], xlabel='Time [ps]', ylabel='Intensity', title='Temporal Response P03')
ax[2, 0].legend()

ax[2, 1].semilogy(frec, abs(F03) ** 2, label='z=0')
ax[2, 1].semilogy(frec, abs(F1) ** 2, label='z=' + str(z_tot) + ' km')
ax[2, 1].set(xlabel='Frequency [GHz]', ylabel='Intensity', title='Frequency Response - P03')
ax[2, 1].legend()

plt.tight_layout()
plt.show()