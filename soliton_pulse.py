
import functions

# Pulses parameters
T0 = 25         # Pulse width [ps]  f0 =40G
ts = 0.5        # Sampling time      V=10000G
n = 14
N = np.power(2,n)
T = N*ts        # Window width [ps]  fs~3G
t = np.arange(-N/2,N/2)*ts  # Time vector [ps]
freq = np.arange(-N/2,N/2)/(ts*N)

z_tot=50        # Total distance traveled by pulses
h=0.05

P=10                      # mW
gamma = 1.8               # 1/W*Km
alph = 0
Ld = 20                   # km
beta_2 = -T0**2/Ld

# Time response of gaussian and sinh pulse to z = 0 km
U0_sh = np.sqrt(1E3*abs(beta_2)/(gamma*T0**2))/np.cosh(t/T0)    # sinh pulse
U0 = np.sqrt(P)*np.exp(-t*t/(2*T0*T0))                          # gaussian pulse

# Frequency response of gaussian and sinh pulse to z = 0 km
F0_sh = np.fft.fftshift(np.fft.fft(U0_sh))                      # sinh pulse
F0 = np.fft.fftshift(np.fft.fft(U0))                            # gaussian pulse

# Time response of gaussian and sinh pulse to z = z_tot km
U_sh = StepFourier_NL_L(U0_sh,h,z_tot,alph,beta_2,gamma,ts,N,t) # sinh pulse
U = StepFourier_NL_L(U0,h,z_tot,alph,beta_2,gamma,ts,N,t)       # gaussian pulse

# Frequency response of gaussian and sinh pulse to z = z_tot km
F_sh = np.fft.fftshift(np.fft.fft(U_sh))                        # sinh pulse
F = np.fft.fftshift(np.fft.fft(U))                              # gaussian pulse

# Plot the frequency and time response of the Gaussian and sinh pulses when they have traveled 0 km and z_tot km.
fig, ax = plt.subplots(2,2)
fig.set_figheight(12)
fig.set_figwidth(14)

ax[0][0].plot(t,np.abs(U0)**2,label='Gaussian pulse $z=0$ km')
ax[0][0].set(xlim=[-500,500])
ax[0][0].plot(t,np.abs(U)**2,label='Gaussian pulse $z=$'+str(z_tot)+' km')
ax[0][0].set(xlim=[-500,500])
ax[0][0].set_title("Time response of gaussian pulse")
ax[0][0].set_ylabel("Amplitude")
ax[0][0].set_xlabel("Time [ps]")
ax[0][0].legend()

ax[1][0].plot(t,np.abs(U0_sh)**2,label = 'Sinh pulse $z=0$ km')
ax[1][0].set(xlim=[-500,500])
ax[1][0].plot(t,np.abs(U_sh)**2, label='Sinh pulse $z=$'+str(z_tot)+' km')
ax[1][0].set(xlim=[-500,500])
ax[1][0].set_title("Time response of sinh pulse")
ax[1][0].set_ylabel("Amplitude")
ax[1][0].set_xlabel("Time [ps]")
ax[1][0].legend()

ax[0][1].plot(freq,np.abs(F0)**2,label='Gaussian pulse spectrum $z=0$ km')
ax[0][1].set(xlim=[-0.05,0.05])
ax[0][1].plot(freq,np.abs(F)**2,label='Gaussian pulse spectrum $z=$'+str(z_tot)+' km')
ax[0][1].set(xlim=[-0.05,0.05])
ax[0][1].set_title("Frequency response of gaussian pulse")
ax[0][1].set_ylabel("Amplitude")
ax[0][1].set_xlabel("Frequency [GHz]")
ax[0][1].legend()

ax[1][1].plot(freq,np.abs(F0_sh)**2,label='Sinh pulse spectrum $z=0$ km')
ax[1][1].set(xlim=[-0.05,0.05])
ax[1][1].plot(freq,np.abs(F_sh)**2,label='Sinh pulse spectrum $z=$'+str(z_tot)+' km')
ax[1][1].set(xlim=[-0.05,0.05])
ax[1][1].set_title("Frequency response of sinh pulse")
ax[1][1].set_ylabel("Amplitude")
ax[1][1].set_xlabel("Frequency [GHz]")
ax[1][1].legend()