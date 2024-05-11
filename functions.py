import numpy as np
import matplotlib.pyplot as plt

# No lineal effects over the pulses
def Eff_NL(U0,z,t,gamma):
  U = U0*np.exp(1j*gamma*(np.abs(U0)**2)*1E-3*z)
  return U

# Lineal effects over the pulses
def Eff_L(U0,z,alph,beta_2,ts,N):
  F_U0 = np.fft.fft(U0)
  freq = np.fft.fftfreq(N,ts)
  F = F_U0*np.exp(2*1j*np.pi*np.pi*beta_2*freq*freq*z-alph*z/2)
  return F

# Split-step Fourier method
def StepFourier_NL_L(U0,h,z_tot,alph,beta_2,gamma,ts,N,t):
  U = U0
  p=0
  while(p<z_tot):
    F = Eff_L(U,h,alph,beta_2,ts,N)
    U = np.fft.ifft(F)
    U = Eff_NL(U,h,t,gamma)
    p+=h
  return U

def superposition_pulse_generator(P0, T0, df):
    E = np.sqrt(P0) * np.exp(-t * t / (2 * T0 ** 2))
    U0 = E * np.exp(1j * 2 * np.pi * df / 2 * t) + E * np.exp(-1j * 2 * np.pi * df / 2 * t)
    return U0

# Function to obtain frequency solution with beta 2 term
def F_disp2(FT_U0, z, beta_2, ts, N):
  frec = np.fft.fftfreq(N, ts)
  F = FT_U0 * np.exp(2 * 1j * np.pi * np.pi * beta_2 * frec * frec * z)
  return F

# Function to obtain frequency solution with beta 2 term and attenuation
def F_disp2_at(FT_U0, z, beta_2, alph, ts, N):
  frec = np.fft.fftfreq(N, ts)
  F = FT_U0 * np.exp(2 * 1j * np.pi * np.pi * beta_2 * frec * frec * z - alph * z / 2)
  return F

# Function to obtain frequency solution with beta 3 term
def F_disp3(FT_U0, z, beta_3, ts, N):
  frec = np.fft.fftfreq(N, ts)
  F = FT_U0 * np.exp(4 * 1j * (np.pi ** 3) * beta_3 * frec * frec * frec * z / 3)
  return F

# Function to obtain frequency solution with both beta 2 and beta 3 terms
def F_disp23(FT_U0, z, beta_2, beta_3, ts, N):
  frec = np.fft.fftfreq(N, ts)
  F = FT_U0 * np.exp(2 * 1j * np.pi * np.pi * beta_2 * frec * frec * z + 2 * 1j * np.pi * np.pi * beta_3 * frec * frec * frec * z / 3)
  return F

# Function to determine the width of the gaussian pulse
def GaussianWidth(G, t):
  Max = np.max(G)
  imax = np.argmax(G)
  i1 = 0
  i2 = 0
  aux = abs(G[0] - Max * (1 / np.e))
  for i in range(0, imax):
    if (abs(G[i] - Max * (1 / np.e))) < aux:
      aux = abs(G[i] - Max * (1 / np.e))
      i1 = i
  aux = abs(G[imax] - Max * (1 - 1 / np.e))
  for i in range(imax, len(G)):
    if (abs(G[i] - Max * (1 / np.e))) < aux:
      aux = abs(G[i] - Max * (1 / np.e))
      i2 = i
  return (t[i2] - t[i1]) / 2
