#!/usr/bin/python
# -*- coding: utf-8 -*-

# Nicolas Bigaouette
# Fall 2010

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy

def fft(t, s, resize_NFFT = True, resize_NFFT_UP = False):
    # t                 : Time  array
    # s                 : Signal array
    # resize_NFFT       : Resize signal and time array for optimal FFT calculation?
    # resize_NFFT_UP    : If resize is to be done, increase the number of points (or decrease it)?
    # Returns:
    # S                 : TF(s)
    # S_Amplitude       : TF amplitude (|S|^2)
    # f                 : Frequency array (time^-1)
    # omega             : Angular frequency

    # Sampling frequency (= dt^-1)
    dt = t[1] - t[0]
    Fs = 1.0 / dt

    NFFT  = len(s)
    if (resize_NFFT):
        # Makes FFT algorithm faster
        # Closest power of two smaller or equal to len(x)
        print "Changing NFFT from", NFFT, "to",
        if (resize_NFFT_UP):
            NFFT = int(2**(numpy.ceil(numpy.log(NFFT) / numpy.log(2))))
        else:
            NFFT = int(2**(numpy.floor(numpy.log(NFFT) / numpy.log(2))))
        #NFFT = 2**6*int(2**(numpy.ceil(numpy.log(NFFT) / numpy.log(2))))
        print NFFT

    #assert(NFFT <= len(s))

    NFFTd = float(NFFT)

    # Frequency vector
    f     = numpy.arange(-NFFTd / 2.0, NFFTd / 2.0) * Fs / NFFTd
    omega = 2.0 * numpy.pi * f

    # S = TF(s)
    S = numpy.fft.fftshift(numpy.fft.fft(s[0:NFFT], NFFT))

    # Keep spectrum for frequencies lower then 1
    #i = numpy.where(omega < 1.0)
    #S = S[i]
    #f = f[i]
    #omega = omega[i]

    #S_Amplitude = (S.real*S.real) + (S.imag*S.imag)
    S_Amplitude = abs(S)**2

    return [S, S_Amplitude, f, omega]
#


def test_fft():
    import matplotlib.pyplot as plt

    # Number of periods
    nT = 10.0

    # Periods
    #T  = numpy.array([5.0, 10.0, 15.0, 20.0])
    #f  = 1.0 / T    # Frequency
    #o  = 2.0 * numpy.pi * f

    # Angulare frequencies
    o  = numpy.array([0.3, 0.35, 0.4, 0.6, 1.25])
    f = o / (2.0 * numpy.pi)
    T = 1.0 / f

    nt = 5000     # Nb of time steps
    tmax = nT * T.max()
    t = numpy.linspace(0.0, tmax, nt)
    dt = t[1] - t[0]

    # Build signal
    s = numpy.zeros(nt, dtype=numpy.complex128)
    for i in xrange(len(T)):
        s += numpy.cos(o[i] * t) + 1.0j*numpy.sin(o[i] * t)

    # Add noise
    noise_amplitude = 0.0
    static = noise_amplitude * (numpy.random.random_sample((nt))-0.5)
    s += (1.0 + 1.0j) * static

    # Calculate FFT
    [FT, FTa, frequencies, angular_frequencies] = fft(t, s)
    FTa_min = (FTa / FTa.max()).min()

    print "Periods:", T
    print "Frequencies:", f
    print "Angular Frequencies:", o

    print "Spectrum's frequency range given by number of points."
    print "Spectrum's frequency precision given by duration."

    omega_range     = (2.0*numpy.pi) / dt
    omega_precision = omega_range / nt
    print "Spectrum's angular frequency range: 2 pi / (nt dt) =", omega_range,     "(" + str(angular_frequencies[-1] - angular_frequencies[0]) + ")"
    print "Spectrum's angular frequency precision:",              omega_precision, "(" + str(angular_frequencies[1]  - angular_frequencies[0]) + ")"

    fig = plt.figure()
    colors = ['b', 'r', 'm', 'c', 'g', 'y', 'k']

    sp1 = plt.subplot(211)
    plt.plot(t, s.real, lw=2, label='Real')
    plt.plot(t, s.imag, lw=2, label='Imaginary')
    plt.xlabel("t [time]")
    plt.legend()

    sp2 = plt.subplot(212)
    plt.plot(angular_frequencies, FTa / FTa.max(), '--xk', lw=2, label=r"|FFT|$^2$")
    for i in xrange(len(o)):
        plt.plot([ o[i],  o[i]], [0.9*FTa_min, 1.1],  '-' + colors[i%len(colors)], lw=2, label=r'$\omega = ' + str('%.4g' %  o[i]) + '$')
        plt.plot([-o[i], -o[i]], [0.9*FTa_min, 1.1], '--' + colors[i%len(colors)], lw=2, label=r'$\omega = ' + str('%.4g' % -o[i]) + '$')
    plt.xlabel(r"$\omega$ [time$^{-1}$]")
    plt.ylabel(r"|FFT|$^2$")
    sp2.set_yscale('log')
    sp2.set_xlim((-1.5*o.max(), 1.5*o.max()))
    sp2.set_ylim((FTa_min, 1.1))
    plt.legend()

    plt.show()
#

if __name__ == "__main__":
    test_fft()

