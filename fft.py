#!/usr/bin/env python2
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

    # Frequencies
    f = numpy.array([1.0, 5.0, 2.5])
    T = 1.0 / f
    o  = 2.0 * numpy.pi * f

    # Angulare frequencies
    #o  = numpy.array([0.3, 0.35, 0.4, 0.6, 1.25])
    #f = o / (2.0 * numpy.pi)
    #T = 1.0 / f

    nt = 5000     # Nb of time steps
    tmax = nT * T.max()
    time = numpy.linspace(0.0, tmax, nt)
    dt = time[1] - time[0]

    # Build signal
    signal = numpy.zeros(nt, dtype=numpy.complex128)
    for i in xrange(len(T)):
        # Real signal: negative and positive frequencies
        #signal += numpy.cos(o[i] * time)
        # Complex signal: only positive frequencies
        signal += numpy.cos(o[i] * time) + 1.0j*numpy.sin(o[i] * time)
    #

    # Add noise
    noise_amplitude = 0.0
    static = noise_amplitude * (numpy.random.random_sample((nt))-0.5)
    signal += (1.0 + 1.0j) * static

    # Calculate FFT
    [FT, FTa, frequencies, angular_frequencies] = fft(time, signal, resize_NFFT = False)

    print "Periods: T =", T
    print "Frequencies: f =", f
    print "Angular Frequencies: ω =", o

    print "Spectrum's frequency range given by number of points."
    print "Spectrum's frequency precision given by duration."

    omega_range     = (2.0*numpy.pi) / dt
    omega_precision = omega_range / nt
    print "Spectrum's angular frequency range:     2 pi / (dt)    =", omega_range,     "(" + str(angular_frequencies[-1] - angular_frequencies[0]) + ")"
    print "Spectrum's angular frequency precision: 2 pi / (nt dt) =", omega_precision, "(" + str(angular_frequencies[1]  - angular_frequencies[0]) + ")"
    print "Spectrum's frequency range:     1 / (dt)    =", omega_range/(2.0*numpy.pi),     "(" + str((angular_frequencies[-1] - angular_frequencies[0])/(2.0*numpy.pi)) + ")"
    print "Spectrum's frequency precision: 1 / (nt dt) =", omega_precision/(2.0*numpy.pi), "(" + str((angular_frequencies[1]  - angular_frequencies[0])/(2.0*numpy.pi)) + ")"

    fig = plt.figure()
    colors = ['b', 'r', 'm', 'c', 'g', 'y', 'k']

    sp1 = plt.subplot(211)
    plt.plot(time, signal.real, lw=2, label='Real')
    plt.plot(time, signal.imag, lw=2, label='Imaginary')
    plt.xlabel("t [time]")
    plt.legend()

    sp2 = plt.subplot(212)
    plt.plot(frequencies, FTa / FTa.max(), '--xk', lw=2, label=r"|FFT|$^2$")
    for i in xrange(len(o)):
        plt.axvline( f[i],  linestyle = '-',  color = colors[i%len(colors)], linewidth=2, label=r'$f = ' + str('%.4g' %  f[i]) + '$')
        plt.axvline(-f[i],  linestyle = '--', color = colors[i%len(colors)], linewidth=2, label=r'$f = ' + str('%.4g' % -f[i]) + '$')
    plt.xlabel(r"frequencies [time$^{-1}$]")
    plt.ylabel(r"|FFT|$^2$")
    sp2.set_yscale('log')
    sp2.set_xlim((-1.5*f.max(), 1.5*f.max()))
    sp2.set_ylim((1.0e-12, 1.1))
    plt.legend()

    plt.show()
#

if __name__ == "__main__":
    test_fft()

