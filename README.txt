This is a small wrapper around Numpy's FFT functions.

Historically, FFT functions takes the signal only as argument,
not the time array. The frequencies associated with the
discrete Fourier Transform are not calculated and left
to the user.

fftshift() also needs to be called on the resulting
Fourier Transform to fhift the zero-frequency component
to the center of the spectrum.

This package does it all: given a signal and its
time array, the discrete Fourier Transform is calculated
and returned, with the corresponding frequencies and angular
frequencies. The signal array can be automatically resized
to a power of 2 for better efficiency of the underlying
FFT algorithm.

For testing purpose, running directly will plot a test
signal and its discrete Fourier Transform:
$ python ./fft.py

Usage:
import fft
[...]
[FT, FTa, frequencies, angular_frequencies] = fft.fft(time, signal)

License: GPL v3
