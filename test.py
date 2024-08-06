import cmath
from math import log
import numpy as np

def bit_reverse_vec(values):
    """Reverses list by reversing the bits of the indices.

    Reverse indices of the given list.
    For example, reversing the list [0, 1, 2, 3, 4, 5, 6, 7] would become
    [0, 4, 2, 6, 1, 5, 3, 7], since 1 = 0b001 reversed is 0b100 = 4,
    3 = 0b011 reversed is 0b110 = 6.

    Args:
        values (list): List of values to be reversed. Length of list must be a power of two. 

    Returns:
        The reversed list based on indices.
    """
    result = [0] * len(values)
    for i in range(len(values)):
        result[i] = values[reverse_bits(i, int(log(len(values), 2)))]
    return result

def reverse_bits(value, width):
    """Reverses bits of an integer.

    Reverse bits of the given value with a specified bit width.
    For example, reversing the value 6 = 0b110 with a width of 5
    would result in reversing 0b00110, which becomes 0b01100 = 12.

    Args:
        value (int): Value to be reversed.   
        width (int): Number of bits to consider in reversal.

    Returns:
        The reversed int value of the input.
    """
    binary_val = '{:0{width}b}'.format(value, width=width)
    return int(binary_val[::-1], 2)

def embedding(coeffs):
    """Computes a variant of the canonical embedding on the given coefficients.

    Computes the canonical embedding which consists of evaluating a given polynomial at roots of unity
    that are indexed 1 (mod 4), w, w^5, w^9, ...
    The evaluations are returned in the order: w, w^5, w^(5^2), ...

    Args:
        coeffs (list): List of complex numbers to transform.

    Returns:
        List of transformed coefficients.
    """
    num_coeffs = len(coeffs)
    result = bit_reverse_vec(coeffs)
    log_num_coeffs = int(log(num_coeffs, 2))
    num_slots = fft_length // 4
    root = np.exp((-1j * np.pi) / fft_length)
    roots_of_unity = [(root) ** (5**j % (2 * fft_length)) for j in range(0,fft_length// 2)]

    rot_group = [1] * num_slots

    for i in range(1, num_slots):
        rot_group[i] = (5 * rot_group[i - 1]) % fft_length

    for logm in range(1, log_num_coeffs + 1):
        idx_mod = 1 << (logm + 2)
        gap = fft_length // idx_mod
        for j in range(0, num_coeffs, (1 << logm)):
            for i in range(1 << (logm - 1)):
                index_even = j + i
                index_odd = j + i + (1 << (logm - 1))

                rou_idx = (rot_group[i] % idx_mod) * gap
                omega_factor = roots_of_unity[rou_idx] * result[index_odd]

                butterfly_plus = result[index_even] + omega_factor
                butterfly_minus = result[index_even] - omega_factor

                result[index_even] = butterfly_plus
                result[index_odd] = butterfly_minus

    return result

def fft(x):
    """
    Compute the Fast Fourier Transform of a 1D list of complex numbers x.

    Parameters:
    x (list of complex): The input time-domain signal.

    Returns:
    list of complex: The FFT of the input signal.
    """
    N = len(x)
    if N <= 1:
        return x

    # Separate even and odd indexed elements
    even = fft(x[0::2])
    odd = fft(x[1::2])

    # Combine results
    T = [cmath.exp(-2j * cmath.pi * k / N) * odd[k] for k in range(N // 2)]
    return [even[k] + T[k] for k in range(N // 2)] + [even[k] - T[k] for k in range(N // 2)]


# Example usage
fft_length = 8
X = [3+4j, 2-1j, 2+1j, 3-4j]  # Example FFT coefficients
# x = fft(X)
x = embedding(X)
print(list(map(lambda x: x/32, x)))
