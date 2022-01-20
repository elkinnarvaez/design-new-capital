"""
    Author: Elkin Jadier Narvaez Paz - 8943358

    Como miembro de la comunidad académica de la Pontifica Universidad Javeriana Cali me comprometo a seguir los más altos estándares
    de integridad académica.

    References:
        - Jia Yan-Bin. Polynomial Multiplication and Fast Fourier Transform. Sep 17, 2020. - FFT
        - Geeks For Geeks
            FFT: url: https://www.geeksforgeeks.org/fast-fourier-transformation-poynomial-multiplication/
        - CP Algorithms
            FFT: https://cp-algorithms.com/algebra/fft.html
            Binomial Coefficients: https://cp-algorithms.com/combinatorics/binomial-coefficients.html#toc-tgt-9
            Modular Multiplicative Inverse: https://cp-algorithms.com/algebra/module-inverse.html#toc-tgt-2 - http://e-maxx.ru/algo/reverse_element
"""

# AC for PR-A, PR-B, PR-C and PR-D

from sys import stdin
import cmath
import math

PI = math.acos(-1)
MAXN = 10**5

mod = 7340033
root = 5 # (1 << 20)-th root of unity 
root_1 = 4404020 # inverse of root modulo mod
root_pw = 1 << 20

fact = None
fact_inv = None
inv = None

a = [0 for _ in range(2**17 + 1)]
b = [0 for _ in range(2**17 + 1)]
na = None; nb = None
n = None

def binpow(a, b):
    ans = None
    if(b == 0): ans = 1
    else:
        ans = binpow(a, b/2)
        if (b % 2 == 0): ans = ans * ans
        else: ans = ans * ans * a
    return ans

def expmod(a, b):
    """
        Description: This function computes pow(a, b) % mod
        Input:
            - a: An integer value
            - b: An integer value
        Output:
            pow(a, b) % mod
    """
    ans = 1
    while(b > 0):
        if(b & 1):
            ans = ans * a % mod
        a = a * a % mod
        b = b >> 1
    return ans

def invmod(a):
    """
        Description: This function calculates inv(a) % mod. mod must be a prime number.
        Input:
            - a: An integer value
        Output:
            inv(a) % mod
    """
    return expmod(a, mod - 2)

def prec_fact():
    """
        Description: This function precomputes the factorial values (and its inverse) in the range 0 .. MAXN
        Input: None
        Output: None
    """
    global fact, fact_inv, inv
    fact = [None for _ in range(MAXN + 1)]
    fact_inv = [None for _ in range(MAXN + 1)]
    inv = [None for _ in range(MAXN + 1)]
    fact[0] = 1; fact[1] = 1
    fact_inv[0] = 1; fact_inv[1] = 1
    inv[1] = 1
    for i in range(2, MAXN + 1):
        inv[i] = (mod - (mod//i) * inv[mod%i] % mod) % mod
        fact[i] = fact[i - 1] * i % mod
        fact_inv[i] = fact_inv[i - 1] * inv[i] % mod


def fft(p, invert):
    """
        Description: This function computes the direct (or inverse) DFT of the given polynomial.
        Input:
            - a: A polynomial in coefficients representation.
            - invert: A flag that indicates whether the direct or the inverse DFT should be computed.
            - p: A flag that indicates what polynomial is being computed. b = True indicates that the polynomial a is computed. b = False indicated that the polynomial b in being computed. 
        Output:
            Direct (or inverse) DFT of the polynomial a.
    """
    global a, b
    i, j = 1, 0
    while(i < n):
        bit = n >> 1
        while(j & bit):
            j = j ^ bit
            bit = bit >> 1
        j = j ^ bit
        if(i < j):
            if(p): a[i], a[j] = a[j], a[i]
            else: b[i], b[j] = b[j], b[i]
        i += 1
    length = 2
    while(length <= n):
        wlen = root_1 if invert else root
        i = length
        while(i < root_pw):
            wlen = wlen * wlen % mod
            i = i << 1
        for i in range(0, n, length):
            w = 1
            for j in range(0, length//2):
                if(p):
                    u = a[i + j]; v = w * a[i + j + length//2] % mod
                    a[i + j] = u + v if u + v < mod else u + v - mod
                    a[i + j + length//2] = u - v if u - v >= 0 else u - v + mod
                else:
                    u = b[i + j]; v = w * b[i + j + length//2] % mod
                    b[i + j] = u + v if u + v < mod else u + v - mod
                    b[i + j + length//2] = u - v if u - v >= 0 else u - v + mod
                w = w * wlen % mod  
        length = length << 1
    if(invert):
        n_1 = inv[n]
        for i in range(n):
            if(p):
                a[i] = a[i] * n_1 % mod
            else:
                b[i] = b[i] * n_1 % mod

def multiply():
    """
        Description: This function computes the multiplication of the polynomials a and b.
        Input: None. Both polynomials are global variables.
        Output: None. The answer is stored in a, which corresponds to a global variable.
    """
    global n, a, b
    n = 1 # Size of the vector after multiplying the two polynomials.
    while(n < na + nb):
        n = n << 1
    # Need complete with zeroes
    for i in range(na, n + 1):
        a[i] = 0
    for i in range(nb, n + 1):
        b[i] = 0
    fft(True, False)
    fft(False, False)
    for i in range(n):
        a[i] = a[i]*b[i]
    fft(True, True)

def C(n, k):
    """
        Description: This function calculates nCk.
        Input:
            - n: An integer value.
            - k: An integer value, such that 0 <= k <=n.
        Output:
            nCk
    """
    return fact[n] * fact_inv[k] * fact_inv[n - k] % mod

def construct_polynomial(c1, c2, p):
    """
        Description: This function constructs a polynomial according to the given number of points on each quadrant.
        Input:
            - c1: An integer value greater or equal to 0.
            - c2: An integer value greater or equal to 0.
            - p: A flag that indicates what polynomial is being computed. b = True indicates that the polynomial a is computed. b = False indicates that the polynomial b in being computed. 
        Output: None. Both polynomials are global variables.
    """
    global a, b, na, nb
    K = min(c1, c2)
    if(p): na = K + 1
    else: nb = K + 1
    for k in range(K + 1):
        if(p): a[k] = (C(c1, k)*C(c2, k))
        else: b[k] = (C(c1, k)*C(c2, k))

def solve(coordinates, num_coordinates):
    ans = list()
    construct_polynomial(coordinates[0], coordinates[2], True)
    construct_polynomial(coordinates[1], coordinates[3], False)
    multiply()
    for i in range(1, num_coordinates//2 + 1):
        ans.append(a[i] % mod)
    return ans

def main():
    """
        Description: This function handles the input and output.
        Input: None.
        Output: None.
    """
    prec_fact()
    lines = stdin.readlines()
    num_line = 0
    ntc = int(lines[num_line].strip()); num_line += 1
    tc = 0
    while(tc < ntc):
        num_coordinates = int(lines[num_line].strip()); num_line += 1
        coordinates = [0, 0, 0, 0] # coordinates[i], 0 < i < 4, represents the number of coordinates located in the (i + 1)th quadrant
        for _ in range(num_coordinates):
            x, y =  map(int, lines[num_line].strip().split()); num_line += 1
            if(x > 0 and y > 0):
                coordinates[0] += 1
            elif(x < 0 and y > 0):
                coordinates[1] += 1
            elif(x < 0 and y < 0):
                coordinates[2] += 1
            else:
                coordinates[3] += 1
        ans = solve(coordinates, num_coordinates)
        print("Case %d:"%(tc + 1)); tc += 1
        print("%d"%(0), end = "")
        i = 0
        for j in range(1, num_coordinates):
            if(j % 2 == 0):
                print(" %d"%(0), end = "")
            else:
                print(" %d"%(ans[i]), end = "")
                i += 1
        print()
    return 0

if __name__ == '__main__':
    main()