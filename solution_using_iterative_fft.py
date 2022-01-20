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

from sys import stdin
from sys import setrecursionlimit
import cmath
import math

MOD = 7340033
PI = math.acos(-1)

def fft(a, invert):
    """
        Description: This function computes the direct (or inverse) DFT of the given polynomial.
        Input:
            - a: A polynomial in coefficients representation.
            - invert: A flag that indicates whether the direct or the inverse DFT should be computed.
        Output:
            Direct (or inverse) DFT of the polynomial a.
    """
    n = len(a)
    i, j = 1, 0
    while(i < n):
        bit = n >> 1
        while(j & bit):
            j = j ^ bit
            bit = bit >> 1
        j = j ^ bit
        if(i < j):
            a[i], a[j] = a[j], a[i]
        i += 1
    length = 2
    while(length <= n):
        alpha = (2 * PI)/length * (-1 if invert else 1)
        wlen = complex(math.cos(alpha), math.sin(alpha)) #e**(i*(2pi/length)) = cos(2pi/length) + isin(2pi/length)
        for i in range(0, n, length):
            w = complex(1, 0)
            for j in range(0, length//2):
                u = a[i + j]; v = w * a[i + j + length//2]
                a[i + j] = u + v
                a[i + j + length//2] = u - v
                w = w * wlen
        length = length << 1
    if(invert):
        for i in range(n):
            a[i] = a[i]/n

def multiply(a, b):
    """
        Description: This function computes the multiplication of the given polynomials.
        Input:
            - a: A polynomial in coefficients representation.
            - b: A polynomial in coefficients representation.
        Output:
            The multiplication of a and b.
    """
    fa, fb = list(), list()
    # The two input vectors are copied into new vectos casting its values as complex numbers
    for i in range(len(a)):
        fa.append(complex(a[i]))
    for i in range(len(b)):
        fb.append(complex(b[i]))
    n = 1 # Size of the vector after multiplying the two polynomials.
    while(n < len(a) + len(b)):
        n = n << 1
    for _ in range(n - len(fa)):
        fa.append(complex(0, 0))
    for _ in range(n - len(fb)):
        fb.append(complex(0, 0))

    fft(fa, False)
    fft(fb, False)
    for i in range(n):
        fa[i] = fa[i]*fb[i]
    fft(fa, True)

    ans = [None for _ in range(n)]
    for i in range(n):
        ans[i] = round(fa[i].real)
    return ans

def C(n, K):
    """
        Description: This function computes the binomial coefficient of C(n, 0), C(n, 1), ... , C(n, K).
        Input:
            - n: An integer value greater than zero.
            - K: An integer value greater than zero and less or equal than n.
        Output:
            An array ans[0 .. K], where ans[k] = C(n, k).
    """
    ans = [0 for _ in range(K + 1)]
    ans[0] = 1; ans[K] = 1
    h = min(K, n//2)
    for k in range(1, h + 1):
        ans[k] = ((n - k + 1)*(ans[k - 1])//k)
        if(n - k <= K):
            ans[n - k] = ans[k]
    return ans

def construct_polynomial(c1, c2):
    """
        Description: This function constructs two polynomials according to the given number of points on each quadrant.
        Input:
            - c1: An integer value greater or equal to 0.
            - c2: An integer value greater or equal to 0.
        Output:
            A polynomial in coefficients representation.
    """
    K = min(c1, c2)
    a = C(c1, K)
    b = C(c2, K)
    p = [None for _ in range(K + 1)]
    for i in range(K + 1):
        p[i] = (a[i]*b[i])
    return p

def solve(coordinates, n):
    ans = list()
    p1 = construct_polynomial(coordinates[0], coordinates[2])
    p2 = construct_polynomial(coordinates[1], coordinates[3])
    p = multiply(p1, p2)
    for i in range(1, n//2 + 1):
        ans.append(p[i] % MOD)
    return ans

def main():
    """
        Description: This function handles the input and output.
        Input: None.
        Output: None.
    """
    lines = stdin.readlines()
    num_line = 0
    ntc = int(lines[num_line].strip()); num_line += 1
    tc = 0
    while(tc < ntc):
        n = int(lines[num_line].strip()); num_line += 1
        coordinates = [0, 0, 0, 0] # coordinates[i], 0 < i < 4, represents the number of coordinates located in the (i + 1)th quadrant
        for _ in range(n):
            x, y =  map(int, lines[num_line].strip().split()); num_line += 1
            if(x > 0 and y > 0):
                coordinates[0] += 1
            elif(x < 0 and y > 0):
                coordinates[1] += 1
            elif(x < 0 and y < 0):
                coordinates[2] += 1
            else:
                coordinates[3] += 1
        ans = solve(coordinates, n)
        print("Case %d:"%(tc + 1)); tc += 1
        print("%d"%(0), end = "")
        i = 0
        for j in range(1, n):
            if(j % 2 == 0):
                print(" %d"%(0), end = "")
            else:
                print(" %d"%(ans[i]), end = "")
                i += 1
        print()
    return 0

if __name__ == '__main__':
    main()