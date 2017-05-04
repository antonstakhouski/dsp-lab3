#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import cmath
import math
import matplotlib.gridspec as gridspec


class Calculator:
    def __init__(self):
        self.n = 32
        self.start = 0.0
        self.end = 2 * math.pi
        self.step = (self.end - self.start) / self.n

        self.x_array = np.arange(self.start, self.end, self.step)
        self.y_array = self.y(self.x_array)

        self.dft_array = list()
        self.dwt_array = list()
        self.count = 0

    def y(self, x):
        return np.cos(2 * x) + np.sin(5 * x)

    def dft(self, k):
        m = 0
        c = 0
        w = cmath.exp(cmath.sqrt(-1) * 2 * cmath.pi / self.n)
        while(m < self.n):
            c += self.y_array[m] * w ** (k * m)
            m += 1
        return (1 / self.n) * c

    def dwt(self, k):
        m = 0
        c = 0
        while(m < self.n):
            c += self.y_array[m] * self.wal(k, m)
            m += 1
        return c

    def wal(self, s, k):
        exp = 0
        n = int(math.log(self.n, 2))
        s = self.gray(s)
        binary = bin(k)[2:]
        blen = len(binary)
        if(blen < n):
            binary = "0" * (n - blen) + binary
        for i in range(0, n):
            exp += int(s[i]) * int(binary[n - 1 - i])
        res = 1/math.sqrt(self.n) * (-1) ** exp
        return res

    def gray(self, s):
        binary = bin(s)[2:]
        n = int(math.log(self.n, 2))
        blen = len(binary)
        if(blen < n):
            binary = "0" * (n - blen) + binary
        gray = ""
        for i in range(1, n):
            gray = str(int(binary[-i]) ^ int(binary[-i - 1])) + gray
        gray = binary[0] + gray
        return gray

    def fft_dit(self, a, direct):
        if len(a) == 1:
            return a
        a_even = list()
        a_odd = list()
        i = 0
        while (i < len(a)):
            if i % 2 == 0:
                a_even.append(a[i])  # четный
            else:
                a_odd.append(a[i])
            i += 1
        b_even = self.fft_dit(a_even, direct)
        b_odd = self.fft_dit(a_odd, direct)

        self.count += 1

        argument = 2 * cmath.pi / len(a)
        wn = cmath.cos(argument) + direct * cmath.sqrt(-1) * cmath.sin(argument)

        w = np.complex(1)

        y = [np.complex(0)] * len(a)
        j = 0
        while(j < len(a) // 2):
            y[j] = b_even[j] + w * b_odd[j]
            y[j + len(a) // 2] = b_even[j] - w * b_odd[j]
            w *= wn
            j += 1
        return y

    def t(self, x):
        print(type(x))

    def idwt(self, k):
        n = int(math.log(self.n, 2))
        binary = bin(k)[2:]
        blen = len(binary)
        if(blen < n):
            binary = "0" * (n - blen) + binary
        exp1 = 0
        for s in range(0, self.n):
            sg = self.gray(s)
            exp = 0
            for i in range(0, n):
                exp += int(sg[i]) * int(binary[n - 1 - i])
            exp1 += self.dwt_array[s] * (-1) ** exp
        res = 1/math.sqrt(self.n) * exp1
        return res

    def idft(self, m):
        k = 0
        x = 0
        w = cmath.exp(cmath.sqrt(-1) * 2 * cmath.pi / self.n)
        while(k < self.n):
            x += self.dft_array[k] * w ** (- (k * m))
            k += 1
        return x

    def draw_signals(self, grid):
        plt.subplot(grid)
        title = "Default signal"
        plt.title(title)
        value = self.y_array
        abs_val = self.get_abs(value)
        plt.vlines(self.x_array, 0, abs_val)

    def draw_dft(self, fignum):
        plt.subplot(fignum)
        plt.title("Discrete Fourier Transform")
        self.dft_array.clear()
        i = 0
        while i < self.n:
            cnum = self.dft(i)
            self.dft_array.append(cnum)
            i += 1
        arr = np.array(self.dft_array)
        abs_y = np.absolute(arr)
        plt.vlines(self.x_array, 0, abs_y)

    def draw_dwt(self, fignum):
        plt.subplot(fignum)
        plt.title("Discrete Walsh Transform")
        self.dwt_array.clear()
        i = 0
        while i < self.n:
            cnum = self.dwt(i)
            self.dwt_array.append(cnum)
            i += 1
        arr = np.array(self.dwt_array)
        abs_y = np.absolute(arr)
        plt.vlines(self.x_array, 0, abs_y)

    def draw_idft(self, fignum):
        plt.subplot(fignum)
        plt.title("Inverse Discrete Fourier Transform")
        idft_array = list()
        i = 0
        while i < self.n:
            idft_array.append(self.idft(i))
            i += 1
        arr = np.array(idft_array)
        abs_y = np.absolute(arr)
        plt.vlines(self.x_array, 0, abs_y)

    def draw_idwt(self, fignum):
        plt.subplot(fignum)
        plt.title("Inverse Discrete Walsh Transform")
        idwt_array = list()
        i = 0
        while i < self.n:
            idwt_array.append(self.idwt(i))
            i += 1
        arr = np.array(idwt_array)
        abs_y = np.absolute(arr)
        plt.vlines(self.x_array, 0, abs_y)

    def get_abs(self, lst):
        abs_lst = list()
        for num in lst:
            abs_lst.append(np.absolute(num.real))
        return abs_lst

    def draw_fft(self, fignum):
        plt.subplot(fignum)
        plt.title("Fast Fourier Transform")
        fft_array = self.fft_dit(self.y_array, -1)
        self.dft_array = fft_array

        abs_y = list()
        for num in fft_array:
            abs_y.append(np.absolute(num.real))
        plt.vlines(self.x_array, 0, abs_y)

    def draw_fft_inv(self, fignum):
        plt.subplot(fignum)
        plt.title("Inverse Fast Fourier Transform")
        fft_array = self.fft_dit(self.dft_array, 1)
        arr = np.array(fft_array)

        fft_real = list()
        abs_y = list()
        for num in arr:
            fft_real.append(num)
            abs_y.append(np.absolute(num))
        plt.vlines(self.x_array, 0, abs_y)

    def draw_corr_conv(self, corr_conv, titles, grid, pos):
        for i in range(0, len(titles)):
            plt.subplot(grid[i, pos])
            plt.title(titles[i])
            value = corr_conv[i]
            abs_val = self.get_abs(value)
            plt.vlines(self.x_array, 0, abs_val)

    def fft_corr_conv(self):
        cy = self.fft_dit(self.y_array, -1)
        cz = self.fft_dit(self.z_array, -1)
        corr = list()
        conv = list()
        for i in range(0, self.n):
            corr.append(np.conjugate(cy[i]) * cz[i])
            conv.append(cy[i] * cz[i])
        corr = self.fft_dit(corr, 1)
        conv = self.fft_dit(conv, 1)
        return (corr, conv)

    def draw_all_corr_conv(self, grid, pos):
        ftype = ["Simple ", "FFT "]
        flst = [self.corr_conv(), self.fft_corr_conv()]
        for i in range(0, len(ftype)):
            titles = [ftype[i] + "Correlation", ftype[i] + "Convolution"]
            corr_conv = flst[i]
            self.draw_corr_conv(corr_conv, titles, grid, pos + i)

    def corr_conv(self):
        corr = list()
        conv = list()
        for m in range(0, self.n):
            corr_sum = 0
            conv_sum = 0
            for h in range(0, self.n):
                corr_sum += self.y_array[h] * self.z_array[(m + h) % self.n]
                conv_sum += self.y_array[h] * self.z_array[m - h]
            mul = 1 / self.n
            corr.append(mul * corr_sum)
            conv.append(mul * conv_sum)
        return (corr, conv)

    def draw(self):
        gs = gridspec.GridSpec(2, 3)
        self.draw_signals(gs[0, 0])
        self.draw_dwt(gs[0, 1])
        self.draw_idwt(gs[1, 1])
        plt.show()


if __name__ == '__main__':
    calc = Calculator()
    calc.draw()
