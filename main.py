#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
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

        self.dwt_array = list()
        self.count = 0

    def y(self, x):
        return np.cos(2 * x) + np.sin(5 * x)

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

    def fwt(self, a):
        if len(a) == 1:
            return a
        c1 = list()
        c2 = list()
        half = len(a) // 2
        for i in range(0, half):
            c1.append(a[i] + a[i + half])
        for i in range(half, len(a)):
            c2.append(a[i - half] - a[i])
        c1 = self.fwt(c1)
        c2 = self.fwt(c2)
        return list(c1 + c2)

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

    def draw_signals(self, grid):
        plt.subplot(grid)
        title = "Default signal"
        plt.title(title)
        value = self.y_array
        abs_val = self.get_abs(value)
        plt.vlines(self.x_array, 0, abs_val)

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

    def draw_fwt(self, fignum):
        plt.subplot(fignum)
        plt.title("Fast Walsh Transform")
        fwt_array = self.fwt(self.y_array)
        self.dwt_array = fwt_array

        abs_y = list()
        #  print(fwt_array)
        for num in fwt_array:
            abs_y.append(np.absolute(num.real))
        plt.vlines(self.x_array, 0, abs_y)

    def draw_fwt_inv(self, fignum):
        plt.subplot(fignum)
        plt.title("Inverse Fast Walsh Transform")
        fwt_array = self.fwt(self.dwt_array)
        arr = np.array(fwt_array)

        fwt_real = list()
        abs_y = list()
        for num in arr:
            fwt_real.append(num)
            abs_y.append(np.absolute(num))
        plt.vlines(self.x_array, 0, abs_y)

    def draw(self):
        gs = gridspec.GridSpec(2, 3)
        self.draw_signals(gs[1, 0])

        self.draw_dwt(gs[0, 1])
        self.draw_idwt(gs[1, 1])

        self.draw_fwt(gs[0, 2])
        self.draw_fwt_inv(gs[1, 2])
        plt.show()


if __name__ == '__main__':
    calc = Calculator()
    calc.draw()
