#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy  as np
from numpy import exp, pi, real, imag, angle, cos, sin, inf
import matplotlib.pyplot as plt
from numpy.fft import fft
from scipy.integrate import quad

from matplotlib import rc, rcParams
rc('font', family='Ubuntu')
#rcParams['toolbar'] = 'None'


def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

I1 = -5
I2 = 5

#Начальные данные
A = 1
alpha = 4
delta = 2/alpha
tau = 1/alpha
j=1j

R1 = 16.875
R2 = 4.32e3
#Функция Хевисайда
def u(t):
	if t > 0:
		return 1
	elif t == 0:
		return 1/2
	else:
		return 0

#Функция сигнала
def s(t):
	return A*exp(-alpha*(abs(t)-tau))*u(abs(t)-tau)-A*exp(-alpha*(tau+delta))*exp(-alpha*(abs(t)-2*tau-delta))*u(abs(t)-tau-delta)#(A*exp(-alpha*(t-tau)))*(u(t-tau)-u(t-(tau+delta)))+(A*exp(alpha*(t+tau)))*(u(t+(tau+delta))-u(t+tau))
	
def S(f):
	return (4*cos(pi*f)-pi*f*sin(pi*f) - 4*exp(-4)*(2*cos(3*pi*f)-pi*f*sin(3*pi*f)))/(3+pi**2*f**2)

def W(f):
	return abs((4*cos(pi*f)-pi*f*sin(pi*f) - 4*exp(-4)*(2*cos(3*pi*f)-pi*f*sin(3*pi*f)))/(3+pi**2*f**2))**2

def Rev(f, t):
	return exp(j*2*pi*f*t)*(4*cos(pi*f)-pi*f*sin(pi*f) - 4*exp(-4)*(2*cos(3*pi*f)-pi*f*sin(3*pi*f)))/(3+pi**2*f**2)



x = [i for i in frange(-5, 5, 0.0001)]
y = [s(i) for i in x]	

fig = plt.figure()

plt.xlim(-5, 5)
plt.ylim(-0.125, 1.125)

plt.title(u'Сигнал')
plt.xlabel(u't, [мс]')
plt.ylabel(u's(t), [В]')
plt.grid(True)

plt.plot(x, y, 'black', linewidth=2)

plt.savefig('SignalPlot.png')
#plt.show()

x = [i for i in frange(-15, 15, 0.0001)]
y = [S(i) for i in x]

fft_fig = plt.figure()

plt.xlim(-10, 10)
plt.ylim(-1, 1.5)

plt.title(u'Спектр')
plt.xlabel(u'$f, [кГц]$')
plt.ylabel(u'$S(f), [В \cdot с]$')
plt.grid(True)

plt.plot(x, y)

plt.savefig('SpectrePlot.png')

Re_fig = plt.figure()

plt.xlim(-10, 10)
plt.ylim(-1, 1.5)

plt.title(u'Действительная часть спектра')
plt.xlabel(u'$f, [кГц]$')
plt.ylabel(u'$|S(f)|, [В \cdot с]$')
plt.grid(True)

ReY = real(y)
plt.plot(x, ReY, 'g-', linewidth=2)

plt.savefig('ReSpectrePlot.png')

#Im_fig = plt.figure()

#plt.xlim(-10, 10)
#plt.ylim(-0.5, 0.5)

#plt.title(u'Мнимая часть спектра')
#plt.xlabel(u'$f, [кГц]$')
#plt.ylabel(u'$Im{S(f)}, [В \cdot с]$')
#plt.grid(True)

#ImY = imag(y)
#plt.plot(x, ImY, 'r-')

#plt.savefig('ImSpectrePlot.png')

#plt.legend(('S(f)','|S(f)|', 'imag{S(f)}'))

#arg_fig = plt.figure()

#plt.xlim(-2, 2)
#plt.ylim(-0.5, 3.5)

#plt.title(u'Аргумент спектра')
#plt.xlabel(u'$f, [кГц]$')
#plt.ylabel(u'arg{S(f)}, $[В \cdot с]$')
#plt.grid(True)

#argY = angle(y)
#plt.plot(x, argY, 'black')

#plt.savefig('ArgSpectrePlot.png')

WY = [abs(w)**2 for w in y]

Wfft_fig = plt.figure()

plt.xlim(-5, 5)
plt.ylim(-0.5, 1.125)

plt.title(u'Энергетический спектр сигнала')
plt.xlabel(u'$f, [кГц]$')
plt.ylabel(u'$W(f), [В^{2} \cdot с]$')
plt.grid(True)

plt.plot(x, WY, 'g-')

plt.savefig('WSpectrePlot.png')
plt.show()


x = [i for i in frange(0, 25, 1e-2)]

print(quad(W, 0, inf))
#E = []
#for f in x:
#	E.append(quad(W, 0, f)[0]/quad(W, 0, inf)[0])

#fig = plt.figure()

#plt.xlim(0, 25)
#plt.ylim(0, 1)

#plt.title(u' Зависимость Доли энергии, попадающей в полосу частот Δf, от ширины полосы.')
#plt.xlabel(u'$f, [кГц]$')
#plt.ylabel(u'$Δf, [кГц]$')
#plt.grid(True)

#plt.plot(x, E, 'black')
#plt.plot(x, [0.75 for i in x], 'r--', linewidth=0.5)
#plt.plot(x, [0.9 for i in x], 'r--', linewidth=0.5)
#plt.plot(x, [0.99 for i in x], 'r--', linewidth=0.5)

#plt.legend(('Δf', '0.75', '0.9', '0.99'))

#plt.show()

#plt.savefig('ESpectrePlot.png')

print('Создаю область определения х для графиков восстановленных сигналов')
x = [i for i in frange(-3, 3, 0.001)]

print('Создаю координатную плоскость')
plt.xlim(-2, 2)
plt.ylim(-0.5, 1.5)

plt.title(u'Восстановленные сигналы')
plt.xlabel(u'$s(t), [В]$')
plt.ylabel(u'$t, [мс]$')
plt.grid(True)

print('Вычисляю интеграл для f75')
F75 = []
for t in x:
	F75.append(quad(Rev, -0.77, 0.77, args=(t))[0])

print('Вычисляю интеграл для f90')
F90 = []
for t in x:
	F90.append(quad(Rev, -1.9, 1.9, args=(t))[0])

print('Вычисляю интеграл для f99')
F99 = []
for t in x:
	F99.append(quad(Rev, -19.5, 19.5, args=(t))[0])

print('Наношу точки на координатную плоскость')
plt.plot(x, F75, 'g--',linewidth=0.5)
plt.plot(x, F90, 'b--',linewidth=1)
plt.plot(x, F99, 'r--',linewidth=1.5)

print('Создаю легенду')
plt.legend(('$s(t)$', '$s\_75 (t)$', '$s\_90 (t)$', '$s\_99 (t)$'))

print('Вывожу график на экран')
plt.savefig('RevSignalPlot.png')
plt.show()
