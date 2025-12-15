import math
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.ticker import MultipleLocator


# n = int(input('размер решётки:'))
# t0 = int(input('начальная температура:'))
# t = int(input('конечная температура:'))
# diff_t = int(input('величина шага по температуре:'))
# J = int(input('константа обменного взаимодействия:'))
# term = int(input('количество шагов Монте-Карло для термолизации:'))
# count = int(input('количество шагов Монте-Карло для измерений:'))
n = 10
t0 = 7
t = 0.1
diff_t = 50
tt = (t0-t) / diff_t
J = 1
term = 4000
count = 15000

shape = (n, n, n)
s = np.random.choice([-1, 1], size = shape) #массив, моделирующий систему


energy = np.zeros(diff_t + 1)
c = np.zeros(diff_t + 1)
m = np.zeros(diff_t + 1)
chi = np.zeros(diff_t + 1)

t_for_plot = np.zeros(diff_t + 1)
l = 0 # индекс для характеристик

while t0 >= t:
	diff_energy = 0 #изменение энергии при одной температуре
	beta = t0 ** (-1)

	a_dict = {12 : np.e ** (- beta * 12), 8 : np.e ** (- beta * 8), 4 : np.e ** (- beta * 4)}


	for v in range(term): #термолизация
		for i in range(n):
			for j in range(n):
				for k in range(n):

					diff_e = 2 * s[i, j, k] * (s[(i - 1) % n, j, k] + s[(i + 1) % n, j, k] + s[i, (j - 1) % n, k] + s[i, (j + 1) % n, k] + s[i, j, (k - 1) % n] + s[i, j, (k + 1) % n])

					if diff_e <= 0:
						a = 1
					else:
						a = a_dict[diff_e]

					r = random.random()
					if r < a:
						s[i, j, k] *= (-1)



	en = 0 #энергия исходной системы (искодного состояния при данном t)
	en += np.sum(s * np.roll(s, shift = -1, axis = 0))
	en += np.sum(s * np.roll(s, shift = -1, axis = 1))
	en += np.sum(s * np.roll(s, shift = -1, axis = 2))
	en *= - J

	mag = s.sum()

	energy_for_one_mcs = np.zeros(1500)
	m_for_one_mcs = np.zeros(1500)

	diff_m = 0

	for v in range(count): #для измерений

		for i in range(n):
			for j in range(n):
				for k in range(n):

					diff_e = 2 * s[i, j, k] * (s[(i - 1) % n, j, k] + s[(i + 1) % n, j, k] + s[i, (j - 1) % n, k] + s[i, (j + 1) % n, k] + s[i, j, (k - 1) % n] + s[i, j, (k + 1) % n])

					if diff_e <= 0:
						a = 1
					else:
						a = a_dict[diff_e]

					r = random.random()
					if r < a:
						s[i, j, k] *= (-1)
						diff_energy += diff_e
						diff_m += 2 * s[i, j, k]

		
		if v % 10 == 0:				
						
			energy_for_one_mcs[(v // 10)] = en + diff_energy
			m_for_one_mcs[(v // 10)] = abs(mag + diff_m)
	

	energy[l] = (en + diff_energy) / (n ** 3)
	m[l] = abs((mag + diff_m)) / (n ** 3)
	c[l] = (beta ** 2) * ( (energy_for_one_mcs ** 2).sum() / (1500)  - ((energy_for_one_mcs.sum() / (1500) ) ** 2) ) / (n ** 3) #1.⟨E²⟩ - среднее квадратов энергии ⟨E²⟩ = (E₁² + E₂² + E₃² + ... + Eₙ²) / n,          2. ⟨E⟩² - квадрат средней энергии ⟨E⟩² = [(E₁ + E₂ + E₃ + ... + Eₙ) / n]²
	chi[l] = beta * ( (m_for_one_mcs ** 2).sum() / (1500)  - ((m_for_one_mcs.sum() / (1500) ) ** 2) ) / (n ** 3)

	t_for_plot[l] = t0
	t0 -= tt
	l +=1

#графики для презентации

style = {'size' : 20, 'weight' : 'bold'}

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0,0].plot(t_for_plot, energy, linewidth=2.5)
axes[0,0].set_xlabel(r'$T$', style) #температура
axes[0,0].set_ylabel(r'$E$', style) #удельная энергия
axes[0,0].yaxis.set_label_coords(-0.03, 1.08)
axes[0,0].yaxis.label.set(rotation=0, ha='right', va='top')
axes[0,0].grid(which='major')
axes[0,0].tick_params(axis='both', which='major', labelsize=16)
axes[0,0].xaxis.set_major_locator(MultipleLocator(1))
axes[0,0].yaxis.set_major_locator(MultipleLocator(0.5))

axes[0,1].plot(t_for_plot, m, linewidth=2.5)
axes[0,1].set_xlabel(r'$T$', style)
axes[0,1].set_ylabel(r'$m$', style) #удельная намагниченность
axes[0,1].yaxis.set_label_coords(-0.03, 1.08)
axes[0,1].yaxis.label.set(rotation=0, ha='right', va='top')
axes[0,1].grid(which='major')
axes[0,1].tick_params(axis='both', which='major', labelsize=16)
axes[0,1].xaxis.set_major_locator(MultipleLocator(1))
axes[0,1].yaxis.set_major_locator(MultipleLocator(0.2))

axes[1,0].plot(t_for_plot, c, linewidth=2.5)
axes[1,0].set_xlabel(r'$T$', style)
axes[1,0].set_ylabel(r'$c$', style) #удельная теплоёмкость
axes[1,0].yaxis.set_label_coords(-0.03, 1.08)
axes[1,0].yaxis.label.set(rotation=0, ha='right', va='top')
axes[1,0].grid(which='major')
axes[1,0].tick_params(axis='both', which='major', labelsize=16)
axes[1,0].xaxis.set_major_locator(MultipleLocator(1))
axes[1,0].yaxis.set_major_locator(MultipleLocator(0.4))

axes[1,1].plot(t_for_plot, chi, linewidth=2.5)
axes[1,1].set_xlabel(r'$T$', style)
axes[1,1].set_ylabel(r'$\chi$', style) # удельная магнитная восприимчивость
axes[1,1].yaxis.set_label_coords(-0.03, 1.08)
axes[1,1].yaxis.label.set(rotation=0, ha='right', va='top')
axes[1,1].grid(which='major')
axes[1,1].tick_params(axis='both', which='major', labelsize=16)
axes[1,1].xaxis.set_major_locator(MultipleLocator(1))
axes[1,1].yaxis.set_major_locator(MultipleLocator(1))

plt.tight_layout()
plt.show()




#графики для отчета
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0,0].plot(t_for_plot, energy)
axes[0,0].set_xlabel(r'$T$') #температура
axes[0,0].set_ylabel(r'$E$') #удельная энергия
axes[0,0].grid(which='major')

axes[0,1].plot(t_for_plot, m)
axes[0,1].set_xlabel(r'$T$')
axes[0,1].set_ylabel(r'$m$') #удельная намагниченность
axes[0,1].grid(which='major')

axes[1,0].plot(t_for_plot, c)
axes[1,0].set_xlabel(r'$T$')
axes[1,0].set_ylabel(r'$c$') #удельная теплоёмкость
axes[1,0].grid(which='major')

axes[1,1].plot(t_for_plot, chi)
axes[1,1].set_xlabel(r'$T$')
axes[1,1].set_ylabel(r'$\chi$') # удельная магнитная восприимчивость
axes[1,1].grid(which='major')

plt.tight_layout()
plt.show()