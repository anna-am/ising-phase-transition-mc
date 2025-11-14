import math
import numpy as np
import matplotlib.pyplot as plt
import random
#from scipy.constants import Boltzmann


#n = int(input('размер решётки:'))
# t0 = int(input('начальная температура:'))
# t = int(input('конечная температура:'))
# diff_t = int(input('количество шагов Монте-Карло:'))
n = 4
t0 = 7
t = 2
diff_t = 15000
tt = (t0-t) / diff_t


shape = (n, n, n)
s = np.random.choice([-1, 1], size = shape) #массив, моделирующий систему

#print(s)

en = 0 #энергия исходной системы
en += np.sum(s * np.roll(s, shift = -1, axis = 0))
en += np.sum(s * np.roll(s, shift = -1, axis = 1))
en += np.sum(s * np.roll(s, shift = -1, axis = 2))
en *= -1
#print(en)

energy = np.zeros(diff_t + 1)
c = np.zeros(diff_t + 1)
m = np.zeros(diff_t + 1)
xi = np.zeros(diff_t + 1)

t_for_plot = np.zeros(diff_t + 1)
l = 0 # индекс для характеристик

energy[0] = en / n
m[0] = s.sum() / n
c[0] = 0
xi[0] = 0



t_for_plot[0] = t0

while t0 <= t:
	#beta = (Boltzmann * t0) ** (-1)
	diff_energy = 0 #изменение энергии на одном шаге монте карло
	beta = (t0 * 1.380649 * (10 ** (-23))) ** (-1)

	a_dict = {12 : np.e ** (- beta * 12), 8 : np.e ** (- beta * 8), 4 : np.e ** (- beta * 4)}

	for i in range(n):
		for j in range(n):
			for k in range(n):
				s[i, j, k] *= (-1)

				diff_e = 2 * s[i, j, k] * (s[(i - 1) % n, j, k] + s[(i + 1) % n, j, k] + s[i, (j - 1) % n, k] + s[i, (j + 1) % n, k] + s[i, j, (k - 1) % n] + s[i, j, (k + 1) % n])

				if diff_e <= 0:
					a = 1
				else:
					a = a_dict[diff_e]

				r = random.random()
				if r >= a:
					s[i, j, k] *= (-1)
				else:
					diff_energy += diff_e

				k += 1
			j += 1
		i += 1
	
	l += 1
	energy[l] = (energy[l - 1] + diff_energy) / (n ** 3)
	m[l] =  s.sum() / (n ** 3)
	c[l] = (beta ** 2) * (((energy.sum() / (l + 1) ) ** 2) - (energy ** 2).sum() / (l + 1)) / (n ** 3) #1. ⟨E⟩² - квадрат средней энергии ⟨E⟩² = [(E₁ + E₂ + E₃ + ... + Eₙ) / n]²,          2. ⟨E²⟩ - среднее квадратов энергии ⟨E²⟩ = (E₁² + E₂² + E₃² + ... + Eₙ²) / n
	xi[l] = beta * (n ** 3)* (((m.sum() / (l + 1) ) ** 2) - (m ** 2).sum() / (l + 1))


	t0 += tt
	t_for_plot[l] = t_for_plot[l-1] + tt

# print(energy)
# print(t_for_plot)


plt.plot(t_for_plot, energy)
plt.xlabel('температура')
plt.ylabel('энергия')
plt.grid(which='major')
plt.show()



plt.plot(t_for_plot, m)
plt.xlabel('температура')
plt.ylabel('намагниченность')
plt.grid(which='major')
plt.show()

plt.plot(t_for_plot, c)
plt.xlabel('температура')
plt.ylabel('удельная теплоёмкость')
plt.grid(which='major')
plt.show()

plt.plot(t_for_plot, m)
plt.xlabel('температура')
plt.ylabel('удельная магнитная восприимчивасть')
plt.grid(which='major')
plt.show()


