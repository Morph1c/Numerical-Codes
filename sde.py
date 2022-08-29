##
# Simulazione di equazioni differenziali stocastiche
#

# librerie
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from math import exp
from math import log
from sklearn.linear_model import LinearRegression

# variabili globali
mu = float(input("mu: "))
sigma = float(input("sigma: "))

# valuta l'errore forte come definito nella tesina
def strong_error(N, X, Y):
	return log(sum(abs(X - Y)) / N)

def weak_error(N, X, Y):
	return log(abs((sum(X)/ N) - (sum(Y)/N)))

# valutazione derivata che nel caso del nostro modello, ho già calcolato esplicitamente
def derivata(scelta, opt, x):
	if scelta ==1:
		if opt == 1:
			return mu
		else:
			return sigma

# scelta modello che nel nostro caso era solo una ma si può estendere a più modelli
def f(scelta, opt, x):
	if scelta ==1:
		if opt == 1:
			return mu*x
		else:
			return sigma*x
## Metodo di eulero-maruyama
def euler_maruyama(N, dt, dW, X, scelta):
	for j in range(1, N):
		X[j] = X[j- 1] + dt*f(scelta, 1, X[j-1]) + f(scelta, 2, X[j-1])*dW[j]

## Metodo di Milstein
def milstein(N, dt, dW, X, scelta):
	for j in range(1, N):
		X[j] = X[j- 1] + dt*f(scelta, 1, X[j-1]) + f(scelta, 2, X[j-1])*dW[j] + 0.5 * f(scelta, 2, X[j-1]) * derivata(scelta, 2, X[j-1]) * (dW[j]**2 - dt)


def main(mu, sigma):
	T = 1
	plot_data_euler = {}
	plot_data_milstein = {}
	scelta = 1
	c = 1
	N = 10
	m = 2000
	n_samples = 5 # numero di test numerici con differenti passi di discretizzazione
	first_X = np.zeros((N), dtype = np.longdouble)
	first_Y = np.zeros((N), dtype = np.longdouble)
	first_Z = np.zeros((N), dtype = np.longdouble)
	type_error = int(input("1)Errore forte\n2)Errore debole\n: "))
	while(c <= n_samples):
		dt = T / N
		dW = np.zeros((N), dtype = np.longdouble)
		W = np.zeros((N), dtype = np.longdouble)
		X = np.zeros((N), dtype = np.longdouble)
		Y = np.zeros((N), dtype = np.longdouble)
		Z = np.zeros((N), dtype = np.longdouble)
		X_T = np.zeros((m), dtype = np.longdouble)
		Y_T = np.zeros((m), dtype = np.longdouble)
		Z_T = np.zeros((m), dtype = np.longdouble)
		
		X[0] = 1 
		Y[0] = 1
		Z[0] = 1


		for t in range(0, m):
			dW[0] = sqrt(dt) * np.random.normal(0, 1, 1)
			# riempio il vettore degli incrementi gaussiani
			for j in range(1, N):
				dW[j] = sqrt(dt) * np.random.normal(0, 1, 1)
				W[j] = W[j-1] + dW[j]
			# valuto sol approx
			euler_maruyama(N, dt, dW, Y, scelta)
			milstein(N, dt, dW, Z, scelta)
			# valuto soluzione esatta
			for j in range(1, N):
				X[j] = X[0]*exp((mu - 0.5*sigma**2)*(dt*j) + sigma * W[j])
				if c == 1:
					first_X[j] = X[j]
					euler_maruyama(N, dt, dW, first_Y, scelta)
					milstein(N, dt, dW, first_Z, scelta)
			X_T[t] = X[N-1]
			Y_T[t] = Y[N-1]
			Z_T[t] = Z[N-1]
			if(type_error == 1):
				plot_data_euler[log(dt)] = strong_error(m, X_T, Y_T)
				plot_data_milstein[log(dt)] = strong_error(m, X_T, Z_T)
			else:
				plot_data_euler[log(dt)] = weak_error(m, X_T, Y_T)
				plot_data_milstein[log(dt)] = weak_error(m, X_T, Z_T)
	
		c += 1
		N *= 2

	## PLOT PER STUDIO CONVERGENZA FORTE (2° plot nella tesina)
	#
	x_data = np.zeros((n_samples), dtype = np.longdouble)
	y_data_euler = np.zeros((n_samples), dtype = np.longdouble)
	y_data_milstein = np.zeros((n_samples), dtype = np.longdouble)
	y_regression_euler= np.zeros((n_samples), dtype = np.longdouble)
	y_regression_milstein = np.zeros((n_samples), dtype = np.longdouble)

	i = 0
	for k in plot_data_euler.items():
		x_data[i] = k[0]
		y_data_euler[i] = k[1]
		i +=1
	i = 0
	for k in plot_data_milstein.items():
		y_data_milstein[i] = k[1]
		i +=1

	# PLOT PER STUDIO ACCURATEZZA DEI METODI (1° PLOT nella tesina)
	figure, axis = plt.subplots(1, 2)
	R = np.arange(0, T, dt)

	axis[0].plot(R, X, color = "r", label = "exact")
	axis[0].plot(R, Y, color = "g", label = "approx euler")
	axis[0].set_title("Euler-Mauryama with dt = %.3f" %(dt))
	axis[0].legend()
	axis[1].plot(R, X, color = "r", label = "exact")
	axis[1].plot(R, Z, color = "b", label = "approx milstein")
	axis[1].set_title("Milstein with dt = %.3f" %(dt))
	axis[1].legend()
	plt.legend()
	plt.show()
	# PLOT SECONDO 
	x_data = x_data.reshape((-1, 1))

	model_euler = LinearRegression()
	model_euler.fit(x_data, y_data_euler)
	print("slope euler: %.3f" % (model_euler.coef_))
	#y_regression_euler[0] = model_euler.intercept_ + model_euler.coef_ * x_data[0]
	#y_regression_euler[1] = model_euler.intercept_ + model_euler.coef_ * x_data[n_samples - 1]
	y_regression_euler = model_euler.predict(x_data)

	model_milstein = LinearRegression()
	model_milstein.fit(x_data, y_data_milstein)
	print("slope milstein: %.3f" % (model_milstein.coef_))
	#y_regression_milstein[0] = model_euler.intercept_ + model_milstein.coef_ * x_data[0]
	#y_regression_milstein[1] = model_euler.intercept_ + model_milstein.coef_ * x_data[n_samples - 1]
	y_regression_milstein = model_milstein.predict(x_data)
	fig, ax = plt.subplots(1, 2)  
	ax[0].plot(x_data, y_data_euler, color = "g", label = "error euler")
	ax[0].plot(x_data, y_regression_euler, color = "r", label = "regression euler")
	ax[0].set_title("Euler-maruyama slope is %.3f" %(model_euler.coef_))
	ax[0].legend()
	ax[1].plot(x_data, y_data_milstein, color = "b", label = "error milstein")
	ax[1].plot(x_data, y_regression_milstein, color = "y", label = "regression milstein")
	ax[1].set_title("Milstein slope is %.3f" %(model_milstein.coef_))
	plt.legend()
	plt.show()

	
	

if __name__ == "__main__":
	main(mu, sigma)