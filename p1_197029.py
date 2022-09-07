#Fernando Teodoro de Cillo
#197029

### IMPORTANDO AS BIBLIOTECAS QUE USAREMOS ###
import scipy.optimize as optimize #métodos numéricos
import numpy as np #funções matemáticas

"""QUESTÃO 1"""
print("\n QUESTÃO 1 \n")

### DEFININDO A FUNÇÃO ###
def func(x):
    return np.sin((x-0.23)**2)-(np.e)**x**2+3*x

### APLICANDO O MÉTODO DA BISSEÇÃO PARA ENCONTRAR AS RAÍZES ###
x1=optimize.bisect(func, 0, 1) #intervalo 1
x2=optimize.bisect(func, 1, 2) #intervalo 2

print("As raízes encontradas pelo método da bisseção foram "+str(x1)+" e "+str(x2)+"\n")

### DEFININDO A DERIVADA ###
def der(x):
  return 2*(-0.23*x)*np.cos((-0.23+x)**2)-2*np.e**(x**2)*x +3
#a derivada da função foi calculada pelo Wolfram Alpha#

### APLICANDO O MÉTODO DE NEWTON-RAPHSON COM X1 E COM X2 ###
t=10**(-15) #tolerância
x1=optimize.newton(func, x1, fprime=lambda x: der(x), tol=t)
x2=optimize.newton(func, x2, fprime=lambda x: der(x), tol=t)

print("As raízes encontradas pelo método de Newton-Raphson foram "+str(x1)+" e "+str(x2)+"\n")

"""QUESTÃO 2"""

### DEFININDO A FUNÇÃO ###
def func(x):
    return (np.e)**(-(x**2))+x

### APLICANDO O MÉTODO DA BISSEÇÃO ###
x1=optimize.bisect(func, -1, -0.5, rtol=10**(-12))

print("\n QUESTÃO 2 \n")
print("A raiz de f(x) é "+str(x1))