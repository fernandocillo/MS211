clear;
clc;
close;

#definindo os parâmetros
a=0;
b=pi/2;
n=100;

#a função que iremos integrar numericamente
f=@(x) abs(sin(x)*cos(x))^(3/2);

#Regra dos Trapézios Repetida
function TrapRep = TrapRep(f, a, b, n)
    h = (b - a) / n; #passo da integração (delta x da int)
    integral = 0.5*(f(a) + f(b)); #valor médio
    x = a + h;
    for i = 1:n-1
        integral += f(x); #somatório - aproxima a integral
        x += h;
    endfor
    TrapRep= h*integral;
endfunction

aprox = TrapRep(f, a, b, n);

printf("Integral Numérica: %g\n", aprox)

#Integral exata calculada pelo Wolfram Alpha
I = 0.309012;
printf("\nIntegral Exata: %g\n", I)
printf("\nErro relativo: %g\n", abs(I-aprox)/abs(I))
