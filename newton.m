%MÉTODO DE NEWTON (MN)
f=@(x) x.^3 -9.*x+3;
phi=@(x) x-(f(x)/(3.*x.^2 -9));
x=2;
E=1e-13;
cont=0;
while abs(f(x))>E
 x=phi(x);
 cont=cont+1;
end
MN=x
contagemMN=cont