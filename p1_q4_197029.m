%Fernando Teodoro de Cillo
%197029
%código feito em conjunto com meu irmão Guilherme

clc
clear
close

%
n=251;
%
% Matriz
m=sparse(n);
for i=1:n
    m(i,i) = 2.35;
endfor
for i=1:n-1
    m(i,i+1) = -0.78;
    m(i+1,i) = -0.41;
endfor
for i=1:n-25
    m(i,i+25) = -0.51;
m(i+25,i) = -0.28;
endfor
%
% Termo independente
for i=1:2:n
b(i) = 1.5;
endfor
for i=2:2:n-1
b(i)=0.75;
endfor
b=b';
A=m;

x0=b\A; 
x0t=x0'; x0gauss=x0'; x0jac=x0';
tau=10^-4; 
kmax=100;

#-------- critério das linhas e de sassenfeld
n=length(A);
alpha=zeros(n,1); bet=zeros(n,1);
for i=1:n
    for j=1:n
        if i<j
            alpha(i)+=(1/abs(A(i,i)))*abs(A(i,j));
            bet(i)+=(1/abs(A(i,i)))*abs(A(i,j))*bet(j);
        elseif i>j
            alpha(i)+=(1/abs(A(i,i)))*abs(A(i,j));
            bet(i)+=(1/abs(A(i,i)))*abs(A(i,j));
            if alpha(i)>=1
                printf("Não converge pelo critério das linhas");
            end
            if bet(i)>=1
                printf("Não converge pelo critério de Sassenfeld");
            end
        end
    end
end

tic; D = diag(diag(A));
M = A-D;
k = 1; Drjac(k) = tau+1;
while (k<=kmax)&&(Drjac(k)>tau)
    k = k+1;
    x = D\(b-M*x0jac);
    Drjac(k) = norm(x-x0jac,inf)/norm(x,inf); 
    x0jac = x;
end; tempo_jacobi=toc;
erro_jacobi=Drjac(k);
residuo_jacobi=norm(b-A*x0jac)/norm(b);


L = tril(A);
U = A-L;
k = 1; Drgauss(k) = tau+1;
while (k<=kmax)&&(Drgauss(k)>tau)
    k = k+1;
    x = L\(b-U*x0gauss);
    Drgauss(k) = norm(x-x0gauss,inf)/norm(x,inf);
    x0gauss = x;
end; tempo_gauss=toc;
erro_gauss= Drgauss(k);
residuo_gauss=norm(b-A*x0gauss)/(norm(b));

tic; x_barra=A\b; tempo_barra=toc;
residuo_barra=norm(b-A*x_barra)/norm(b);