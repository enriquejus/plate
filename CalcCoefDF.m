function Coef = CalcCoefDF(x,k)
% Calcula coeficientes diferencias finitas de la derivada k, con los puntos
% definidos en el vector x (la distancia real es multiplicando las x por el 
% intervalo tipo d
% x es un vector fila de dimension k+1
% c es la solucion del sistema, luego se le aÃ±ade un 1 al final, son los
% coef por los que se multiplica cada ecuacion.
% El último coeficiente de cada derivada es el correspondiente al punto de
% cálculo (a)
% Res1 = coeficientes con decimales
% Res2, denom = coeficientes factorizados

A = zeros (k);
b = zeros (k,1);

for i = 1:k;
    i1 = i;
    if i == k
        i1 = k+1;
    end % if
    b(i,1) = -x(1)^i1/factorial(i1);
    for j = 1:k;
        A(i,j) = x(j+1)^i1/factorial(i1);
    end % j
end % i

c = A\b;
c = [1;c]

k1 = sum(c);
k2 = (x.^k)*c;

Res2 = c;
Res2 (k+2) = -k1;
denom = k2/factorial(k);
Coef = Res2 / denom;
end