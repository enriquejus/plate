function Coef = FiniteDif(x,clave)
% Calcula matriz de coeficientes para todas las derivadas
% Clave, casos en que no se cumple el caso general
% Clave = 0, caso general
% Clave = 1. i = n-2. La derivada primera y la segunda se hacen simétricas.
% Clave = 2. i = n-1. Bound = 6

Coef = zeros(4,6);

for k = 1:4  % Orden derivada
    
    ini = 5-k;
    x1 = x(ini:5);
    Zini = zeros(1,4-k);
    Res = CalcCoefDF(x1,k);
    Res_t = transpose(Res);
    Coef(k,:) = [Zini,Res_t];
end % for  

    if clave == 1
        Coef(1,:) = [0,0,-.5,.5,0,0];
        Coef(2,:) = [0,0,1,1,0,-2];
    elseif clave == 2
        Coef(1,:) = [0,0,-1/3,1.3333,0,-1];
        Coef(2,:) = [0,-0.2,2,3.2,0,-5];
    end % if clave
    
end