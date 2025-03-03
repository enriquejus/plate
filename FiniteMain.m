clear classes
clear
clc

x = [-3,-2,-1,0.5,1];

Coef = zeros(4,6);

for k = 1:4  % Orden derivada

    ini = 5-k;
    x1 = x(ini:5);
    Res = CalcCoefDF(x1,k);

    Zini = zeros(1,4-k)
    Res_t = transpose(Res)
    Coef(k,:) = [Zini,Res_t]

end % for  

Coef

 xprim = [-.5,0.5];
 xseg = [-1.5,-.5,0.5,1.5];
 Resprim = CalcCoefDF(xprim,1)
 Resseg = CalcCoefDF(xseg,2)

xter = [-1.5,-.5,.5,1.5];
ResTer = CalcCoefDF(xter,3)

