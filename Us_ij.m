function dRes = Us_ij(t,r1,ri,z,nu,E)
% Ip 

%format shortG
r = sqrt (vpa(ri).^2 + vpa(r1).^2 - 2*vpa(r1).*vpa(ri).*cos(vpa(t)));

Res = Umind(r,z,nu,E);
Res = Res.*r1;
dRes = double(Res);

end % Ss_ij