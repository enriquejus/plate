function dRes = Ss_ij(t,r1,ri,z,nu,E)
% Ip 

%format shortG
r = sqrt (vpa(ri).^2 + vpa(r1).^2 - 2*vpa(r1).*vpa(ri).*cos(vpa(t)));

% rd = double(r);
% escero = find (r < 1e-7);
% if (not(isempty(escero)))
%     t
%     r1
%     ri
%     r
% end % if

Res = Smind(r,z,nu,E);
Res = Res.*r1;
dRes = double(Res);

end % Ss_ij