function Res = Smind(r,z,nu,E)
% Factor de influencia para la tensión vertical debido a 
% una carga vertical P (Mindlin)

    R0 = sqrt (r.^2 + z.^2);
    Res = -3*z.^3./(2*pi.*R0.^5);
    
end % Smind


    
    
    
    