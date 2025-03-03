function Res = Umind(r,z,nu,E)
% Factor de influencia para la tensión vertical debido a 
% una carga vertical P (Mindlin)

    R0 = sqrt (r.^2 + z.^2);
    Res = r.*(1+nu).*(4.*z./(R0.^3) - 4*(1-2*nu)./(R0.*(R0+z)))./(8*pi.*E);
    
end % Umind


    
    
    
    