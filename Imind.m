function Res = Imind(r,z,nu,E)
% Factor de influencia para el desplazamiento vertical debido a 
% una carga vertical P (Mindlin)

    R0 = sqrt (r.^2 + z.^2);
    Res = 8*(1-nu)./R0 + 4*z.^2./R0.^3;

    Res = Res.* (1 + nu)./ (8*pi.*E);
    
end % Imind


    
    
    
    