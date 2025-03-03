classdef cSoilLayer < handle
    %cSoilLayer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private , GetAccess = private)
        
    Nombre % Nombre de la capa 
    Tipo % Tipo de terreno (granular,cohesivo)
    nLayer % numero de la capa
    Zi % Profundidad inicial
    Zf % Profundidad final
    UnitWeight % Peso específico
    Cu % Cohesión sin drenaje
    FrictionAngle % Angulo de rozamiento interno
    E % Modulo de Elasticidad 
    Nu % Coeficiente de Poisson
    d % Diam pilote en este estrato
    Resist % Resistencia del terreno
        
    end
    
    methods
        
        function this = cSoilLayer (IdSoil,nLayer,Zi,SoilData);
            
        this.Nombre = IdSoil(1);
        this.Tipo = IdSoil(2);
        this.nLayer = nLayer;
        this.Zi = Zi;
        this.Zf = SoilData(1);
        this.UnitWeight = SoilData(2);
        this.Cu = SoilData(3);
        this.FrictionAngle = SoilData(4);
        this.E = SoilData(5);
        this.Nu = SoilData(6);
        this.Resist = SoilData(7);
        
        end % Constructor
        
        %function this = cSoilLayer % constructor sin argumentos
        
        %end %Constructor
        
        function disp_ (this)
        
        IdSoil = [this.Nombre,this.Tipo];
        
        disp('IdSoil')
        disp(IdSoil)
        
        DataSoil = [this.Zi,this.Zf,this.UnitWeight, ...
        this.Cu,this.FrictionAngle,this.E,this.Nu];
               
        disp('DataSoil')
        disp(DataSoil)      
        
        end
        
        function draw (this)
        
            y1 = this.Zi;
            y2 = this.Zf;
        
            line ([-2;2],[-y1;-y1],'color','b'); %Dibuja principio capa
            line ([-2;2],[-y2;-y2],'color','b'); %Dibuja fin capa
            
        end % draw
                
        function Sf = Set_Sv (this,Si)
        
            e = this.Zf-this.Zi;
            Sf = Si + this.UnitWeight*e;
            this.Si = Si;
            this.Sf = Sf;
            
            Elements = this.Elements;
            ispart = not(isempty(Elements));
            
        end % set Sv
        
        function Soiltype = get_Soiltype (this)
            
            Soiltype = this.Tipo;
            
        end % get_Soiltype 
        

        function FrictionAngle = get_FrictionAngle (this)
            
            FrictionAngle = this.FrictionAngle;
            
        end % get_FrictionAngle
 

        function Z = get_Z (this)
        
            Z = [this.Zi,this.Zf];
        
        end % get_z
        

        function Nombre = get_Nombre (this)
        
            Nombre = this.Nombre;
        
        end % get_Nombre

        
        function Cu = get_Cu (this)
            
            Cu = this.Cu;
            
        end % get_Cu

        function E = get_E (this)
        
            E = this.E;
            
        end % get_E
        
        function Nu = get_Nu (this)
        
            Nu = this.Nu;
            
        end % get_Nu      
                
        function nLayer = get_nLayer (this)
        
            nLayer = this.nLayer;
            
        end % get_nLayer      

        function Resist = get_Resist(this)
            
            Resist = this.Resist;
            
        end % get_Resist

    
    end % Methods
    
end
