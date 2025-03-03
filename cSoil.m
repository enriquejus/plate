classdef cSoil < handle
    % cSoil class
    
    properties (SetAccess = private , GetAccess = private)
    
    nLayers    % nº de capas 
    ResFile; % Archivo resultados
    SoilLayers % Capas de terreno
    soildata % Parámetros del terreno (para escritura)
%    writeout = 1; % ¿Calcula y escribe los asientos fuera de la placa?
    Calc_Tens = 0; % Calcula tensiones verticales?
    Calc_AsLayer = 0; % Calcula asientos en inicio capa?
    Calc_SinPlaca = 0; % Calcula los asientos del terreno con la carga sin placa
    Calc_DespHor = 1 ; % Calcula desplazamientos horizontales en superficie
    q = 0;     % Carga superficial del líquido (kPa)
    nelem = 0; % nº elementos placa
    Rp % Radio de la placa circular
    Ep % Módulo de elasticidad de la placa
    h % Canto de la placa
    Nu_p % Coef Poisson de la placa
    PlacData % Vector de datos de la placa
    nPointsS % nº ptos cálculo de asiento incluyendo dentro y fuera placa
    howfarS  % nº radios para cálculo de asiento fuera placa
    d % Ancho del elemento
    Elem_e;  % vector de elementos que aún no han plastificado
    
    end %Properties

    properties 
        
        ELEMENTCHANGE 
        
    end % Public Properties
    
    methods
        
    % Constructor
            
        function this = cSoil(DataFile,ResFile)
        
        % Lee datos de la placa
        this.ResFile = ResFile;
        [PlacData] = this.readPlacData (DataFile);
        this.Rp = PlacData(1);
        this.h = PlacData(2);
        this.Ep = PlacData(3);
        this.Nu_p = PlacData(4);
        this.nelem = PlacData(7);
        this.q = PlacData(8);
        this.nPointsS = PlacData(9);
        this.howfarS = PlacData(10);
        this.d = (this.Rp)/(this.nelem); % Ancho del elem.
        this.Elem_e = linspace(1,this.nelem,this.nelem); % Todos elásticos inicialmente

        % Lee datos del terreno

        [IdSoil,SoilData] = this.readSoilData (DataFile);
        
        this.soildata = SoilData;
    
        % Crea SoilLayers
        
        zi = 0;  % Init zi
        nLayers = size(IdSoil,1);
        this.nLayers = nLayers;
        
            for i = 1:nLayers;
            
            SoilLayers(i) = cSoilLayer(IdSoil(i,:),i,zi,SoilData(i,:));
                    
            zi = SoilData(i,1);
            
            end  % for i
            
        this.SoilLayers = SoilLayers;
        this.Resist = get_Resist(SoilLayers(i)); % Resistencia del terreno

        end % Constructor
 

        function [idsoil,soildata] = readSoilData (this,file_name) 
       
            [soildata,idsoil] = xlsread (file_name,'Terreno','A2:I20');
               
        end % readSoilData
        
        
        function [PlacData] = readPlacData (this,file_name) 
       
            [PlacData] = xlsread (file_name,'Placa', ...
            'A2:J2');
            
        end % readPlacData

        function disp_ (this)
            
            for i = 1:this.nLayers;
                        
            disp_(this.SoilLayers(i));
        
            end % for
            
        end % disp      
   
        function Isij = Calc_Is (this,rini,rfin,ri,zi,Nu,E)

            f = @(t,r1)Is_ij(t,r1,ri,zi,Nu,E);
            Isij = integral2(f,0,2*pi,rini,rfin,'RelTol',1e-3);
                        
        end % Calc_Us      

        function Usij = Calc_Us (this,rini,rfin,ri,zi,Nu,E)

            f = @(t,r1)Us_ij(t,r1,ri,zi,Nu,E);
            Usij = integral2(f,0,2*pi,rini,rfin,'RelTol',1e-3);
                        
        end % Calc_Us    

        function [Is,Is_Capa,Ss,Us] = AsientosTerreno (this) %-----------------------------
        
        % Cálculo matriz asientos terreno
        options = optimset('MaxFunEvals',50000);    
        SoilLayers = this.SoilLayers;  
        n = this.nelem;
        Is(n,n)=0; % Preallocate memory for Is
        Is_Capa(n,n,2*(this.nLayers)) = 0;
        Ss(n,n,this.nLayers)=0; % Preallocate memory for Ss 
        Us(n,n) = 0; % Preallocate memory for Us
        Rp = this.Rp;
        h = this.h;        
        d = this.d; % Ancho del elem.
        
        for j=1:n % ojo n

                % Datos de la placa en elem. i
                
                rini = d*(j-1);
                rfin = d*j;
                
            for i=1:n % ojo n
                
                ri = i*d - .5*d;
                
                Is(i,j) = 0; %inicializar Isij
                                                    
                % Para todas las capas de terreno
                    
                for k = 1:this.nLayers;
                    
                    Layer_k = SoilLayers(k);                    
                    Zk = get_Z (Layer_k); 
                    zi_layerk = Zk(1); % Prof inicial y final de capa
                    zf_layerk = Zk(2);
                    Ek = get_E (Layer_k);
                    Nuk = get_Nu (Layer_k);            
                    Ishi = Calc_Is(this,rini,rfin,ri,zi_layerk,Nuk,Ek);
                    Ishf = Calc_Is(this,rini,rfin,ri,zf_layerk,Nuk,Ek);
            
                    Is(i,j) = Is(i,j) + Ishi-Ishf;  

           % Asientos en inicio y final de capas                    
                    if(this.Calc_AsLayer == 1)
                      Is_Capa(i,j,2*k-1) = Ishi;
                      Is_Capa(i,j,2*k) = Ishf;
                    end % if Calc_AsLayer
                    
            % Tensión vertical a principio de capa
                    if(this.Calc_Tens == 1)
                      Ss(i,j,k) = Calc_Ss(this,rini,rfin,ri,zi_layerk,Nuk,Ek); 
                    end % if Calc_Tens

            % Tensión vertical a principio de capa
                    if((this.Calc_DespHor == 1) && (k == 1))
                      Us(i,j) = Calc_Us(this,rini,rfin,ri,0,Nuk,Ek); 
                    end % if Calc_DespHor
                    
                end % for k
            end % for i                        
            
        end % for j            

        end % AsientosTerreno

        function [Is_out,Ss_out] = AsientosTerrenoOut (this) %-----------------------------
        
        % Cálculo matriz asientos terreno
        options = optimset('MaxFunEvals',50000);    
%        Elements = this.get_Elements;
        SoilLayers = this.SoilLayers;        
        n = this.nelem;
        nPointsS = this.nPointsS;
        Is_out(nPointsS,n)=0; % Preallocate memory for Is
        Ss_out(nPointsS,n,this.nLayers)=0; %...and for Ss

        Ep = this.Ep;
        Nu_p = this.Nu_p;
        Rp = this.Rp;
        h = this.h;        
        D = Ep*h^3/(12*(1-Nu_p^2)); % Rigidez a flexion
        d = this.d; % Ancho del elem.
        pointsS = linspace(Rp,this.howfarS*Rp,nPointsS+1);
        pointsS(1) = []; % Eliminar Rp de los ptos de calc.
        
        for j=1:n                 
                
                rini = d*(j-1);
                rfin = d*j;
                
            for i=1:nPointsS
                
                ri = pointsS(i);
                Is_out(i,j) = 0; %inicializar Is_outij
                                                    
                % Para todas las capas de terreno
                                    
                for k = 1:this.nLayers;
                    
                    Layer_k = SoilLayers(k);                    
                    Zk = get_Z (Layer_k); 
                    zi_layerk = Zk(1); % Prof inicial y final de capa
                    zf_layerk = Zk(2);
                    Ek = get_E (Layer_k);
                    Nuk = get_Nu (Layer_k);            
                    
                    Ishi = Calc_Is(this,rini,rfin,ri,zi_layerk,Nuk,Ek);
                    Ishf = Calc_Is(this,rini,rfin,ri,zf_layerk,Nuk,Ek);
                    
                    Is_out(i,j)=Is_out(i,j) + Ishi-Ishf;                    

            % Tensión vertical a principio de capa
                    if(this.Calc_Tens == 1)
                      Ss_out(i,j,k) = Calc_Ss(this,rini,rfin,ri,zi_layerk,Nuk,Ek); 
                    end % if Calc_Tens

                end % for k
                                
            end % for i                        
            
        end % for j            

        end % AsientosTerrenoOut
    

        function [Ip,Y] = AsientosPlaca (this) %----------------------

        % Crea la matriz de asientos de la placa
        
        n = this.nelem;
        Ip = zeros (n); 
        Ep = this.Ep;
        Nu = this.Nu_p;
        Rp = this.Rp;
        h = this.h;        
        d = this.d; % Ancho del elem.
        D = Ep*h^3/(12*(1-Nu^2)); % Rigidez a flexion

            for i = 1:n % ojo n
                
                % Datos de la placa en elem. i
                
                ri = d*(i-1);
                rf = d*i;
                r = .5*(ri+rf);
                                    
                if(i == 1) % Primer elemento
                    
                    Ip_2 = r/d^4 - 1/d^3;
                    Ip_1 = -4*r/d^4 + 2/d^3 - 1/(r*d^2) - 1/(2*r^2*d);
                    Ip(i,i) = (6*r/d^4 + 2/(r*d^2) + Ip_1)/r;
                    Ip(i,i+1) = (-4*r/d^4 - 2/d^3 - 1/(r*d^2) + 1/(2*r^2*d)...
                    + Ip_2)/r;
                    Ip(i,i+2) = (r/d^4 + 1/d^3)/r;

                elseif (i == 2) % Segundo elemento

                    Ip_2 = r/d^4 - 1/d^3;
                    Ip(i,i-1) = (-4*r/d^4 + 2/d^3 - 1/(r*d^2) - 1/(2*r^2*d)...
                    +Ip_2)/r;
                    Ip(i,i) = (6*r/d^4 + 2/(r*d^2))/r;
                    Ip(i,i+1) = (-4*r/d^4 - 2/d^3 - 1/(r*d^2) + 1/(2*r^2*d))/r;
                    Ip(i,i+2) = (r/d^4 + 1/d^3)/r;                                      
                
                elseif (i == n-1) % Penultimo elemento

                    A = 1 + 2*Nu*d/Rp;
                    B = 1 - 2*Nu*d/Rp;
                    C = (1 + Nu)*d^2/Rp^2 + 3;
                    Ipi2 = 1/d^4 + 1/(r*d^3);
                    Ip(i,n-3) = 1/d^4 - 1/(r*d^3);
                    Ip(i,n-2) = -4/d^4 + 2/(r*d^3) - 1/(r^2*d^2) - 1/(2*r^3*d);
                    Ip(i,n-1) = 6/d^4 + 2/(r^2*d^2) + Ipi2*(-2/(C-B));
                    Ip(i,n) = -4/d^4 - 2/(r*d^3) - 1/(r^2*d^2) + 1/(2*r^3*d)...
                    + Ipi2*(C+A)/(C-B);

                elseif (i == n) % Ultimo elemento
               
                    Ipi1 = -4/d^4 - 2/(r*d^3) - 1/(r^2*d^2) + 1/(2*r^3*d);
                    Ipi2 = 1/d^4 + 1/(r*d^3);                    
                    Ip(i,n-2) = 1/d^4 - 1/(r*d^3);
                    Ip(i,n-1) = -4/d^4 + 2/(r*d^3) - 1/(r^2*d^2) - 1/(2*r^3*d)...
                    + Ipi1*(-2/(C-B)) + Ipi2*(-1-2*B/(C-B));
                    Ip(i,n) = 6/d^4 + 2/(r^2*d^2) + Ipi1*(C+A)/(C-B) + ...
                    Ipi2*(A+B*(C+A)/(C-B));

                else % Elemento generico
                                                           
                    Ip(i,i-2) = 1/d^4 - 1/(r*d^3);
                    Ip(i,i-1) = -4/d^4 + 2/(r*d^3) - 1/(r^2*d^2) - 1/(2*r^3*d);
                    Ip(i,i) = 6/d^4 + 2/(r^2*d^2);
                    Ip(i,i+1) = -4/d^4 - 2/(r*d^3) - 1/(r^2*d^2) + 1/(2*r^3*d);
                    Ip(i,i+2) = 1/d^4 + 1/(r*d^3);
                    
                end % if
                
            end % for i
            
        Ip = D*Ip;

        end % AsientosPlaca

% --------------------------------------------------------------------

        function P = Solve_pl_inc (this,Q0) 
            % Calculo con plastificacion  
            % P0 = Carga en cabeza pilote;    
            
            % inicializar variables
            % q = carga superficial en la placa
            bound = this.bound;
            ResFile = this.ResFile;
            digits(24);
            q = this.q;            
            n = this.nelem;
            Y = zeros (n,1);
            Pf = zeros(n,1);
            q = this.q;
            Rp = this.Rp;
            Nu = this.Nu_p;
            d = this.d;
            Resist = this.Resist;
            num_it = length(Q0)
            Elem_e = this.Elem_e; % vector con los num. de elem. no plastif.
            Pt = zeros(n,num_it);
            Wt = zeros(n,num_it);
            Wcapa = zeros(n,2*(this.nLayers));
            Scapa = zeros(n,2*(this.nLayers));
            Usup = zeros(n);
            Sout = zeros(this.nPointsS,2*(this.nLayers));
            
            Ip = AsientosPlaca (this); 
            [Is,Is_Capa,Ss,Us] = AsientosTerreno (this);            
            
            Is_inv = inv(Is);
            M = Ip + Is_inv; % Matriz inicial del sistema
            
            for it = 1:num_it;
                
                for i = 1:n
                    Y(i) = Q0(it); % incremento de carga en término independiente
                end % for i
             W = M \ Y; % Calc. desplazamientos vert
             if (it == 1)
                 Wt(:,it) = W;
             else
                 Wt(:,it) = Wt(:,it-1) + W; % Acumular asientos de la iteración             
             end % if
             
             Pel = Is \ W; % Calc. presiones placa-terreno (solo elem no plast)
             Elem_p = []; % Inicializar elementos que plastifican en la iteración
             i_plast = []; % y sus índices
             
              for i = 1:length(Elem_e) % para elem. no plastif, comprobar plastificacion
                  
                if (it == 1)
                    Pt(Elem_e(i),1) = Pel(i);
                else
                    Pt(Elem_e(i),it) = Pt(Elem_e(i),it-1)+ Pel(i); % Actualizar presión
                end % if

                 if(Pt(Elem_e(i),it) >= Resist) 
                    disp('plastifica elemento')
                    i
                    Pt (Elem_e(i),it) = Resist; % Igualar presión a la resist.
                    M(Elem_e(i),:) = Ip(Elem_e(i),:); % cambiar la fila de la matriz M
                    i_plast = [i_plast,i];
                    Elem_p = [Elem_p,i];
                 end % if
              end % for i

            Elem_e(Elem_p) = []; % eliminar los elem. plastificados en la iteración
            Is(i_plast,:) = []; %  quitar fila y col. correspond. de Is
            Is(:,i_plast) = [];
            
            end % for it
             
            Pf = Pt(:,it)
            [Is_out,Ss_out] = AsientosTerrenoOut(this); % Matriz para calc asientos terreno out            
            Wout = Is_out*Pf;

            if(this.Calc_SinPlaca == 1)
              q0 = ones(n,1);
              Pf = q0 * q;
            end % if sin placa

%    Presiones finales con o sin placa            
      % Calc asientos en cada capa
            if(this.Calc_AsLayer == 1)              
              for k = 1:2*(this.nLayers);
                Wcapa(:,k) = Is_Capa(:,:,k)*Pf;
              end % next k              
            end % if Calc_AsLayer

    % Calc tensiones verticales en inicio capas          
            if(this.Calc_Tens == 1)
              for k = 1:this.nLayers;
                Scapa(:,k) = Ss(:,:,k)*Pf;
                Sout(:,k) = Ss_out(:,:,k)*Pf;
              end % next k
            end % if Calc_Tens

    % Calc asientos horizontales en superficie          
            if(this.Calc_DespHor == 1)
                Usup = Us*Pf;
            end % if Calc_DespHor 

      % Calc desplazam. en borde placa (r = Rp)          
            A = 1 + 2*Nu*d/Rp;
            B = 1 - 2*Nu*d/Rp;
            C = (1 + Nu)*d^2/Rp^2 + 3;
            F = .5 + Nu*d/(8*Rp);
            G = .5 - Nu*d/(8*Rp);
            Wb = -2*F*Wt(n-1,it)/(C-B) + Wt(n,it)*(G+F*(C+A)/(C-B));
                                            
            % Escribe resultados en archivo resul.xlsx
                        
            fu_write{this.nelem+1,5} = [];
                        
            for i = 1:n; 

               ri = d*i - d/2;

               fu_write{i,1} = i;
               fu_write{i,2} = ri;
               fu_write{i,3} = ri/Rp;
               fu_write{i,4} = -1000*Wt(i,it);
               fu_write{i,5} = Pt(i,it);
 
            end % for i

            fu_write{n+1,1} = n+1;
            fu_write{n+1,2} = Rp;
            fu_write{n+1,3} = 1;
            fu_write{n+1,4} = -1000*Wb;
            fu_write{n+1,5} = 0;

            Sheet = num2str(now);
            Header = {'Elemento','ri','ri/Rp','W (mm)', 'P (kPa)','Rp (m)','q (kPa)'...
            ,'Nº elem','Ep','Nu_p','h'};
            xlswrite(ResFile,Header,Sheet);
            xlswrite(ResFile,fu_write,Sheet,'A2');
            xlswrite(ResFile,[Rp,q,n,this.Ep,this.Nu_p,this.h],Sheet,'F2');
            PosNum = n+5;
                        
%            if (this.writeout == 1)
                
               PosString = int2str(PosNum);
               Pos = ['A' PosString];             
               xlswrite(ResFile,'Asientos fuera',Sheet,Pos);

               PosNum = PosNum + 1;
               PosString = int2str(PosNum);
               Pos = ['A' PosString];             
               Header = {'Punto','ri','ri/Rp','W (mm)'};
               xlswrite(ResFile,Header,Sheet,Pos); 
                
               fout_write{this.nPointsS,4} = [];    
               pointsS = linspace(Rp,this.howfarS*Rp,this.nPointsS+1);
               pointsS(1) = [];

                for l = 1:this.nPointsS;

                   rl = pointsS(l);
                   fout_write{l,1} = l;
                   fout_write{l,2} = rl;
                   fout_write{l,3} = rl/Rp;
                   fout_write{l,4} = -1000*Wout(l);                
 
                end % for l                  
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];            
            xlswrite(ResFile,fout_write,Sheet,Pos);
            PosNum = PosNum + this.nPointsS + 5;            

%            end % if writeout
            
% Escribe asientos en capas

            if (this.Calc_AsLayer == 1)
                
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,'Asientos en capas',Sheet,Pos);        
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];            
            xlswrite(ResFile,-1000*Wcapa,Sheet,Pos);
            PosNum = PosNum +n+5;
 
            end % if Calc_AsLayer
            
            if (this.Calc_Tens == 1)
                
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,'Tensión Vert en capas',Sheet,Pos);        
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];            
            xlswrite(ResFile,Scapa,Sheet,Pos);               
            PosNum = PosNum + n+5;

            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,'Tensiones en capas fuera',Sheet,Pos);        
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];            
            xlswrite(ResFile,Sout,Sheet,Pos);
            PosNum = PosNum + this.nPointsS + 5;                        

            end % if Calc_Tens

            if (this.Calc_DespHor == 1)
                
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,'Desplazamientos horizontales en sup',Sheet,Pos);        
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];            
            xlswrite(ResFile,Usup,Sheet,Pos);               
            PosNum = PosNum + n+5;

            end % if Calc_Desphor

            it1 = [1:it];
            Q1 = 0;
            
            for i = 1:it
                Q1 = Q1 + Q0(i);
                Qtot (i) = Q1;
            end % for i
            Qtot
            
            PosString = int2str(PosNum);
            Pos = ['A' PosString];     
            Header = {'Asientos it'};
            xlswrite(ResFile,Header,Sheet,Pos); 
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];     
            xlswrite(ResFile,Qtot(it1),Sheet,Pos);
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,Wt(:,it1),Sheet,Pos);

            PosNum = PosNum + this.nelem + 4;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];                
            Header = {'Presiones it'};
            xlswrite(ResFile,Header,Sheet,Pos); 
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];     
            xlswrite(ResFile,Qtot(it1),Sheet,Pos);
            
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,Pt(:,it1),Sheet,Pos);
                      
            PosNum = PosNum + this.nelem + 4;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,Ip,Sheet,Pos);

            PosNum = PosNum + this.nelem + 4;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,Y,Sheet,Pos);            

            PosNum = PosNum + this.nelem + 4;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,Is,Sheet,Pos);

            PosNum = PosNum + this.nelem + 4;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,Is_inv,Sheet,Pos);

            PosNum = PosNum + this.nelem + 4;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            xlswrite(ResFile,M,Sheet,Pos);

            PosNum = PosNum + this.nelem + 4;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];             
            Header2 = {'zf','Peso Esp','Cu','Fi','E','Nu', 'Resist'};
            xlswrite(ResFile,Header2,Sheet,Pos);
            PosNum = PosNum + 1;
            PosString = int2str(PosNum);
            Pos = ['A' PosString];
            xlswrite(ResFile,this.soildata,Sheet,Pos);

        end % Solve_pl_inc            

    end % methods
    
end % class

