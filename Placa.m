clear classes
clear
clc

%FilesDat = ["Datos-1.xlsx" "Datos-2.xlsx" "Datos-3.xlsx"];
%FilesRes = ["Resultados.xlsx" "Resultados-2.xlsx" "Resultados-3.xlsx"];

%FilesDat = ["Datos-1.xlsx"];
DataFile = 'Datos-1.xlsx';

FilesRes = 'Resultados.xlsx';

%for i=1:size(FilesDat,2)
    
%DataFile = FilesDat(i)
soil = cSoil (DataFile,FilesRes);

Q0 = linspace(10,10,6);

Solve_pl_inc (soil,60) 

%end % for i
