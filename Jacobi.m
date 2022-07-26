clc;
clear;

%%%%%%%%%   METODO JACOBI   %%%%%%%%%%%%
%%%%%%%%%   MATRIZ ECUACIONES   %%%%%%%%
ecuaciones = [4 -1 1; 4 -8 1; -2 1 5];    % 3x3

%%%%%%%%%   RESULTADOS INDEPENDIENTES DE MATRIZ ECUACIONES   %%%%%%%%
resultados = [7; -21; 15];               % 3x1

%%%%%%%%%   VALORES INICIALES   %%%%%%%%
iniciales = [1; 2; 2];

%%%%%%%%%   TOLERANCIA   %%%%%%%%
tolerancia = 0.000001;

%%%%%%%%%   TAMAÑO DE LA MATRIZ   %%%%%%%%
[filasEcuaciones, columnasEcuaciones] = size(ecuaciones);
[filasResultados, columnasResultados] = size(resultados);
[filasIniciales, columnasIniciales] = size(iniciales);

%%%%%%%%%   SE VERIFICA QUE SEA DEL MISMO TAMAÑO LA MATRIZ   %%%%%%%%
    if (filasEcuaciones ~= columnasEcuaciones)
        error('No es una matriz del mismo tamaño');
    else if(filasResultados ~= filasEcuaciones)
            error('Los resultados no son del mismo tamaño que la matriz');
        else if(filasIniciales ~= filasEcuaciones)
                error('Los datos iniciales no son del mismo tamaño que la matriz');
            end
        end
    end

%%%%%%%%%   SE VERIFICA QUE SEA DIAGONALMENTE DOMINANTE   %%%%%%%%

%  for i = 1: filasEcuaciones
%      extraer = abs( ecuaciones(i, :) ); % los dos puntos extraen la columna completa
%      sumatoriaExtraer = sum(extraer) - extraer(i); %sumatoria de los diagonales menos la diagonal
%      
%        if( extraer(i) < sumatoriaExtraer )
%            error('No es una matriz diagonalmente dominante');
%        end    
%  end

for i = 1: filasEcuaciones
    
sumatoriaExtraer = 0;
    for j = 1: columnasEcuaciones
        sumatoriaExtraer = abs( ecuaciones(i, j) ) + sumatoriaExtraer;
    end
    
    sumatoriaExtraer = sumatoriaExtraer - abs( ecuaciones(i, i) );
    
    if( abs( ecuaciones(i, i) ) < sumatoriaExtraer )
        error('No es una matriz diagonalmente dominante');
    end  
end

%%%%%%%%%   SE CALCULA LA DIAGONAL DE LAS ECUACIONES   %%%%%%%%
 
diagonal = diag(diag(ecuaciones) );    % Obtengo la diagonal de las ecuaciones

%%%%%%%%%   SE CALCULA EL L   %%%%%%%%

L = tril(ecuaciones, -1);

%%%%%%%%%   SE CALCULA EL U   %%%%%%%%

U = triu(ecuaciones, 1);

%%%%%%%%%   DECLARAN LAS VARIABLES DE ERROR, DE ITERACIONES,    %%%%%%%% 
iteraciones = 0;
error = inf;

%%%%%%%%%   SE UTILIZA EL MÉTODO DE JACOBI   %%%%%%%%
while error > tolerancia
 iteraciones = iteraciones + 1;
 
 jacobi = ( ( diagonal )\ resultados ) - ( ( diagonal )\( ( U + L )*iniciales ) );
 
 error = max( abs( (jacobi - iniciales) ) ); %se necesita el cambio maximo de los parametros junto con el valor absoluto
 iniciales = jacobi;
end