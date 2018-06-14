function [posicion_scatt]=Uniforme(parametros,Ro,Num_Disp)

    %INICIALIZACION DE VARIABLES	
    count = 0;
    radio = parametros;
    pos_vec = randint(3,Num_Disp,[0 radio]);
    
    posicion_scatt = pos_vec+kron(Ro,ones(1,Num_Disp));
