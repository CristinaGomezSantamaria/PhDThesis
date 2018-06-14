function [posicion_scatt]=gauss(parametros,Ro,Num_Disp)

radio_cluster = parametros(1);
h1 = parametros(2);
h2 = parametros(3);

phi = randint(1,Num_Disp,[0 360]);             %V.A.Uniforme entre 0 y 360
r = radio_cluster*randn(1,Num_Disp);        %V.A.Gaussiana con desviación estándar radio_cluster
z = randint(1,Num_Disp,[h1 h2]);               %V.A.Uniforme entre h1 y h2

[x,y] = pol2cart(phi*pi/180,r);
posicion_scatt = [x;y;z]+kron(Ro,ones(1,Num_Disp));

