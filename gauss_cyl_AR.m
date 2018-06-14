function [posicion_scatt]=gauss_cyl_AR(parametros,Ro,Num_Disp)

    %INICIALIZACION DE VARIABLES	
    count=0;
	r_vec=zeros(1,Num_Disp);
    lim_h=parametros(3:4);
    k1=parametros(1);
    k2=parametros(2);
       
    %GENERACION POSICION DISPERSORES SEGUN PARAMETROS DE ENTRADA
	%Metodo de Aceptacion/Rechazo para x:  Generacion de x~U[0,1/k2] y de amplitud k1
%     while count < sum(Num_Disp)
%         x=rand;                                                                                            
%         u=rand;
%         f_x=k1*exp(-k2*x^2);
%         M_g_x=k1;
%         if u <= f_x/M_g_x
%             r=x;
%             r_vec(count+1) = r;
%             count=count+1;
%         end
% 	end
while count < sum(Num_Disp)
    x = randint(1,1,[0 1000]);                                                                                            
    u = k1*rand;
    f_x = k1*exp(-k2*x^2);
    if u <= f_x
        r=x;
        r_vec(count+1) = r;
        count=count+1;
    end
end
    
    %Distribucion Uniforme para z y teta
	z=randint(1,Num_Disp,[lim_h(1) lim_h(2)]);
	teta=randint(1,Num_Disp,[0,360]);
	teta=teta*pi/180;
	
    %Conversion de Coordenadas Cilindricas a Cartesianas
	[xvec,yvec]=pol2cart(teta,r_vec);
	xvec=xvec+Ro(1);
	yvec=yvec+Ro(2);
	zvec=z+Ro(3);
    posicion_scatt = [xvec;yvec;zvec];
