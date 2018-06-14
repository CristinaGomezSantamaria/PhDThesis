%Plantilla Configuración Escenario
%Esta plantilla configura el escenario de simulación, la geometría de los
%arreglos de antena y los parámetros del canal a simular.
%Para el caso del canal físico genera los dispersores SOLO con fines de
%GRAFICACIÓN de una muestra del Escenario en General!!!
%OJO! Es recomendable poner un número exacto de BS para que la grilla de
%cuadrada, tipo 4, 9, 16, 25, 36... pues si no, funciona pero quedan unos
%MS generados lejos de las BS, por ej si son 7 BS, igual la grilla sería
%3x3 osea con 9 cuadros digamos, de los cuales se usan 7 para las BS, pero
%en el cuadro 8 y 9 igual se distribuyen usuarios... entonces eso queda
%como raro ahí...!!!

function [Escenario]=Plantilla()

%PARAMETROS GENERALES
Escenario.f = 0.9e9;                                                                %Frecuencia de operación COST231Hata 900 MHz, para COST231HataModified puede hasta 2GHz 
Escenario.lamda = 3e8/Escenario.f;                                                  %Longitud de onda según frecuencia, queda en mts
Escenario.num_BS = 1;                                                               %Número estaciones base             
Escenario.num_MS = 1;                                                               %Número de móviles
if Escenario.num_MS==1 && Escenario.num_BS==1
    Escenario.id='SU';
else
    Escenario.id='MU';
end

%GENERACIÓN POSICIONES BS Y MS
Escenario.dref = 1;                                                                               %Distancia de referencia típica 1 mt, según Molisch
Escenario.radio_cell = 3000;                                                                      %Radio de cada celda en mts para GBU-COST259
Escenario.h_BS = 50;                                                                              %Altura típica de las estaciones base en mts para GBU-COST259
Escenario.h_MS = 1.5;                                                                             %Altura típica de los móviles en mts
grilla = ceil(sqrt(Escenario.num_BS));
tam_grilla = grilla*Escenario.radio_cell*2;
[x_grid,y_grid] = meshgrid(Escenario.radio_cell:2*Escenario.radio_cell:grilla*2*Escenario.radio_cell);                            %Generación grilla del escenario
Escenario.BS_posicion = [x_grid(1:Escenario.num_BS); y_grid(1:Escenario.num_BS); ones(1,Escenario.num_BS)*Escenario.h_BS];        %Posición 3D de las BS

if strcmp(Escenario.id,'SU')==1
    Escenario.MS_posicion=[497.6422; 0; Escenario.h_MS];                        %calculado tal que la distancia en mts entre tx y rx de 500m, según GBU-COST259
    Escenario.MS_posicion(1:2,:)=Escenario.MS_posicion(1:2)+Escenario.BS_posicion(1:2,:);
else
    Escenario.MS_posicion=[randint(2,Escenario.num_MS,[0 tam_grilla]); ones(1,Escenario.num_MS)*Escenario.h_MS];
end

%ASIGNACIÓN ESTACIONES BASE A CADA MÓVIL
if Escenario.num_BS>1
    for i=1:Escenario.num_BS
        Escenario.distancia_BS_MS(i,:)=sqrt(sum((repmat(Escenario.BS_posicion(:,i),1,Escenario.num_MS)-Escenario.MS_posicion).^2));
    end
    index_dist_min=find(Escenario.distancia_BS_MS==repmat(min(Escenario.distancia_BS_MS),Escenario.num_BS,1));
    dist_min=Escenario.distancia_BS_MS(index_dist_min);
    dist_min_final=dist_min;
    if length(dist_min)>Escenario.num_MS
        for i=1:length(dist_min)
            [tempfila,tempcolum]=find(dist_min_final==dist_min(i));
            if length(tempfila)>1
                dist_min_final(i)=[];
            end
        end
    end
    [Escenario.BS_asig,Escenario.MS_asig]=find(Escenario.distancia_BS_MS==repmat(dist_min_final',Escenario.num_BS,1));
else
    Escenario.distancia_BS_MS=sqrt(sum((repmat(Escenario.BS_posicion,1,Escenario.num_MS)-Escenario.MS_posicion).^2));
    Escenario.BS_asig=ones(Escenario.num_MS,1);         %índice de BS asignada a cada MS
    Escenario.MS_asig=(1:Escenario.num_MS)';            %índice de cada MS al cual se le asignó una BS
end
Escenario.BS_asignada=Escenario.BS_posicion(:,Escenario.BS_asig);       %posición de BS asignada a cada MS
Escenario.MS_asignado=Escenario.MS_posicion(:,Escenario.MS_asig);       %posición correspondiente de cada MS asignado a cada BS
Escenario.rij=Escenario.MS_asignado-Escenario.BS_asignada;              %vector de BS a sus MS asignados
for bs=1:Escenario.num_BS
    Escenario.deseados{bs}=find(Escenario.BS_asig==bs);
    Escenario.interferentes{bs}=find(Escenario.BS_asig~=bs);
    Escenario.N_d{bs}=length(Escenario.deseados{bs});
    Escenario.N_i{bs}=length(Escenario.interferentes{bs});
end

%CONFIGURACIÓN ARREGLO DE ANTENAS EN BS (ASUMIENDO TODAS LAS BS CON EL MISMO TIPO DE
%ARREGLO DE ANTENAS)
Escenario.d_BS = .25;                                                        %Espaciamiento entre elementos de antena en BS
Escenario.N = 4;                                                            %Número elementos de antena en BS
% BS_array=[zeros(Escenario.N,1) (0*Escenario.d_BS:Escenario.d_BS:(Escenario.N-1)*Escenario.d_BS)' zeros(Escenario.N,1)];  %Posición elementos de antena arreglo en BS, en eje y para resolución azimutal!!!
BS_array=[zeros(Escenario.N,1) ((-Escenario.N+1)/2*Escenario.d_BS:Escenario.d_BS:(Escenario.N-1)/2*Escenario.d_BS)' zeros(Escenario.N,1)];  %Posición elementos de antena arreglo en BS, en eje y para resolución azimutal!!!
Escenario.BS_array=BS_array.';

%CONFIGURACIÓN ARREGLO DE ANTENAS EN MS (ASUMIENDO TODOS LOS MS CON EL MISMO TIPO DE
%ARREGLO DE ANTENAS)
Escenario.d_MS = .25;                                                        %Espaciamiento entre elementos de antena en MS
Escenario.P = 2;                                                            %Número elementos de antena en MS
% MS_array=[zeros(Escenario.P,1) (0*Escenario.d_MS:Escenario.d_MS:(Escenario.P-1)*Escenario.d_MS)' zeros(Escenario.P,1)];    %Posición elementos de antena arreglo en MS, en eje y para resolución azimutal!!!
MS_array=[zeros(Escenario.P,1) ((-Escenario.P+1)/2*Escenario.d_MS:Escenario.d_MS:(Escenario.P-1)/2*Escenario.d_MS)' zeros(Escenario.P,1)];    %Posición elementos de antena arreglo en MS, en eje y para resolución azimutal!!!
Escenario.MS_array=MS_array.';   

%PARÁMETROS GENERALES DEL ESQUEMA DE CAPA FÍSICA
Escenario.tam_trama=10^4;
Escenario.Modulacion='QAM';                                                 %QAM con M=2,4 corresponde exactamente a BSPK y QPSK respectivamente

%PARÁMETROS GENERALES DEL CANAL
Escenario.channel.ID='Physical';                                            %Puede ser 'IIDRayleigh','CorrRayleig' o 'Physical'                                         
switch Escenario.channel.ID
    case 'CorrRayleigh'     %Introducir parámetros de este tipo de canal
        Escenario.A_S=25;                                                   %Angle Spread
    case 'Physical'        %Introducir todos los parámetros de entrada para el simulador físico según A. Molisch
        %Parámetros Globales restantes
        Escenario.channel.Type = 'MacrocellBU';                            %Para otro tipo como microcelda o macrocelda TU, cambian los parámetros y las nubes existentes
        Escenario.K_Rice = 0.1;                                              %Factor de Rice según Molisch
        Escenario.F_Ill = 360;                                               %Función de Iluminación para DS (entre FC y MS), según Molisch
        Escenario.Var_Shadowing_dB = 9;                                     %Desviación estándar de la V.A. lognormal para Shadowing, según COST259
        Escenario.LoS = 'si';
        Escenario.Single_Scatt = 'si';
        Escenario.Double_Scatt = 'si';
        %Tipos de nubes dispersoras
        Escenario.BS_cluster.id='no';                                           %Nube dispersora alrededor de la BS
        Escenario.MS_cluster.id='si';                                           %Nube dispersora alrededor de la MS
        Escenario.FC_cluster.id='si';                                           %Nubes dispersoras lejanas
        
        %Parámetros de cada nube dispersora
        if strcmp(Escenario.MS_cluster.id,'si')
            Escenario.MS_cluster.radio_clusterMS = 50; %[3 2 1 4 3 6 1];       %radio nube alrededor del MS en mts según Molisch, dimensión (1xnum_MS)
            Escenario.MS_cluster.num_scattMS = 50; %[10 3 6 9 2 1 7];           %número de dispersores en nube alrededor del MS según Molisch, dimensión (1xnum_MS)
            Escenario.MS_cluster.MS_Power_dB = -2.2; %-[2 2 3 3 3 4 4];                            %Potencia en dB relativa a PL según Molisch
            Escenario.MS_cluster.MS_Power = 10.^(Escenario.MS_cluster.MS_Power_dB/10);      %Potencia anterior lineal
            Escenario.MS_cluster.pdf_crosssection = 1;                          %Coeficiente de sección transversal, Constante delta(sigma-1) según Molisch
%            Escenario.MS_cluster.pdf_scatt = repmat('gauss_cyl',Escenario.num_MS,1);                 %'uniform  '        %id de pdf de dispersores
            Escenario.MS_cluster.pdf_scatt = repmat('gauss',Escenario.num_MS,1);
            
            posicion_scattMS = zeros(3,sum(Escenario.MS_cluster.num_scattMS));
            prev_MS = 0;
            for ii=1:Escenario.num_MS
                switch Escenario.MS_cluster.pdf_scatt(ii,:)
                    case 'gauss_cyl'
                        Escenario.MS_cluster.parametros{ii} = [1/15 1/(2*Escenario.MS_cluster.radio_clusterMS(ii)^2) -5 10];        %pdf según Molisch es 1/15*exp(-r^2/(2*100^2)) para -5<z<10
                    case 'gauss'
                        Escenario.MS_cluster.parametros{ii} = [Escenario.MS_cluster.radio_clusterMS(ii) -5 10];
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Esta generación de los dispersores se hace SOLO CON FINES DE GRAFICACIÓN!!!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch Escenario.MS_cluster.pdf_scatt(ii,:)
                case 'gauss_cyl'
                    posicion_scattMS(:,prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(ii)) = ...
                        gauss_cyl_AR(Escenario.MS_cluster.parametros{ii},Escenario.MS_posicion(:,ii),Escenario.MS_cluster.num_scattMS(ii));
                case 'gauss'
                    posicion_scattMS(:,prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(ii)) = ...
                        gauss(Escenario.MS_cluster.parametros{ii},Escenario.MS_posicion(:,ii),Escenario.MS_cluster.num_scattMS(ii));
            end
                prev_MS = prev_MS+Escenario.MS_cluster.num_scattMS(ii);
            end
            posicion_scattMS(3,find(posicion_scattMS(3,:)<0)) = 0;            
        end
        
        if strcmp(Escenario.BS_cluster.id,'si')
            Escenario.BS_cluster.radio_clusterBS = 100; %[3 2 6 5];              %radio nubes alrededor de las BS en mts según Molisch, dimensión (1xnum_BS)
            Escenario.BS_cluster.num_scattBS = 4; %[10 5 7 2];                   %número de dispersores en nube alrededor del BS según Molisch, dimensión (1xnum_BS)
            Escenario.BS_cluster.BS_Power_dB = -Inf; %-[5 3 6 1];                            %Potencia en dB relativa a PL según Molisch
            Escenario.BS_cluster.BS_Power = 10.^(Escenario.BS_cluster.BS_Power_dB/10);      %Potencia anterior lineal
            Escenario.BS_cluster.pdf_crosssection = 1;                          %Coeficiente de sección transversal, Constante delta(sigma-1) según Molisch
%             Escenario.BS_cluster.pdf_scatt = repmat('gauss_cyl',Escenario.num_BS,1);   %'uniform  '        %id de pdf de dispersores            
            Escenario.BS_cluster.pdf_scatt = repmat('gauss',Escenario.num_BS,1);   %'uniform  '        %id de pdf de dispersores            

            posicion_scattBS = zeros(3,sum(Escenario.BS_cluster.num_scattBS));            
            prev_BS = 0;
            for ii=1:Escenario.num_BS
                switch Escenario.BS_cluster.pdf_scatt(ii,:)
                    case 'gauss_cyl'
                        Escenario.BS_cluster.parametros{ii} = [1/10 1/(2*Escenario.BS_cluster.radio_clusterBS(ii)^2) -5 5];     %pdf según Molisch es 1/10*exp(-r^2/(2*100^2)) para -5<z<5
                    case 'gauss'
                        Escenario.BS_cluster.parametros{ii} = [Escenario.BS_cluster.radio_clusterBS(ii) -5 5];     %pdf según Molisch es 1/10*exp(-r^2/(2*100^2)) para -5<z<5                        
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Esta generación de los dispersores se hace SOLO CON FINES DE GRAFICACIÓN!!!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 switch Escenario.BS_cluster.pdf_scatt(ii,:)
                    case 'gauss_cyl'
                        posicion_scattBS(:,prev_BS+1:prev_BS+Escenario.BS_cluster.num_scattBS(ii)) = ...
                            gauss_cyl_AR(Escenario.BS_cluster.parametros{ii},Escenario.BS_posicion(:,ii),Escenario.BS_cluster.num_scattBS(ii));                
                     case 'gauss'
                        posicion_scattBS(:,prev_BS+1:prev_BS+Escenario.BS_cluster.num_scattBS(ii)) = ...
                            gauss(Escenario.BS_cluster.parametros{ii},Escenario.BS_posicion(:,ii),Escenario.BS_cluster.num_scattBS(ii));                
                 end                         
                prev_BS = prev_BS+Escenario.BS_cluster.num_scattBS(ii);
            end
            posicion_scattBS(3,find(posicion_scattBS(3,:)<0)) = 0;
        end
        
        if strcmp(Escenario.FC_cluster.id,'si')
            Escenario.FC_cluster.num_FC = 4; %[4 2 3 1];                                             %número de nubes FC presentes en cada celda, esto para el sistema diseñado con OSTBC 3/4, dimensión (1xnum_BS)
            Escenario.FC_cluster.radio_clusterFC = [100 100 100 100]; %[3 1 2 1 4 4 2 9 1 3];       %radio de cada nube FC en mts según Molisch, dimensión (1xsum(num_FC))
            Escenario.FC_cluster.num_scattFC = [10 10 10 10]; %[10 3 2 5 4 3 1 2 5 2];              %número de dispersores en cada nube FC según Molisch, dimension (1xsum(num_FC))
            Escenario.FC_cluster.FC_Power_dB = [-5 -11 -13 -10]; %repmat(-10,1,sum(Escenario.FC_cluster.num_FC));                 %Potencia en dB relativa a PL según Molisch para cada nube FC
            Escenario.FC_cluster.FC_Power = 10.^(Escenario.FC_cluster.FC_Power_dB/10);                                    %Potencia anterior lineal            
            Escenario.FC_cluster.FC_PowerDS_dB = -7;                %Potencia en dB para Double Scattering con MS, relativa a Potencia FC, asumida IGUAL PARA TODOS LOS DS!!!
            Escenario.FC_cluster.FC_PowerDS = 10.^(Escenario.FC_cluster.FC_PowerDS_dB/10);                             %Potencia anterior de DS lineal
            Escenario.FC_cluster.pdf_crosssection = 1;                          %Coeficiente de sección transversal, Constante delta(sigma-1) según Molisch            
%            Escenario.FC_cluster.pdf_scatt = repmat('gauss_cyl',sum(Escenario.FC_cluster.num_FC),1);   %'uniform  '        %id de pdf de dispersores            
            Escenario.FC_cluster.pdf_scatt = repmat('gauss',sum(Escenario.FC_cluster.num_FC),1);   %'uniform  '        %id de pdf de dispersores            

            posicion_scattFC = zeros(3,sum(Escenario.FC_cluster.num_scattFC));
            Escenario.FC_cluster.posicion_FC = zeros(3,sum(Escenario.FC_cluster.num_FC));
            prev_FC = 0;
            for bs=1:Escenario.num_BS
                BS(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_FC(bs)) = kron(Escenario.BS_posicion(:,bs),ones(1,Escenario.FC_cluster.num_FC(bs)));
                if strcmp(Escenario.id,'SU')
                    Escenario.FC_cluster.posicion_FC = [250 -750 0;250 2000 50; 400 -3000 100; 400 1500 50]'+BS(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_FC(bs));
                else
                    Escenario.FC_cluster.posicion_FC(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_FC(bs)) = [randint(2,Escenario.FC_cluster.num_FC(bs),...
                    [-Escenario.radio_cell Escenario.radio_cell]); randint(1,Escenario.FC_cluster.num_FC(bs),[0 Escenario.h_BS*3])]...
                    +BS(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_FC(bs));
                end
                Escenario.FC_cluster.num_scattFC_celda(bs) = sum(Escenario.FC_cluster.num_scattFC(prev_FC+1:prev_FC+Escenario.FC_cluster.num_FC(bs)));
                prev_FC = prev_FC+Escenario.FC_cluster.num_FC(bs);
            end
            prev_FC = 0;
            for ii=1:sum(Escenario.FC_cluster.num_FC)
                switch Escenario.FC_cluster.pdf_scatt(ii,:)
                    case 'gauss_cyl'
                        Escenario.FC_cluster.parametros{ii} = [1 1/(2*Escenario.FC_cluster.radio_clusterFC(ii)^2) 0 0];         %pdf según Molisch es 1*exp(-(r^2-ro)/(2*100^2))
                    case 'gauss'
                        Escenario.FC_cluster.parametros{ii} = [Escenario.FC_cluster.radio_clusterFC(ii) 0 0];         %pdf según Molisch es 1*exp(-(r^2-ro)/(2*100^2))
                end                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Esta generación de los dispersores se hace SOLO CON FINES DE GRAFICACIÓN!!!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                switch Escenario.FC_cluster.pdf_scatt(ii,:)
                    case 'gauss_cyl'
                    posicion_scattFC(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC(ii)) = ...
                        gauss_cyl_AR(Escenario.FC_cluster.parametros{ii},Escenario.FC_cluster.posicion_FC(:,ii),...
                        Escenario.FC_cluster.num_scattFC(ii));  
                    case 'gauss'
                    posicion_scattFC(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC(ii)) = ...
                        gauss(Escenario.FC_cluster.parametros{ii},Escenario.FC_cluster.posicion_FC(:,ii),...
                        Escenario.FC_cluster.num_scattFC(ii));                      
                end
                Escenario.FC_cluster.BS_FC(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC(ii)) = ...
                    kron(BS(:,ii),ones(1,Escenario.FC_cluster.num_scattFC(ii)));
                Escenario.FC_cluster.FC_Power_scatt(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC(ii)) = ...
                    kron(Escenario.FC_cluster.FC_Power(ii),ones(1,Escenario.FC_cluster.num_scattFC(ii)));
                
                prev_FC = prev_FC+Escenario.FC_cluster.num_scattFC(ii);
            end
        end
end
%GRÁFICAS DEL ESCENARIO
figure;
for i=1:Escenario.num_MS
    temp_grafica=[Escenario.BS_asignada(:,i) Escenario.MS_asignado(:,i)];
    plot3(temp_grafica(1,:),temp_grafica(2,:),temp_grafica(3,:),'-k')
    hold on; grid on;
end
hBS=plot3(Escenario.BS_posicion(1,:),Escenario.BS_posicion(2,:),Escenario.BS_posicion(3,:),'*k');
hMS=plot3(Escenario.MS_posicion(1,:),Escenario.MS_posicion(2,:),Escenario.MS_posicion(3,:),'ok');
xlabel('Distance mts');
ylabel('Distance mts');
zlabel('Distance mts');
for i=1:Escenario.num_BS
    [x,y,z]=sphere(10);
    x=1.05*Escenario.radio_cell*x+Escenario.BS_posicion(1,i);
    y=1.05*Escenario.radio_cell*y+Escenario.BS_posicion(2,i);
    z=zeros(size(y));
    surfl(x,y,z);
    shading interp
    colormap(pink)
end
if strcmp(Escenario.channel.ID,'Physical')
    hBScluster=[]; hMScluster=[]; hFCcluster=[];
    if strcmp(Escenario.BS_cluster.id,'si')
        hBScluster=plot3(posicion_scattBS(1,:),posicion_scattBS(2,:),posicion_scattBS(3,:),'+k');
    end
    if strcmp(Escenario.MS_cluster.id,'si')
        hMScluster=plot3(posicion_scattMS(1,:),posicion_scattMS(2,:),posicion_scattMS(3,:),'^k');
    end
    if strcmp(Escenario.FC_cluster.id,'si')
        hFCcluster=plot3(posicion_scattFC(1,:),posicion_scattFC(2,:),posicion_scattFC(3,:),'sk');
    end    
    legend([hBS hMS hBScluster hMScluster hFCcluster],'BS','MS','BScluster','MScluster','FCcluster');
else
    legend([hBS hMS],'BS','MS');
end
end
