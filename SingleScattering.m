%Esta función calcula los parámetros de los rayos que pasan por un proceso
%de Single Scattering


function [SS]=SingleScattering(pos_scatt,nube,LoS,NLoS)

global Escenario

switch nube    
    case 'MS'
        prev_NR = 0;
        prev_MS = 0;
        for bs=1:Escenario.num_BS
            SS.num_rayos(bs) = sum(Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}));
        end
        rij = zeros(3,sum(SS.num_rayos));
        rjk = zeros(3,sum(SS.num_rayos));
        Tao0 = zeros(1,sum(SS.num_rayos));
        Potencia = zeros(1,sum(SS.num_rayos));
        for bs=1:Escenario.num_BS
            for ms=1:Escenario.N_d{bs}
                MS = pos_scatt(:,prev_NR+1:prev_NR+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
                rij(:,prev_NR+1:prev_NR + Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = MS - ...
                    kron(Escenario.BS_posicion(:,bs),ones(1,Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))));

                rjk(:,prev_NR+1:prev_NR + Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = MS - ...
                    kron(Escenario.MS_posicion(:,Escenario.deseados{bs}(ms)),ones(1,Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))));
                
                Tao0(prev_NR+1:prev_NR + Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                    kron(LoS.Tao(Escenario.deseados{bs}(ms)),ones(1,Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))));

                Potencia(prev_NR+1:prev_NR + Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                    kron(ones(1,Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))),...
                    Escenario.MS_cluster.MS_Power(Escenario.deseados{bs}(ms))*NLoS.P(Escenario.deseados{bs}(ms)));

                prev_MS = prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));
                prev_NR = prev_NR + Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));

            end
        end
                   
    case 'BS'
        %Como en el caso de FC toca primero calcular el número de rayos lo cual permitirá hacer las inicializaciones correspondientes, lo
        %cual justifica esta otra iteración
        for bs=1:Escenario.num_BS
            SS.num_rayos(bs) = Escenario.BS_cluster.num_scattBS(bs)*Escenario.N_d{bs};
        end
        %Preasignación de Memoria        
        rij = zeros(3,sum(SS.num_rayos));
        rjk = zeros(3,sum(SS.num_rayos));
        Tao0 = zeros(1,sum(SS.num_rayos));
        Potencia = zeros(1,sum(SS.num_rayos));
        prev_BS = 0;
        prev_NR = 0;
        prev_NR2 = 0;
        for bs=1:Escenario.num_BS
            BS = pos_scatt(:,prev_BS+1:prev_BS+Escenario.BS_cluster.num_scattBS(bs));                         %extrae dispersores de BS de la celda bs
            rij_temp = BS - kron(Escenario.BS_posicion(:,bs),ones(1,Escenario.BS_cluster.num_scattBS(bs)));        %vector de BS a cada scattFC de esta celda
            rij(:,prev_NR+1:prev_NR+SS.num_rayos(bs)) = kron(ones(1,Escenario.N_d{bs}),rij_temp);             %se repite para cada MS asignado a la BS de esta celda bs
            
            rjk(:,prev_NR+1:prev_NR+SS.num_rayos(bs)) = kron(ones(1,Escenario.N_d{bs}),BS) - ...
                kron(Escenario.MS_posicion(:,Escenario.deseados{bs}),ones(1,Escenario.BS_cluster.num_scattBS(bs)));
                                                                            %vector de MS a cada scattBS de esta celda, en esa dirección para calcular 
                                                                            %directamente los ángulos de interés            
            Tao0(prev_NR+1:prev_NR+SS.num_rayos(bs)) = kron(LoS.Tao(Escenario.deseados{bs}),ones(1,Escenario.BS_cluster.num_scattBS(bs)));
            for ms=1:Escenario.N_d{bs}
                Potencia_temp = Escenario.BS_cluster.BS_Power(bs)*NLoS.P(Escenario.deseados{bs}(ms));
                Potencia(prev_NR2+1:prev_NR2+Escenario.BS_cluster.num_scattBS(bs)) = kron(ones(1,Escenario.BS_cluster.num_scattBS(bs)),Potencia_temp);
                prev_NR2 = prev_NR2+Escenario.BS_cluster.num_scattBS(bs);          
            end
                        
            prev_BS = prev_BS+Escenario.BS_cluster.num_scattBS(bs);
            prev_NR = prev_NR+SS.num_rayos(bs);
        end        
        
    case 'FC'
        %Se trata como si fuera 1 sola nube en cada celda con el total de dispersores dado por la suma de dispersores que hay en todas
        %las nubes lejanas ubicadas en cada celda.
        %IMPORTANTISIMO!!! Este primer paso permitirá preasignar memoria a grandes variables, por eso se justifica utilizar doble
        %iteración!!! y como es variable el # de rayos según los MS asignados a cada celda y sus FC, toca calcularlo antes!
        for bs=1:Escenario.num_BS
            %Número de rayos
            SS.num_rayos(bs) = Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.N_d{bs};        %es el número de dispersores lejanos en esta celda bs por el número de MS asignados a la celda
        end
        %Preasignación de Memoria        
        rij = zeros(3,sum(SS.num_rayos));
        rjk = zeros(3,sum(SS.num_rayos));
        Tao0 = zeros(1,sum(SS.num_rayos));
        Potencia = zeros(1,sum(SS.num_rayos));
        prev_FC = 0;
        prev_NR = 0;
        prev_pot = 0;
        for bs=1:Escenario.num_BS
            FC = pos_scatt(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC_celda(bs));                         %extrae dispersores de FC de la celda bs
            rij_temp = FC - kron(Escenario.BS_posicion(:,bs),ones(1,Escenario.FC_cluster.num_scattFC_celda(bs)));        %vector de BS a cada scattFC de esta celda
            rij(:,prev_NR+1:prev_NR+SS.num_rayos(bs)) = kron(ones(1,Escenario.N_d{bs}),rij_temp);                                                           %se repite para cada MS asignado a la BS de esta celda bs
            
            rjk(:,prev_NR+1:prev_NR+SS.num_rayos(bs)) = kron(ones(1,Escenario.N_d{bs}),FC) - ...
                kron(Escenario.MS_posicion(:,Escenario.deseados{bs}),ones(1,Escenario.FC_cluster.num_scattFC_celda(bs)));
                                                                            %vector de MS a cada scattFC de esta celda, en esa dirección para calcular 
                                                                            %directamente los ángulos de interés            
            Tao0(prev_NR+1:prev_NR+SS.num_rayos(bs)) = kron(LoS.Tao(Escenario.deseados{bs}),ones(1,Escenario.FC_cluster.num_scattFC_celda(bs)));
            
            for ms=1:Escenario.N_d{bs}
                Potencia(prev_pot+1:prev_pot+Escenario.FC_cluster.num_scattFC_celda(bs)) = Escenario.FC_cluster.FC_Power_scatt(prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC_celda(bs))*NLoS.P(Escenario.deseados{bs}(ms));
                prev_pot = prev_pot+Escenario.FC_cluster.num_scattFC_celda(bs);
            end
            prev_FC = prev_FC+Escenario.FC_cluster.num_scattFC_celda(bs);
            prev_NR = prev_NR+SS.num_rayos(bs);
        end
end

%Distancia
distancia_rij = sqrt(rij(1,:).^2+rij(2,:).^2+rij(3,:).^2);
distancia_rjk = sqrt(rjk(1,:).^2+rjk(2,:).^2+rjk(3,:).^2);

%Retardos
Tao = (distancia_rij+distancia_rjk)/3e8;                            %retardo absoluto de cada rayo
SS.Tao = Tao - Tao0;                                                %retardo de exceso de cada rayo
%Fases
Phase = 2*pi*(distancia_rij+distancia_rjk)/Escenario.lamda;         %Fase calculada por Phase=2*pi*dtotal/lamda
SS.Phase = mod(Phase,2*pi);                                         %Para mantener la fase entre [0,2*pi]

%Angulos de Salida
[alpha,phi,radious] = cart2sph(rij(1,:),rij(2,:),rij(3,:));         %por eficiencia es mejor hacerlo así
SS.AoD_azim = alpha*180/pi;
SS.AoD_elev = phi*180/pi;

%Angulos de Llegada
[alphaA,phiA,radiousA] = cart2sph(rjk(1,:),rjk(2,:),rjk(3,:));      %por eficiencia es mejor hacerlo así
SS.AoA_azim = alphaA*180/pi;
SS.AoA_elev = phiA*180/pi;

%Potencia
SS.Potencia = Potencia;
SS.Potencia_dB = 10*log10(Potencia);
%Amplitud
SS.Amplitud = sqrt(SS.Potencia);   

%Para efectos de Double Scattering es más eficiente también extraer las distancias y vectores:
SS.rij = rij;
SS.rjk = rjk;
SS.distancia_rij = distancia_rij;
SS.distancia_rjk = distancia_rjk;

%Cálculo Angle_Spread
switch nube
    case 'MS'
        prev_MS = 0;
        for bs=1:Escenario.num_BS
            for ms=1:Escenario.N_d{bs}
                Pot = SS.Potencia(prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));                
                %AoD Azimutal
                AoD_az = SS.AoD_azim(prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
                temp1 = sum(Pot.*AoD_az.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoD_az)/sum(Pot)).^2;
                SS.A_S_AoD_az(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                %AoD Elevación
                AoD_el = SS.AoD_elev(prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
                temp1 = sum(Pot.*AoD_el.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoD_el)/sum(Pot)).^2;
                SS.A_S_AoD_el(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                %AoA Azimutal
                AoA_az = SS.AoA_azim(prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
                temp1 = sum(Pot.*AoA_az.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoA_az)/sum(Pot)).^2;
                SS.A_S_AoA_az(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                %AoD Elevación
                AoA_el = SS.AoA_elev(prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
                temp1 = sum(Pot.*AoA_el.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoA_el)/sum(Pot)).^2;
                SS.A_S_AoA_el(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);

                prev_MS = prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));
            end
        end
        
    case 'BS'
        prev_BS = 0; prev_MS = 0;
        for bs=1:Escenario.num_BS
            AoD_az = SS.AoD_azim(prev_BS+1:prev_BS+Escenario.BS_cluster.num_scattBS(bs));
            AoD_el = SS.AoD_elev(prev_BS+1:prev_BS+Escenario.BS_cluster.num_scattBS(bs));
            prev_BS = prev_BS+SS.num_rayos(bs);
            for ms=1:Escenario.N_d{bs}
                Pot = SS.Potencia(prev_MS+1:prev_MS+Escenario.BS_cluster.num_scattBS(bs));
                %AoD Azimutal
                temp1 = sum(Pot.*AoD_az.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoD_az)/sum(Pot)).^2;
                SS.A_S_AoD_az(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                %AoD Elevación
                temp1 = sum(Pot.*AoD_el.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoD_el)/sum(Pot)).^2;
                SS.A_S_AoD_el(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                %AoA Azimutal
                AoA_az = SS.AoA_azim(prev_MS+1:prev_MS+Escenario.BS_cluster.num_scattBS(bs));
                temp1 = sum(Pot.*AoA_az.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoA_az)/sum(Pot)).^2;
                SS.A_S_AoA_az(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                %AoD Elevación
                AoA_el = SS.AoA_elev(prev_MS+1:prev_MS+Escenario.BS_cluster.num_scattBS(bs));
                temp1 = sum(Pot.*AoA_el.^2)/sum(Pot);
                temp2 = (sum(Pot.*AoA_el)/sum(Pot)).^2;
                SS.A_S_AoA_el(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);

                prev_MS = prev_MS+Escenario.BS_cluster.num_scattBS(bs);
            end
        end
        
    case 'FC'
        prev_FC = 0; prev_MS = 0; count = 0;
        %Para cada celda!!!
        for bs=1:Escenario.num_BS
            if ~isempty(Escenario.deseados{bs})
                AoD_az = SS.AoD_azim(prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
                AoD_el = SS.AoD_elev(prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
                prev_FC = prev_FC+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.N_d{bs};
                for ms=1:Escenario.N_d{bs}
                    Pot = SS.Potencia(prev_MS+1:prev_MS+Escenario.FC_cluster.num_scattFC_celda(bs));
                    %Para todos los clusters en cada celda!!!                    
                    %AoD Azimutal
                    temp1 = sum(Pot.*AoD_az.^2)/sum(Pot);
                    temp2 = (sum(Pot.*AoD_az)/sum(Pot)).^2;
                    SS.A_S_AoD_az(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                    %AoD Elevación
                    temp1 = sum(Pot.*AoD_el.^2)/sum(Pot);
                    temp2 = (sum(Pot.*AoD_el)/sum(Pot)).^2;
                    SS.A_S_AoD_el(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                    %AoA Azimutal
                    AoA_az = SS.AoA_azim(prev_MS+1:prev_MS+Escenario.FC_cluster.num_scattFC_celda(bs));
                    temp1 = sum(Pot.*AoA_az.^2)/sum(Pot);
                    temp2 = (sum(Pot.*AoA_az)/sum(Pot)).^2;
                    SS.A_S_AoA_az(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);
                    %AoD Elevación
                    AoA_el = SS.AoA_elev(prev_MS+1:prev_MS+Escenario.FC_cluster.num_scattFC_celda(bs));
                    temp1 = sum(Pot.*AoA_el.^2)/sum(Pot);
                    temp2 = (sum(Pot.*AoA_el)/sum(Pot)).^2;
                    SS.A_S_AoA_el(Escenario.deseados{bs}(ms)) = sqrt(temp1-temp2);

                    prev_MS = prev_MS+Escenario.FC_cluster.num_scattFC_celda(bs);
                    
                    %Para cada cluster lejano
                    prev_fc = 0;
                    for fc=1:Escenario.FC_cluster.num_FC(bs)
                        count = count+1;
                        index_fc = sum(Escenario.FC_cluster.num_FC(1:bs))-Escenario.FC_cluster.num_FC(bs);
                        Pot_fc = Pot(prev_fc+1:prev_fc+Escenario.FC_cluster.num_scattFC(index_fc+fc));
                        %AoD Azimutal
                        AoD_az_fc = AoD_az(prev_fc+1:prev_fc+Escenario.FC_cluster.num_scattFC(index_fc+fc));
                        temp1_fc = sum(Pot_fc.*AoD_az_fc.^2)/sum(Pot_fc);
                        temp2_fc = (sum(Pot_fc.*AoD_az_fc)/sum(Pot_fc)).^2;
                        SS.A_S_AoD_az_fc(count) = sqrt(temp1_fc-temp2_fc);
                         %AoD Elevación
                        AoD_el_fc = AoD_el(prev_fc+1:prev_fc+Escenario.FC_cluster.num_scattFC(index_fc+fc));
                        temp1_fc = sum(Pot_fc.*AoD_el_fc.^2)/sum(Pot_fc);
                        temp2_fc = (sum(Pot_fc.*AoD_el_fc)/sum(Pot_fc)).^2;
                        SS.A_S_AoD_el_fc(count) = sqrt(temp1_fc-temp2_fc);
                        %AoA Azimutal
                        AoA_az_fc = AoA_az(prev_fc+1:prev_fc+Escenario.FC_cluster.num_scattFC(index_fc+fc));
                        temp1_fc = sum(Pot_fc.*AoA_az_fc.^2)/sum(Pot_fc);
                        temp2_fc = (sum(Pot_fc.*AoA_az_fc)/sum(Pot_fc)).^2;
                        SS.A_S_AoA_az_fc(count) = sqrt(temp1_fc-temp2_fc);
                         %AoD Elevación
                        AoA_el_fc = AoA_el(prev_fc+1:prev_fc+Escenario.FC_cluster.num_scattFC(index_fc+fc));
                        temp1_fc = sum(Pot_fc.*AoA_el_fc.^2)/sum(Pot_fc);
                        temp2_fc = (sum(Pot_fc.*AoA_el_fc)/sum(Pot_fc)).^2;
                        SS.A_S_AoA_el_fc(count) = sqrt(temp1_fc-temp2_fc);
                      
                        prev_fc = prev_fc+Escenario.FC_cluster.num_scattFC(index_fc+fc);
                    end                    
                end
            end
        end
end

end