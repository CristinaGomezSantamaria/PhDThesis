%Esta función calcula los parámetros de los rayos que pasan por un proceso
%de Single Scattering

function [DS]=DoubleScattering(SS1,SS2,LoS,NLoS,h)

global Escenario    

for bs=1:Escenario.num_BS
    DS.num_rayos(bs) = Escenario.FC_cluster.num_scattFC_celda(bs)*sum(Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}));
end

%Preasignación de Memoria        
rij = zeros(3,sum(DS.num_rayos));
rjk = zeros(3,sum(DS.num_rayos));
Tao0 = zeros(1,sum(DS.num_rayos));
Potencia = zeros(1,sum(DS.num_rayos));
prev_NR = 0;
prev_NR2 = 0;
prev_NR3 = 0;
prev_NR4 = 0;
prev_int = 0;
prev_pos1 = 0;
prev_pos2 = 0;
prev_pot = 0;
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
        rij_temp = SS1.rij(:,prev_NR2+1:prev_NR2+SS1.num_rayos(bs));
        length_tramobasico = Escenario.FC_cluster.num_scattFC_celda(bs);
        pos_FC_temp = h.posicion_scattFC(:,prev_pos1+1:prev_pos1+Escenario.FC_cluster.num_scattFC_celda(bs));
        num_repeticiones = 0;
        length_MS = sum(Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}));
        rjk(:,prev_NR4+1:prev_NR4+length_MS*Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
            kron(SS2.rjk(:,prev_NR3+1:prev_NR3+length_MS),ones(1,Escenario.FC_cluster.num_scattFC_celda(bs)));
        pos_MS(:,prev_NR4+1:prev_NR4+length_MS*Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
            kron(h.posicion_scattMS(:,prev_NR3+1:prev_NR3+length_MS),ones(1,Escenario.FC_cluster.num_scattFC_celda(bs)));
        for ms=1:Escenario.N_d{bs}
            num_repeticiones = num_repeticiones+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));
            Tao0(prev_int+1:prev_int+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))*Escenario.FC_cluster.num_scattFC_celda(bs)) =...
                kron(LoS.Tao(Escenario.deseados{bs}(ms)),ones(1,Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))*...
                Escenario.FC_cluster.num_scattFC_celda(bs)));
            prev_int = prev_int+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))*Escenario.FC_cluster.num_scattFC_celda(bs);
        end
        rij(:,prev_NR+1:prev_NR+num_repeticiones*length_tramobasico) = kron(ones(1,num_repeticiones),rij_temp(:,1:length_tramobasico));    
        pos_FC(:,prev_pos2+1:prev_pos2+Escenario.FC_cluster.num_scattFC_celda(bs)*num_repeticiones) = ...
            kron(ones(1,num_repeticiones),pos_FC_temp);
        
        for ms=1:Escenario.N_d{bs}
            Potencia_temp = Escenario.FC_cluster.FC_Power_scatt(prev_pos1+1:prev_pos1+Escenario.FC_cluster.num_scattFC_celda(bs))*NLoS.P(Escenario.deseados{bs}(ms));
            Potencia(prev_pot+1:prev_pot+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                kron(ones(1,Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))),Potencia_temp)*Escenario.FC_cluster.FC_PowerDS;
            prev_pot = prev_pot+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));
        end
        prev_NR = prev_NR+num_repeticiones*length_tramobasico; 
        prev_NR2 = prev_NR2+SS1.num_rayos(bs);
        prev_NR3 = prev_NR3+length_MS;
        prev_NR4 = prev_NR4+length_MS*Escenario.FC_cluster.num_scattFC_celda(bs);
        prev_pos1 = prev_pos1+Escenario.FC_cluster.num_scattFC_celda(bs);
        prev_pos2 = prev_pos2+Escenario.FC_cluster.num_scattFC_celda(bs)*num_repeticiones;
    end
end

%Vector entre dispersores del FC y dispersores del MS
rjj = pos_MS-pos_FC;
distancia_rjj = sqrt(rjj(1,:).^2+rjj(2,:).^2+rjj(3,:).^2);

%Distancia
distancia_rij = sqrt(rij(1,:).^2+rij(2,:).^2+rij(3,:).^2);
distancia_rjk = sqrt(rjk(1,:).^2+rjk(2,:).^2+rjk(3,:).^2);

%Retardos
Tao = (distancia_rij+distancia_rjj+distancia_rjk)/3e8;              %retardo absoluto de cada rayo
DS.Tao = Tao - Tao0;                                                %retardo de exceso de cada rayo
%Fases
Phase = 2*pi*(distancia_rij+distancia_rjj+distancia_rjk)/Escenario.lamda;         %Fase calculada por Phase=2*pi*dtotal/lamda
DS.Phase = mod(Phase,2*pi);                                         %Para mantener la fase entre [0,2*pi]

%Angulos de Salida
[alpha,phi,radious] = cart2sph(rij(1,:),rij(2,:),rij(3,:));         %por eficiencia es mejor hacerlo así
DS.AoD_azim = alpha*180/pi;
DS.AoD_elev = phi*180/pi;

%Angulos de Llegada
[alphaA,phiA,radiousA] = cart2sph(rjk(1,:),rjk(2,:),rjk(3,:));      %por eficiencia es mejor hacerlo así
DS.AoA_azim = alphaA*180/pi;
DS.AoA_elev = phiA*180/pi;

%Potencia
DS.Potencia = Potencia;
DS.Potencia_dB = 10*log10(Potencia);
%Amplitud
DS.Amplitud = sqrt(DS.Potencia);   

end