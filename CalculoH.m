%Esta función realiza lo siguiente:
%- Genera la posición de los Dispersores en cada nube dispersora, según los parámetros generales dados en la Plantilla de Configuración
%- Calcula los parámetros de los rayos
%- Normaliza la contribución de cada rayo y suma para calcular h(tao,AoD,AoA)
%- Calcula la matriz del canal H
%Su entrada es LoS pues esta se calcula 1 sola vez por escenario para optimizar,
%pero debe tenerse en cuenta en todos los cálculos de H para sumar sus contribuciones

function [H,h]=CalculoH(LoS)

global Escenario snr
tic;

%PATH LOSS CALCULATION
[NLoS]=NonLineofSight();

%SINGLE SCATTERING
if strcmp(Escenario.MS_cluster.id,'si')
    %Generating the specific positions of scatterers
    h.posicion_scattMS = zeros(3,sum(Escenario.MS_cluster.num_scattMS));
    prev_MS = 0;
    for bs=1:Escenario.num_BS
        for ii=1:Escenario.N_d{bs}
            switch Escenario.MS_cluster.pdf_scatt(ii,:)
                case 'gauss_cyl'
                    h.posicion_scattMS(:,prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ii))) = ...
                        gauss_cyl_AR(Escenario.MS_cluster.parametros{Escenario.deseados{bs}(ii)},Escenario.MS_posicion(:,Escenario.deseados{bs}(ii)),...
                        Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ii)));
                case 'gauss'
                    h.posicion_scattMS(:,prev_MS+1:prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ii))) = ...
                        gauss(Escenario.MS_cluster.parametros{Escenario.deseados{bs}(ii)},Escenario.MS_posicion(:,Escenario.deseados{bs}(ii)),...
                        Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ii)));
            end
            prev_MS = prev_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ii));
        end
    end
    h.posicion_scattMS(3,find(h.posicion_scattMS(3,:)<0)) = 0;
    %Calculating physical parameters for Single Scattering
    [SS_MS]=SingleScattering(h.posicion_scattMS,'MS',LoS,NLoS);    
end

if strcmp(Escenario.BS_cluster.id,'si')
    %Generating the specific positions of scatterers
    h.posicion_scattBS = zeros(3,sum(Escenario.BS_cluster.num_scattBS));            
    prev_BS = 0;
    for ii=1:Escenario.num_BS
        switch Escenario.BS_cluster.pdf_scatt(ii,:)
            case 'gauss_cyl'
            h.posicion_scattBS(:,prev_BS+1:prev_BS+Escenario.BS_cluster.num_scattBS(ii)) = ...
                gauss_cyl_AR(Escenario.BS_cluster.parametros{ii},Escenario.BS_posicion(:,ii),Escenario.BS_cluster.num_scattBS(ii));                
            case 'gauss'
               h.posicion_scattBS(:,prev_BS+1:prev_BS+Escenario.BS_cluster.num_scattBS(ii)) = ...
                gauss(Escenario.BS_cluster.parametros{ii},Escenario.BS_posicion(:,ii),Escenario.BS_cluster.num_scattBS(ii));                             
        end
        prev_BS = prev_BS+Escenario.BS_cluster.num_scattBS(ii);
    end
    h.posicion_scattBS(3,find(h.posicion_scattBS(3,:)<0)) = 0;
    %Calculating physical parameters for Single Scattering
    [SS_BS]=SingleScattering(h.posicion_scattBS,'BS',LoS,NLoS);
end

if strcmp(Escenario.FC_cluster.id,'si')
    %Generating the specific positions of scatterers    
    h.posicion_scattFC = zeros(3,sum(Escenario.FC_cluster.num_scattFC));
    prev_FC = 0;
    for ii=1:sum(Escenario.FC_cluster.num_FC)
        switch Escenario.FC_cluster.pdf_scatt(ii,:)
            case 'gauss_cyl'
                h.posicion_scattFC(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC(ii)) = ...
                    gauss_cyl_AR(Escenario.FC_cluster.parametros{ii},Escenario.FC_cluster.posicion_FC(:,ii),...
                    Escenario.FC_cluster.num_scattFC(ii));
            case 'gauss'
                h.posicion_scattFC(:,prev_FC+1:prev_FC+Escenario.FC_cluster.num_scattFC(ii)) = ...
                    gauss(Escenario.FC_cluster.parametros{ii},Escenario.FC_cluster.posicion_FC(:,ii),...
                    Escenario.FC_cluster.num_scattFC(ii));                
        end
        prev_FC = prev_FC+Escenario.FC_cluster.num_scattFC(ii);
    end
    %Calculating physical parameters for Single Scattering
    [SS_FC]=SingleScattering(h.posicion_scattFC,'FC',LoS,NLoS);    
end

%DOUBLE SCATTERING
%So far it is implemented only between Far Scatterers and Scatterers around MS because this is the process we are interested on... in the future
%extend to include the double scattering between Far Scatterers and Scatterers around BS, and between scatterers around BS and MS
[DS_FC_MS]=DoubleScattering(SS_FC,SS_MS,LoS,NLoS,h);

%CALCULATION OF h(Tao,AoD,AoA) FOR EACH LINK BETWEEN A BS AND ITS ASSIGNED MS's --> IT DOESN´T INCLUDE YET THE LoS COMPONENT
%This calculation still doesn´t collect all the rays with similar properties (angles and delays), precisely because it is interesting to observe
%its variations

%Each Double Directional Channel Impulse Response h{bs,ms}(Tao,AoD,AoA) between a BS and its assigned MS's is
%collected in a matrix where each row corresponds to a physical parameter of the channel
count_MS = 0; count_BS = 0; count_FC = 0; count_DS_FC_MS= 0;
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
    for ms=1:Escenario.N_d{bs}
        count = 0;
        %Single Scattering Contributions
        if strcmp(Escenario.Single_Scatt,'si')
        if strcmp(Escenario.MS_cluster.id,'si')
            h.DDCIR{bs,ms}(1,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.Tao(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            h.DDCIR{bs,ms}(2,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.Phase(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            h.DDCIR{bs,ms}(3,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.AoD_azim(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            h.DDCIR{bs,ms}(4,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.AoD_elev(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            h.DDCIR{bs,ms}(5,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.AoA_azim(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            h.DDCIR{bs,ms}(6,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.AoA_elev(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            h.DDCIR{bs,ms}(7,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.Potencia(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            h.DDCIR{bs,ms}(8,count+1:count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
                SS_MS.Amplitud(count_MS+1:count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
            
            count = count+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));
            count_MS = count_MS+Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));
        end
        if strcmp(Escenario.BS_cluster.id,'si')        
            h.DDCIR{bs,ms}(1,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.Tao(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            h.DDCIR{bs,ms}(2,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.Phase(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            h.DDCIR{bs,ms}(3,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.AoD_azim(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            h.DDCIR{bs,ms}(4,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.AoD_elev(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            h.DDCIR{bs,ms}(5,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.AoA_azim(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            h.DDCIR{bs,ms}(6,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.AoA_elev(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            h.DDCIR{bs,ms}(7,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.Potencia(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            h.DDCIR{bs,ms}(8,count+1:count+Escenario.BS_cluster.num_scattBS(bs)) = ...
                SS_BS.Amplitud(count_BS+1:count_BS+Escenario.BS_cluster.num_scattBS(bs));
            
            count = count+Escenario.BS_cluster.num_scattBS(bs);
            count_BS = count_BS+Escenario.BS_cluster.num_scattBS(bs);
        end
        if strcmp(Escenario.FC_cluster.id,'si')        
            h.DDCIR{bs,ms}(1,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.Tao(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            h.DDCIR{bs,ms}(2,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.Phase(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            h.DDCIR{bs,ms}(3,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.AoD_azim(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            h.DDCIR{bs,ms}(4,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.AoD_elev(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            h.DDCIR{bs,ms}(5,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.AoA_azim(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            h.DDCIR{bs,ms}(6,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.AoA_elev(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            h.DDCIR{bs,ms}(7,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.Potencia(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            h.DDCIR{bs,ms}(8,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)) = ...
                SS_FC.Amplitud(count_FC+1:count_FC+Escenario.FC_cluster.num_scattFC_celda(bs));
            
            count = count+Escenario.FC_cluster.num_scattFC_celda(bs);
            count_FC = count_FC+Escenario.FC_cluster.num_scattFC_celda(bs);
        end
        end
        %Double Scattering Contributions
        if strcmp(Escenario.Double_Scatt,'si')
        h.DDCIR{bs,ms}(1,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.Tao(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
        h.DDCIR{bs,ms}(2,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.Phase(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
        h.DDCIR{bs,ms}(3,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.AoD_azim(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
        h.DDCIR{bs,ms}(4,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.AoD_elev(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
        h.DDCIR{bs,ms}(5,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.AoA_azim(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
        h.DDCIR{bs,ms}(6,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.AoA_elev(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
        h.DDCIR{bs,ms}(7,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.Potencia(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));
        h.DDCIR{bs,ms}(8,count+1:count+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms))) = ...
            DS_FC_MS.Amplitud(count_DS_FC_MS+1:count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms)));

        count_DS_FC_MS = count_DS_FC_MS+Escenario.FC_cluster.num_scattFC_celda(bs)*Escenario.MS_cluster.num_scattMS(Escenario.deseados{bs}(ms));
        end
        %Calculation of General Angular Dispersion for each link
        %AoD Azimuthal
        a_2 = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(3,:).^2)./sum(h.DDCIR{bs,ms}(7,:));
        a = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(3,:))./sum(h.DDCIR{bs,ms}(7,:));
        h.A_S_AoD_az(Escenario.deseados{bs}(ms)) = sqrt(a_2-a.^2);
        %AoD Elevation
        a_2 = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(4,:).^2)./sum(h.DDCIR{bs,ms}(7,:));
        a = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(4,:))./sum(h.DDCIR{bs,ms}(7,:));
        h.A_S_AoD_el(Escenario.deseados{bs}(ms)) = sqrt(a_2-a.^2);
        %AoA Azimuthal
        a_2 = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(5,:).^2)./sum(h.DDCIR{bs,ms}(7,:));
        a = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(5,:))./sum(h.DDCIR{bs,ms}(7,:));
        h.A_S_AoA_az(Escenario.deseados{bs}(ms)) = sqrt(a_2-a.^2);
        %AoA Elevation
        a_2 = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(6,:).^2)./sum(h.DDCIR{bs,ms}(7,:));
        a = sum(h.DDCIR{bs,ms}(7,:).*h.DDCIR{bs,ms}(6,:))./sum(h.DDCIR{bs,ms}(7,:));
        h.A_S_AoA_el(Escenario.deseados{bs}(ms)) = sqrt(a_2-a.^2);
    end
    end
end
%Extracting the results of Angle Dispersion for each cluster in each link
if strcmp(Escenario.MS_cluster.id,'si')
    h.A_S_AoD_az_ms = SS_MS.A_S_AoD_az;
    h.A_S_AoD_el_ms = SS_MS.A_S_AoD_el;
    h.A_S_AoA_az_ms = SS_MS.A_S_AoA_az;
    h.A_S_AoA_el_ms = SS_MS.A_S_AoA_el;
end
if strcmp(Escenario.BS_cluster.id,'si')
    h.A_S_AoD_az_bs = SS_BS.A_S_AoD_az;
    h.A_S_AoD_el_bs = SS_BS.A_S_AoD_el;
    h.A_S_AoA_az_bs = SS_BS.A_S_AoA_az;
    h.A_S_AoA_el_bs = SS_BS.A_S_AoA_el;
end
if strcmp(Escenario.FC_cluster.id,'si')
    h.A_S_AoD_az_fc = SS_FC.A_S_AoD_az_fc;
    h.A_S_AoD_el_fc = SS_FC.A_S_AoD_el_fc;
    h.A_S_AoA_az_fc = SS_FC.A_S_AoA_az_fc;
    h.A_S_AoA_el_fc = SS_FC.A_S_AoA_el_fc;
end

%CALCULATION OF H(tao) AND PtotalNLoS, IT CALCULATES WITH PtotalNLoS AND K_RICE, THE PotLoS TO GENERATE THE LINE OF SIGHT COMPONENT
H.Tsymbol = 5e-6;           %Parámetro por definir según resolución y tipo del canal, para asegurar que sea Banda Ancha, y que no haya demasiados
                            %taps en el procesamiento requerido
H.Max_Retardo = 20e-6;      %Debe ser seleccionado también según el tipo de canal, tal que se asegure que se tendrán en cuenta todos los posibles
                            %retardos del canal
H.Limite_Taps = H.Max_Retardo/H.Tsymbol;
Pot_totalNLoS = zeros(1,Escenario.num_MS);
num_disp_lt = zeros(H.Limite_Taps);
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
    for ms=1:Escenario.N_d{bs}
        H.Tap{bs,ms}=zeros(Escenario.P,Escenario.N*H.Limite_Taps);
        for lt=1:H.Limite_Taps
            ind_tap_lt = find(h.DDCIR{bs,ms}(1,:)>=((lt-1)*H.Tsymbol) & h.DDCIR{bs,ms}(1,:)<(lt*H.Tsymbol));
            num_disp_lt(lt) = length(ind_tap_lt);
            k_BS = 2*pi*[cos(h.DDCIR{bs,ms}(3,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(4,ind_tap_lt)*pi/180); ...
                      sin(h.DDCIR{bs,ms}(3,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(4,ind_tap_lt)*pi/180); ...
                      sin(h.DDCIR{bs,ms}(4,ind_tap_lt)*pi/180)].';
            AMV_BS = kron(ones(Escenario.N,1),h.DDCIR{bs,ms}(8,ind_tap_lt)).*exp(-j*k_BS*Escenario.BS_array)';
            k_MS = 2*pi*[cos(h.DDCIR{bs,ms}(5,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(6,ind_tap_lt)*pi/180); ...
                      sin(h.DDCIR{bs,ms}(5,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(6,ind_tap_lt)*pi/180); ...
                      sin(h.DDCIR{bs,ms}(6,ind_tap_lt)*pi/180)].';
            AMV_MS = kron(ones(Escenario.P,1),exp(j*h.DDCIR{bs,ms}(2,ind_tap_lt))).*exp(-j*k_MS*Escenario.MS_array)';
            AMV_BS_temp = kron(ones(Escenario.P,1),reshape(AMV_BS,1,[]));
            AMV_MS_temp = kron(AMV_MS,ones(1,Escenario.N));
            ch_temp = AMV_BS_temp.*AMV_MS_temp;
            suma = zeros(Escenario.P,Escenario.N);
            for nn=1:Escenario.N
                index = nn:Escenario.N:length(ind_tap_lt)*Escenario.N;
                for pp=1:Escenario.P
                    suma(pp,nn) = sum(ch_temp(pp,index));
                end
            end
            H.Tap{bs,ms}(:,(lt-1)*Escenario.N+1:lt*Escenario.N) = suma;
            %Calculation of the Total Power provided by the scatterers (SS & DS)
            P_lt = (norm(H.Tap{bs,ms}(:,(lt-1)*Escenario.N+1:lt*Escenario.N),'fro'))^2;
            Pot_totalNLoS(Escenario.deseados{bs}(ms)) = Pot_totalNLoS(Escenario.deseados{bs}(ms))+P_lt;
            %For calculating the Transmit Covariance Matrix
            AMVBS{lt} = (sum(AMV_BS.')).';
            P_Rtx_lt(lt) = (norm(AMVBS{lt},'fro'))^2;
        end
        for lt=1:H.Limite_Taps
            H.AMVBS_norm{Escenario.deseados{bs}(ms),lt} = AMVBS{lt}/sqrt(sum(P_Rtx_lt));
        end
        %Calculation of the Contribution provided by the Line of Sight Component and inclusion in the Double Directional Channel Impulse
        %Response h{bs,ms}(Tao,AoD,AoA)
        if strcmp(Escenario.LoS,'si')
        LoS_temp = [0; LoS.Phase(Escenario.deseados{bs}(ms)); LoS.AoD_azim(Escenario.deseados{bs}(ms)); LoS.AoD_elev(Escenario.deseados{bs}(ms));...
            LoS.AoA_azim(Escenario.deseados{bs}(ms)); LoS.AoA_elev(Escenario.deseados{bs}(ms)); ...
            Escenario.K_Rice*Pot_totalNLoS(Escenario.deseados{bs}(ms)); sqrt(Escenario.K_Rice*Pot_totalNLoS(Escenario.deseados{bs}(ms)))];
        h.DDCIR{bs,ms} = [LoS_temp h.DDCIR{bs,ms}];
        end
    end
    end
end

%Including the LoS component in the final calculation of the first tap of H(tao)
if strcmp(Escenario.LoS,'si')
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
    for ms=1:Escenario.N_d{bs}
        pot = h.DDCIR{bs,ms}(7,1);
        k_BS = 2*pi*[cos(h.DDCIR{bs,ms}(3,1)*pi/180).*cos(h.DDCIR{bs,ms}(4,1)*pi/180); ...
                     sin(h.DDCIR{bs,ms}(3,1)*pi/180).*cos(h.DDCIR{bs,ms}(4,1)*pi/180); ...
                     sin(h.DDCIR{bs,ms}(4,1)*pi/180)].';
        AMV_BS = kron(ones(Escenario.N,1),sqrt(pot)).*exp(-j*k_BS*Escenario.BS_array)';
        k_MS = 2*pi*[cos(h.DDCIR{bs,ms}(5,1)*pi/180).*cos(h.DDCIR{bs,ms}(6,1)*pi/180); ...
                     sin(h.DDCIR{bs,ms}(5,1)*pi/180).*cos(h.DDCIR{bs,ms}(6,1)*pi/180); ...
                     sin(h.DDCIR{bs,ms}(6,1)*pi/180)].';
        AMV_MS = kron(ones(Escenario.P,1),exp(j*h.DDCIR{bs,ms}(2,1))).*exp(-j*k_MS*Escenario.MS_array)';      
        ch_temp = AMV_MS*AMV_BS.';
        H.Tap{bs,ms}(:,1:Escenario.N) = H.Tap{bs,ms}(:,1:Escenario.N)+ch_temp;
        AMVBS{1} = AMVBS{1}+AMV_BS;
        P_Rtx_lt(1) = (norm(AMVBS{1},'fro'))^2;
        for lt=1:H.Limite_Taps
            H.AMVBS_norm{Escenario.deseados{bs}(ms),lt} = AMVBS{lt}/sqrt(sum(P_Rtx_lt));
        end
    end
    end
end
end

%Calculating the Delay Spread for each link between a BS and its assigned MS's
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
    for ms=1:Escenario.N_d{bs}
        %Sin Normalización
        P_lt = zeros(1,H.Limite_Taps); 
        num_t_2 = zeros(1,H.Limite_Taps); 
        num_t = zeros(1,H.Limite_Taps);
        for lt=1:H.Limite_Taps
            P_lt(lt) = (norm(H.Tap{bs,ms}(:,(lt-1)*Escenario.N+1:lt*Escenario.N),'fro'))^2;
            num_t_2(lt) = P_lt(lt)*(lt*H.Tsymbol)^2;
            num_t(lt) = P_lt(lt)*(lt*H.Tsymbol);
        end
        t_2 = sum(num_t_2)/sum(P_lt);
        t = sum(num_t)/sum(P_lt);
        H.D_S(Escenario.deseados{bs}(ms)) = sqrt(t_2-t.^2);
        %Con Normalización en Potencia
        num_t_2_norm = zeros(1,H.Limite_Taps); 
        num_t_norm = zeros(1,H.Limite_Taps);
        P_lt_norm = P_lt./sum(P_lt);
        for lt=1:H.Limite_Taps
            num_t_2_norm(lt) = P_lt_norm(lt)*(lt*H.Tsymbol)^2;
            num_t_norm(lt) = P_lt_norm(lt)*(lt*H.Tsymbol);
        end
        t_2_norm = sum(num_t_2_norm)/sum(P_lt_norm);
        t_norm = sum(num_t_norm)/sum(P_lt_norm);
        H.D_S_norm(Escenario.deseados{bs}(ms)) = sqrt(t_2_norm-t_norm.^2);
        %Normalización del Canal
        H.Tap_norm{bs,ms}=H.Tap{bs,ms}/sqrt(sum(P_lt));
        %Con Normalización en Canal
        P_lt_norm2 = zeros(1,H.Limite_Taps); 
        num_t_2_norm2 = zeros(1,H.Limite_Taps); 
        num_t_norm2 = zeros(1,H.Limite_Taps);
        for lt=1:H.Limite_Taps
            P_lt_norm2(lt) = (norm(H.Tap_norm{bs,ms}(:,(lt-1)*Escenario.N+1:lt*Escenario.N),'fro'))^2;
            num_t_2_norm2(lt) = P_lt_norm2(lt)*(lt*H.Tsymbol)^2;
            num_t_norm2(lt) = P_lt_norm2(lt)*(lt*H.Tsymbol);
        end
        t_2_norm2 = sum(num_t_2_norm2)/sum(P_lt_norm2);
        t_norm2 = sum(num_t_norm2)/sum(P_lt_norm2);
        H.D_S_norm2(Escenario.deseados{bs}(ms)) = sqrt(t_2_norm2-t_norm2.^2);
        %Introduciendo el AWGN por cada tap en los elementos de antena receptores
        n = kron(1/sqrt(2)*(randn(Escenario.P,H.Limite_Taps) + j*randn(Escenario.P,H.Limite_Taps)),ones(1,Escenario.N));
        SNR = 10^(snr/10); % SNR in linear scale
        variance = 1/(2*SNR); % Variance
        sigma = sqrt(variance); % Standard Deviation
        H.Tap_norm_awgn{bs,ms} = H.Tap_norm{bs,ms} + sigma*n;
    end
    end
end
%toc;
end


