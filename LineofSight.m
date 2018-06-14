%Esta función calcula los parámetros de los rayos de Línea de Vista

function [LoS]=LineofSight()

global Escenario

distancias_rij = sqrt(Escenario.rij(1,:).^2+Escenario.rij(2,:).^2+Escenario.rij(3,:).^2);   %distancias de cada BS a sus MS asignados
LoS.num_rayos = Escenario.num_MS;
% LoS.Amplitud = 1./(4*pi*distancias_rij./Escenario.lamda);                    %Cálculo como a=1/(4*pi*d/lamda)
% LoS.Potencia = LoS.Amplitud.^2;
LoS.Tao = distancias_rij/3e8;                                               %Retardo de cada BS a sus MS asignados
LoS.Phase = 2*pi*distancias_rij/Escenario.lamda;                            %Fase calculada por Phase=2*pi*d/lamda
LoS.Phase = mod(LoS.Phase,2*pi);                                            %Para mantener la fase entre [0,2*pi]
LoS.AoD_elev = asin(Escenario.rij(3,:)./distancias_rij)*180/pi;             %Ang Elevación en extremo BS, no require mapeo, da entre [-90°,90°]
LoS.AoA_elev = -LoS.AoD_elev;                                               %Ang Elevación en extremo MS, como r_bs_ms=-r_ms_bs, este ángulo es sólo la negación
LoS.AoD_azim = atan(Escenario.rij(2,:)./Escenario.rij(1,:))*180/pi;         %Ang Azimutal en extremo BS, si requiere mapeo para que de entre [0°,360°]
for ms=1:Escenario.num_MS
    if Escenario.rij(1,ms)<0 && Escenario.rij(2,ms)<0
        LoS.AoD_azim(ms) = LoS.AoD_azim(ms)+180;
        LoS.AoA_azim(ms) = LoS.AoD_azim(ms)-180;                            %Como r_bs_ms=-r_ms_bs, este ángulo queda rotado 180° según el cuadrante                       
    elseif Escenario.rij(1,ms)<0 && Escenario.rij(2,ms)>0
        LoS.AoD_azim(ms) = LoS.AoD_azim(ms)+180;
        LoS.AoA_azim(ms) = LoS.AoD_azim(ms)+180;
    elseif Escenario.rij(1,ms)>0 && Escenario.rij(2,ms)<0
        LoS.AoD_azim(ms) = 360+LoS.AoD_azim(ms);
        LoS.AoA_azim(ms) = LoS.AoD_azim(ms)-180;
    else
        LoS.AoA_azim(ms) = LoS.AoD_azim(ms)+180;        
    end
end
%[alpha,phi,radious] = cart2sph(Escenario.rij(1,:),Escenario.rij(2,:),Escenario.rij(3,:));
end