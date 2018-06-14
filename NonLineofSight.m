%Esta función calcula Ptot o PNLoS, que corresponde a la potencia por Path
%Loss y además incluye también el fenómeno de Shadowing

function [NLoS]=NonLineofSight()

global Escenario

distancias_rij = sqrt(Escenario.rij(1,:).^2+Escenario.rij(2,:).^2+Escenario.rij(3,:).^2);   %distancias de cada BS a sus MS asignados

if Escenario.f <1e9
    %COST231-Hata
    %Aplica para 150<f<1000 MHz, 1<d<20 Km, 30<hbs<200 m, 1<hms<10m.
    NLoS.id = 'COST231-Hata';
    a_hmobile = (1.1*log10(Escenario.f/1e6)-0.7)*Escenario.h_MS-(1.56*log10(Escenario.f/1e6)-0.8);
    NLoS.P_dB = -(69.55+26.16*log10(Escenario.f/1e6)-13.82*log10(Escenario.h_BS)-a_hmobile+(44.9-6.55*log10(Escenario.h_BS))*log10(distancias_rij/1000));
elseif 1.5e9<Escenario.f<=2e9
    %COST-Hata
    %Aplica para 1500<f<2000 MHz, 1<d<20 Km, 30<hbs<200 m, 1<hms<10m.
    NLoS.id = 'COST-Hata';
    a_hmobile = (1.1*log10(Escenario.f/1e6)-0.7)*Escenario.h_MS-(1.56*log10(Escenario.f/1e6)-0.8);
    Cm = 3;     %0 dB para ciudades medianas & suburbanas, 3 dB para centros metropolitanos
    NLoS.P_dB = -(46.3+33.9*log10(Escenario.f/1e6)-13.82*log10(Escenario.h_BS)-a_hmobile+(44.9-6.55*log10(Escenario.h_BS))*log10(distancias_rij/1000)+Cm);    
end

%Shadowing
shad_dB = Escenario.Var_Shadowing_dB*randn;

%Pérdidas Totales en dB
NLoS.P_dB = NLoS.P_dB + shad_dB;
NLoS.P = 10.^((NLoS.P_dB)/10);