%MÓDULO TRANSMISOR

function [TX]=TransmisorSTCyBF()

global Escenario CH;

%GENERACIÓN SÍMBOLOS ENTEROS ALEATORIOS
TX.s_int=randint(Escenario.num_MS,Escenario.tam_trama,[0 Escenario.M-1]);  %Generación de datos a transmitir a cada usuario, 
                                                                            %cada fila contiene la información a tx a cada móvil
%MODULACIÓN
switch Escenario.Modulacion
    case 'PSK'
        TX.Modulador = modem.pskmod('M',Escenario.M,'SYMBOLORDER','gray','INPUTTYPE','integer');    %Objeto Modulador
        if Escenario.M==1
            g=1;
        else
            g=(sin(pi/Escenario.M))^2;
        end
    case 'QAM'
        TX.Modulador = modem.qammod('M',Escenario.M,'SYMBOLORDER','gray','INPUTTYPE','integer');    %Objeto Modulador
        g=3/(2*(Escenario.M-1));
end
TX.modSignal = modulate(TX.Modulador,TX.s_int);                             %Señal Modulada

%CODIFICACIÓN ESPACIO TEMPORAL
[TX.codSignal]=STBCcoder(TX.modSignal);                                     %STBC

%BEAMFORMER
[Uh,Dh,Vh]=svd(CH.Rtx);
Dhp=sum(Dh);

%CARGA EN POTENCIA
for rr=1:Escenario.N
   SNR_Threshold(rr)=(1/g)*(rr/Dhp(rr)-sum(1./Dhp(1:rr)));
end
TX.SNR_Threshold=[SNR_Threshold Inf];

f_2=zeros(1,Escenario.N);
for rr=1:Escenario.N
   if (10^(Escenario.snr/10))>TX.SNR_Threshold(rr) && (10^(Escenario.snr/10))<TX.SNR_Threshold(rr+1)
       Ntp=rr;
       for rrp=1:Ntp
           Dh_temp=1./(Dhp(1:rrp))-1./(Dhp(rrp));
           f_2(rrp)=(1/Ntp+(1/(g*(10^(Escenario.snr/10))))*((1/Ntp)*sum(Dh_temp)));
       end
   end
end
% f_2=[1 zeros(Escenario.N-1)];                                               %1D EigenBF                                                               
TX.Uh=Uh;
TX.f=sqrt(f_2);
TX.codSignal=conj(TX.codSignal')*diag(sqrt(f_2))*Uh';
TX.codSignal=conj(TX.codSignal');

%FILTRADO
[TX.param_tx,s_tx,t_tx]=filtrorootraisedcosineTx(TX.codSignal);             %Filtro Transmisor
TX.s_tx=conj(s_tx');