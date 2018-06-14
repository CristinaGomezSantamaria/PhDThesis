%MÓDULO RECEPTOR

function [RX]=Receptor()

global Escenario CH TX

% %FILTRO RECEPTOR
% [RX.s_out_filt,s_defilt,t_defilt]=filtrorootraisedcosineRx();               %Filtro Receptor

RX.s_out_filt = CH.x;
%Scale factor for fear comparison with constellation symbols
RX.s_out_filt = RX.s_out_filt/(TX.k_QAM*TX.PL);

%DECODIFICADOR
[RX.decodSignal]=STBCdecoder(RX.s_out_filt);                                %STBC

%DEMODULACIÓN
%Demodulator
RX.Demodulador = modem.qamdemod('M',Escenario.M,'SYMBOLORDER','gray','OUTPUTTYPE','integer');  %Objeto Demodulador
RX.demodSignal = demodulate(RX.Demodulador,RX.decodSignal);                 %Señal Demodulada

end