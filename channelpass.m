%CANAL

function [CH]=channelpass()

global Escenario CH TX bs

for ms=1:Escenario.N_d{bs}
    sig = sqrt(0.5/(10^(Escenario.snr/10)));                                              %Varianza del Ruido
    CH.n = sig*(randn(Escenario.P,size(TX.s_tx,2))+ j*randn(Escenario.P,size(TX.s_tx,2)));                    %Ruido
    CH.x=CH.hS{ms,1}*TX.s_tx + CH.n;                                                         %Señal recibida
end

