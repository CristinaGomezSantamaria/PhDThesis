%BF+STBC+PL SU FOR FREQUENCY SELECTIVE MACROCELL CHANNELS
% Simulator with beamforming, power loading, space-time coding and preequalization
% for macrocell outdoor channels with frequency selectivity
% This simulator includes the following cases:

% 1- Channel with no Far clusters, flat fading wireless channel. It collapses to the proposal
% presented in:
% Shengli Zhou and G. B. Giannakis, "Optimal transmitter eigen-beamforming and space-time block coding based on channel correlations," 
% in IEEE Transactions on Information Theory, vol. 49, no. 7, pp.
% 1673-1690, July 2003. doi: 10.1109/TIT.2003.813565

% 2- Channel with same conditions as 1 and additionally with low
% correlation at the transmitter, it collapses to the only STC and Equal
% Power Loading (EPL) case

% 3- Channel with same conditions as 1 and high correlation at the
% transmitter, it collapses to the 1Dimensional Beamforming (1D-BF) case

% 4- Channel as 1 and with high SNR, it collapses to STC (CSITx not required) and EPL

% 5- Channel as 1 and low SNR, it collapses to the 1D-BF case

% 6- Analysis of the equivalent cases when the channel has frequency selectivity

% The goal in this simulations is actually to analyze the case 6

close all; clear all; clc;

global Escenario TX CH RX bs
tic;
kk=4;
[Escenario]=Plantilla();                                                    %Scenario setup, antenna arrays and general parameters of the channel
[LoS]=LineofSight();                                                        %Cálculo de rayos que van por via directa, se hace en este punto pues
                                                                            %solo depende de Escenario, por tanto debe calcularse 1 sola vez
kk=4;
while(kk)
    THRESHOLD = 100;                                                        % Minimum number of errors before proceeding to next SNR
    MAXITER = 100;                                                          % Maximum number of iterrations before quiting
    snr_vec=0:2:30;
    BitErr=zeros(1,length(snr_vec));                                        % Array for counting biterrors
    SerErr=zeros(1,length(snr_vec));                                        % Array for counting sererrors
    bs=1;
    filevar = ['Eb_N0_dB'; 'maxEb   '; 'Ber     '; 'Ber4    '; 'Ber16   '; 'Ber64   '];
    
    PL=zeros(length(snr_vec),Escenario.N);
    for ii=1:length(snr_vec);
        Escenario.snr=snr_vec(ii);
        countiter = 0;
        fprintf('\n At SNR = %d of %s \n',snr_vec(ii),strtrim(filevar(kk+2,:)));    
        while(THRESHOLD>BitErr(ii) && countiter<MAXITER)
            countiter = countiter+1;
            [H,h]=CalculoH(LoS);                                                 % Generates Scatterers, and Calculates Channel Matrix
            
            [CH]=channel_gen();
            [TX]=TransmisorSTCyBF();                                                  %Módulo generador de señales a transmitir
            [CH]=channelpass();
            [RX]=Receptor();                                                    %Módulo procesador de señales recibidas
            PL(snr_it,:)=PL(snr_it,:)+TX.f;
            if it==1
                BF=TX.Uh;
                SNR_Threshold=TX.SNR_Threshold;
            end
            %Counting Bit Error
            errors = sum(sum(abs(dec2bin(TX.s_int)-dec2bin(RX.demodSignal))));
            BitErr(ii) = BitErr(ii) + errors;
            errors_ser = length(find(TX.s_int~=RX.demodSignal));
            SerErr(ii) = SerErr(ii) + errors_ser;
        end
        if(BitErr(ii)==0)
            maxEb = ii;
            break;
        end
        BitErr(ii) = BitErr(ii)./(Escenario.tam_trama*log2(Escenario.M)*countiter);
        SerErr(ii) = SerErr(ii)./(Escenario.tam_trama*countiter);    
    end
    nombre=strcat(Escenario.tx,Escenario.Codec,num2str(Escenario.N),'x',num2str(Escenario.P),Escenario.channel,'d_BS',num2str(Escenario.d_BS),'d_MS',num2str(Escenario.d_MS),num2str(TX.Modulador.M),'.mat');
    save(nombre,'BitErr','SerErr','snr_vec','PL','BF','SNR_Threshold');
    kk=kk-1;
end
toc;