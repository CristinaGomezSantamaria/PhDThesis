%plot of powers for our proposal PreRakeCris
close all; clc; clear all;

nTx = 4;
nRx = 1;
d_vec = [0.25 0.5];
QAM = ['2 '; '4 '; '16'; '64' ];
type = ['k.-'; 'k: '; 'k.-'; 'k:.'];
nubes = 'LoSnoSSnoDSsi';
kk=4;
while(kk)
    figure;
    H = [];
    cd GraficaPower;
    for d = 1:length(d_vec)
        name = ['CrispreRakePower',num2str(nTx),'x',num2str(nRx),num2str(d_vec(d)),num2str(d_vec(d)),nubes,'-',strtrim(QAM(kk,:)),'.mat'];
        load(name);
        Eb_N0_dB = Power.Eb_N0_dB;
        for ii=1:length(Eb_N0_dB)
            Distr_Power = Power.f_2_output(ii,:);
            SNR = 10^(Eb_N0_dB(ii)/10); % SNR in linear scale
            variance = 1/(2*SNR); % Variance
            Distr_Power_dB(ii,:) = 10*log10(Distr_Power/variance);
            %Distr_Power_dB(ii,find(Distr_Power_dB(ii,:)==-Inf)) = -5;
        end
        for ntx=1:nTx
            h = plot(Eb_N0_dB,Distr_Power_dB(:,ntx),strtrim(type(d,:)),'LineWidth',1);
            hold on;
        end
        H = [H h];
    end
    xlabel('Total Power/Noise Power (dB)');
    ylabel('Distributed Power/Noise Power (dB)');
    legend(H,'0.25 \lambda','0.5 \lambda');
    kk=kk-1;
    cd ..;
end
