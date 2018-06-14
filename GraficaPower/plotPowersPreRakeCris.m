%plot of powers for our proposal PreRakeCris
close all; clc; clear all;

%PRIMERA GRÁFICA: Comparando entre diferentes separaciones de elementos de
%antena, para un sistema fijo
nTx = 4;
nRx = 1;
d_vec = [0.25 0.5];
num_MPC = 4;
QAM = ['2 '; '4 '; '16'; '64' ];
type = ['k- '; 'k--'];
nubes = 'LoSnoSSnoDSsi';
kk=4;
for qam=1:size(QAM,1)
    for nn=1:length(nTx)
        figure;
        cont = 0;
        H = [];
        for pp=1:length(nRx)
            for d = 1:length(d_vec)
                name = ['CrispreRakePower',num2str(nTx(nn)),'x',num2str(nRx(pp)),num2str(d_vec(d)),num2str(d_vec(d)),nubes,'-',strtrim(QAM(qam,:)),'.mat'];
                load(name);
                Eb_N0_dB = Power.Eb_N0_dB;
                for ii=1:length(Eb_N0_dB)
                    Distr_Power = Power.f_2_output(ii,:);
                    SNR = 10^(Eb_N0_dB(ii)/10); % SNR in linear scale
                    variance = 1/(2*SNR); % Variance
                    Distr_Power_dB(ii,:) = 10*log10(Distr_Power/variance);
                end
                SNR_Threshold{qam,d} = 10*log10(Power.SNR_Threshold_output);
                cont = cont+1;
                for ntx=1:num_MPC
                    h = plot(Eb_N0_dB,Distr_Power_dB(:,ntx),strtrim(type(cont,:)),'LineWidth',1);
                    hold on; grid on;
                end
                H = [H h];
            end
        end
        xlabel('Total Power/Noise Power (dB)');
        ylabel('Distributed Power/Noise Power (dB)');
        legend(H,'0.25\lambda','0.5\lambda','location','northwest');
    end
end

%SEGUNDA GRÁFICA:  comparando entre diferentes sistemas con separación
%entre elementos de antena fijos
nTx = [2 3 4];
nRx = [1 2 4];
d_vec = 0.5;
QAM = ['2 '; '4 '; '16'; '64' ];
type = ['ko-'; 'ks-'; 'k^-'];
nubes = 'LoSnoSSnoDSsi';
kk=4;
for qam=1:size(QAM,1)
    figure;
    cont = 0;
    H = [];
    for nn=1:length(nTx)
        cont = cont+1;
        for pp=1:length(nRx)
            for d = 1:length(d_vec)
                name = ['CrispreRakePower',num2str(nTx(nn)),'x',num2str(nRx(pp)),num2str(d_vec(d)),num2str(d_vec(d)),nubes,'-',strtrim(QAM(qam,:)),'.mat'];
                load(name);
                Eb_N0_dB = Power.Eb_N0_dB;
                for ii=1:length(Eb_N0_dB)
                    Distr_Power = Power.f_2_output(ii,:);
                    SNR = 10^(Eb_N0_dB(ii)/10); % SNR in linear scale
                    variance = 1/(2*SNR); % Variance
                    Distr_Power_dB(ii,:) = 10*log10(Distr_Power/variance);
                end
                SNR_ThresholdII{qam,nn,pp} = 10*log10(Power.SNR_Threshold_output);
                for ntx=1:num_MPC
                    h = plot(Eb_N0_dB,Distr_Power_dB(:,ntx),strtrim(type(cont,:)),'LineWidth',1);
                    hold on; grid on;
                end
            end
        end
        H = [H h];
    end
    xlabel('Total Power/Noise Power (dB)');
    ylabel('Distributed Power/Noise Power (dB)');
    legend(H,'0.5\lambda 2','0.5\lambda 3','0.5 \lambda 4','location','northwest');
end
