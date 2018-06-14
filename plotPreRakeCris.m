%plot results for our proposal PreRakeCris
close all; clear all; clc;

nTx = 4;
nRx = 2;
d_BS = 0.5;
d_MS = 0.5;
nubes = 'LoSnoSSnoDSsi';
QAM = ['Ber  '; 'Ber4 '; 'Ber16'; 'Ber64'];
QAM2 = ['Ser  '; 'Ser4 '; 'Ser16'; 'Ser64'];
H_ber = []; H_ser = [];
fig_ber = figure;
fig_ser = figure;
for qam = 1:size(QAM,1)
    filename = ['CrispreRakeOSTBCthreequarterPhy' num2str(nTx) 'x' num2str(nRx) num2str(d_BS) num2str(d_MS) nubes '-' strtrim(QAM) '.mat'];
    load(filename);
    figure(fig_ber);
    switch QAM
        case 'Ber64'
            h_ber = semilogy(Eb_N0_dB,Ber64,'-:');
        case 'Ber16'
            h_ber = semilogy(Eb_N0_dB,Ber16,'--');
        case 'Ber4 '
            h_ber = semilogy(Eb_N0_dB,Ber4,'-.');
        case 'Ber  '
            h_ber = semilogy(Eb_N0_dB,Ber,'..');
    end
    H_ber = [H_ber h_ber];
    hold on; grid on;
    figure(fig_ser);
    switch QAM
        case 'Ber64'
            h_ser = semilogy(Eb_N0_dB,Ser64,'-:');
        case 'Ber16'
            h_ser = semilogy(Eb_N0_dB,Ser16,'--');
        case 'Ber4 '
            h_ser = semilogy(Eb_N0_dB,Ser4,'-.');
        case 'Ber  '
            h_ser = semilogy(Eb_N0_dB,Ser,'..');
    end
    H_ser = [H_ser h_ser];
    hold on; grid on;
end
figure(fig_ber);
title('System with BF+PL+OSTBC for typical Macrocell Channels')
xlabel('SNR(dB)')
ylabel('BER')
legend(H_ber,'64QAM','16QAM','QPSK','BPSK');

figure(fig_ser);
title('System with BF+PL+OSTBC for typical Macrocell Channels')
xlabel('SNR(dB)')
ylabel('SER')
legend(H_ser,'64QAM','16QAM','QPSK','BPSK')

