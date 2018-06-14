%This function calculates the CDF and PDF of the PDP for the systems
%simulated and compares between them, is useful for the analysis of
%Bandwith Coherence...

clear all; close all; clc;

%OJO!!! CONFIGURAR YA SEA PARA GRAFICAR PARA VARIOS SISTEMAS (NXP) CON NUBE
%FIJA O PARA GRAFICAR VARIAS NUBES PARA UN SISTEMA FIJO!!!
nTx = 4; %[2 3 4];
nRx = 2; %[1 2 4];
d_vec = [0.5];
nubes = ['LoSnoSSnoDSsi'; 'LoSsiSSnoDSsi'; 'LoSnoSSsiDSsi'; 'LoSsiSSsiDSsi']; %'LoSnoSSnoDSsi';

MAXITER = 1000;
snr_vec = 15;
cont = 0;
for nn=1:length(nTx)
    for pp=1:length(nRx)
        for d=1:length(d_vec)
            for nb=1:size(nubes,1)
                if strcmp(nubes(nb,:),'LoSnoSSnoDSsi') || strcmp(nubes(nb,:),'LoSsiSSnoDSsi')
                    num_scatt = 12;
                elseif strcmp(nubes(nb,:),'LoSnoSSsiDSsi') || strcmp(nubes(nb,:),'LoSsiSSsiDSsi')
                    num_scatt = 50;
                end
                cd 'RhhDEFSistemas'/'Gauss';
                %PARA RADIO 50 O 50FC
                cd 'Radio50';
                nombre=strcat(num2str(nTx(nn)),'x',num2str(nRx(pp)),'Physical',num2str(d_vec(d)),num2str(d_vec(d)),'SUMacrocellBU','scatt',...
                    num2str(num_scatt),'r',num2str(50),nubes(nb,:),'.mat');
                load(nombre);
                cd ..; cd ..;
                %PARA RADIO 50 O 50FC
                cd ..;

                %DELAY SPREAD: GENEAL FOR EACH LINK
                cont = cont+1;
                x_delay = 0:H_Total(1).Tsymbol/10:H_Total(1).Max_Retardo;
                Delay_temp = zeros(MAXITER,Escenario.num_MS);
                for countiter=1:MAXITER
                    Delay_temp(countiter,:) = H_Total(countiter).D_S_norm2;
                end
                pdf_delay(cont,:) = hist(Delay_temp,x_delay);
                cdf_delay(cont,:) = cumsum(pdf_delay(cont,:))./sum(pdf_delay(cont,:));
            end
        end
    end
end

h1 = bar(x_delay,cdf_delay',1,'grouped'); grid on;
xlabel('Delay Spread (seg)'); ylabel('Prob A_S<Abscissa');
legend(h1,'DS','LoSDS','SSDS','LoSSSDS','Location','NorthWest');

% legend(h1,'2x1','2x2','2x4',...
%           '3x1','3x2','3x4',...
%           '4x1','4x2','4x4','Location','NorthWest');

figure;
h2 = plot(x_delay,cdf_delay); grid on;
xlabel('Delay Spread (seg)'); ylabel('Prob A_S<Abscissa');
legend(h2,'DS','LoSDS','SSDS','LoSSSDS','Location','NorthWest');
% legend(h2,'2x1','2x2','2x4',...
%           '3x1','3x2','3x4',...
%           '4x1','4x2','4x4','Location','NorthWest');

figure;
h3 = bar(x_delay,pdf_delay'/MAXITER,0.5,'grouped'); grid on;
xlabel('Delay Spread (seg)'); ylabel('Prob A_S<Abscissa');
legend(h3,'DS','LoSDS','SSDS','LoSSSDS','Location','NorthWest');
% legend(h3,'2x1','2x2','2x4',...
%           '3x1','3x2','3x4',...
%           '4x1','4x2','4x4');