%This function extracts only the transmitter taken a previously calculated
%Rtx, and executes the Power Loading Algorithm

nTx = 4;
nRx = 1;
d_BS = 0.25;
d_MS = 0.25;

%Fixed scenario with fixed global parameters (FC centers,...), and variations on the local parameters (scatterers positions,...)
pdf = 'GAUSS';
nubes = 'LoSnoSSnoDSsi';
cd RhhDEFSistemas;
cd(pdf);
name = ['Rhh',pdf,num2str(nTx),'x',num2str(nRx),num2str(d_BS),num2str(d_MS),nubes,'.mat'];
load(name);
cd ..; cd ..;

Escenario = Rhh.Escenario{1};       %{1} Pues este es el caso 'SU'
LoS = Rhh.LoS{1};

%OTROS PARÁMETROS DE CONFIGURACIÓN
num_MPC = Rhh.num_MPC{1};
long_ch = Rhh.num_MPC{1};  %OJO!!! que en este caso está configurado EL CANAL para que coincidan, pero no siempre coincidirán!!! si los retardos
                        %No son cada Tsym, num_MPC indica esa enumeración primer, segundo, tercero, etc MPC, pero long_ch
                        %indica el máximo retardo encontrado en múltiplos de Tsym!!!
Eb_N0_dB = 0:5:40;
QAM = ['2 '; '4 '; '16'; '64' ];
kk=4;
while(kk) %For each different QAM
    filename = ['CrispreRakePower' num2str(nTx) 'x' num2str(nRx) num2str(d_BS) num2str(d_MS) nubes '-' strtrim(QAM(kk,:)) '.mat'];
    if(kk==1)
        qam = 2;
        k_QAM = 1;
        g=1; 
    elseif(kk==2)
        qam = 4;
        k_QAM = 1;
        g=3/(2*(qam-1));
    elseif(kk==3)
        qam = 16;
        k_QAM = 1;
        g=3/(2*(qam-1));
    elseif(kk==4)
        qam = 64;
        k_QAM = 1;
        g=3/(2*(qam-1));
    end
    for ii = 1:length(Eb_N0_dB)
        Dhp = Rtx.EigValues{1};                                %1 value for each tap for one block
        [Dhp,index_Dhp] = sort(Dhp,'descend');                 %To make sure they are in descending order
        Ntp = 0;
        for rr=1:num_MPC
            SNR_Threshold(rr) = (1/g)*(rr/Dhp(rr)-sum(1./Dhp(1:rr)));
            if 10^(Eb_N0_dB(ii)/10) <= SNR_Threshold(rr)
                Ntp = rr-1;
                break;
            end
            if rr==num_MPC
                Ntp = num_MPC;
            end
        end
        f_2 = zeros(num_MPC,1);
        Dh_temp = (1/Ntp)*sum(1./Dhp(1:Ntp));
        for rr=1:Ntp
            f_2(rr) = (1/Ntp)+(1/(g*10^(Eb_N0_dB(ii)/10)))*(Dh_temp-1/Dhp(rr)); 
            if f_2(rr) < 0
                f_2(rr) = 0;
            end
        end 
        Power.f_2_output(ii,:) = f_2;
        Power.SNR_Threshold_output = SNR_Threshold;
        Power.Ntp_output(ii) = Ntp;
        Power.Eb_N0_dB = Eb_N0_dB;
    end
    cd GraficaPower;
    save(filename,'Power');
    cd ..;
    kk=kk-1;
end