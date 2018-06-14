%transmitter and receiver for an BF+PL+OSTBC with rate 3/4 with Physical
%Frequency Selective Fading Channel for several transmitted blocks

%PARA ENSAYAR!!! Rhh de GaussCyl 4x2both0.50.5scatt12r100SS.mat

function done = OSTBCthreequarterCrispreRake1DBF()

global Escenario snr

done=0;
%ESTO SE PUEDE CONFIGURAR COMO ARGUMENTOS DE ENTRADA A LA FUNCIÓN!!!
%PARÁMETROS DE ENTRADA SEGÚN SISTEMA Y PARÁMETROS GLOBALES DEL CANAL!!!
nTx = 4;
nRx = 2;
d_BS = 0.5;
d_MS = 0.5;

%Fixed scenario with fixed global parameters (FC centers,...), and variations on the local parameters (scatterers positions,...)
pdf = 'GAUSS';
nubes = 'LoSsiSSnoDSsi';
cd RhhDEFSistemas;
cd(pdf);
%CON RADIO 50
cd Radio50;
name = ['Rhh',pdf,num2str(nTx),'x',num2str(nRx),num2str(d_BS),num2str(d_MS),nubes,'.mat'];
load(name);
cd ..; cd ..;
%CON RADIO 50
cd ..;

Escenario = Rhh.Escenario{1};       %{1} Pues este es el caso 'SU'
LoS = Rhh.LoS{1};

%OTROS PARÁMETROS DE CONFIGURACIÓN
num_MPC = Rhh.num_MPC{1};
long_ch = Rhh.num_MPC{1};  %OJO!!! que en este caso está configurado EL CANAL para que coincidan, pero no siempre coincidirán!!! si los retardos
                        %No son cada Tsym, num_MPC indica esa enumeración primer, segundo, tercero, etc MPC, pero long_ch
                        %indica el máximo retardo encontrado en múltiplos de Tsym!!!

Eb_N0_dB = 0:5:40;
N = 10^4;
blockSym = 3;
length_block = 4;
THRESHOLD = 100; % Minimum number of errors before proceeding to next SNR
MAXITER=100; % Maximum number of iterrations before quiting
% The total number of transmitted symbols must be adjusted to be an
% exact multiple of blockSym
while(mod(N,blockSym))
    N=N+1;
end
N_blocks=N/blockSym;    %Number of code blocks transmitted

kk=4;
%if kk==3
while(kk) %For each different QAM
    filename = ['CrispreRakeOSTBCthreequarterPhy' num2str(nTx) 'x' num2str(nRx) num2str(d_BS) num2str(d_MS) nubes];
    filevar = ['Eb_N0_dB'; 'maxEb   '; 'Ber     '; 'Ber4    '; 'Ber16   '; 'Ber64   ';'Ser     ';'Ser4    ';'Ser16   ';'Ser64   '];

    BitErr=zeros(1,length(Eb_N0_dB)); % Array for counting biterrors
    SerErr=zeros(1,length(Eb_N0_dB)); % Array for counting sererrors
    %k_QAM normalizes transmit power to unity
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
    Modulador = modem.qammod('M',qam,'SYMBOLORDER','gray','INPUTTYPE','integer');    %Objeto Modulador
    Demodulador = modem.qamdemod('M',qam,'SYMBOLORDER','gray','OUTPUTTYPE','integer');    %Objeto Demodulador

    % Begin simulation
    for ii = 1:length(Eb_N0_dB)
        countiter = 0; % Number of iterations
        fprintf('\n At SNR = %d of %s \n',Eb_N0_dB(ii),strtrim(filevar(kk+2,:)));

        while(THRESHOLD>BitErr(ii) && countiter<MAXITER)
            countiter = countiter+1;

            %TRANSMISOR
            %Generating the random symbols to be transmitted
            transmitted = randint(N,1,[0,qam-1]);
            %Modulation
            bitmod = modulate(Modulador,transmitted);
            s_precoded = bitmod *k_QAM; % scaling the modulated signal according to its constellation
            s_precoded = conj(s_precoded');
            %Coding
            X1=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            X2=[0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0];
            X3=[0 0 0 -1; 0 0 1 0; 0 -1 0 0; 1 0 0 0];
            Y1=[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1];
            Y2=[0 0 1 0; 0 0 0 -1; 1 0 0 0; 0 -1 0 0];
            Y3=[0 0 0 -1; 0 0 -1 0; 0 -1 0 0; -1 0 0 0];
            X11=repmat(X1,1,N_blocks);
            X22=repmat(X2,1,N_blocks);
            X33=repmat(X3,1,N_blocks);
            Y11=repmat(Y1,1,N_blocks);
            Y22=repmat(Y2,1,N_blocks);
            Y33=repmat(Y3,1,N_blocks);

            S_1=kron(s_precoded(1:3:end),ones(num_MPC,length_block));
            S_2=kron(s_precoded(2:3:end),ones(num_MPC,length_block));
            S_3=kron(s_precoded(3:3:end),ones(num_MPC,length_block));
            codSignal=X11.*real(S_1)+X22.*real(S_2)+X33.*real(S_3)+j*Y11.*imag(S_1)+j*Y22.*imag(S_2)+j*Y33.*imag(S_3);

            %Channel correlations at transmitter of the channel known at the transmitter, TAKEN FROM A PREVIOUSLY SIMULATED PHYSICAL CHANNEL
            %Calculating the Beamformer
            Dhp = Rtx.EigValues{1};                                %1 value for each tap for one block
            [Dhp,index_Dhp] = sort(Dhp,'descend');                 %To make sure they are in descending order
            Uh = Rtx.EigVectors{1}(:,index_Dhp);                   %1 vector (nTx x 1) for each tap for one block (nTx x num_MPC,1), assuring
                                                                   %they follow the same order as the Eigenvalues
            delayS_matrix = index_Dhp-1;                           %OJO!!! in this case because num_MPC = long_ch, thus the index is also the
                                                                   %delay, but is not always the case
            Uh2 = zeros(nTx*num_MPC,num_MPC);                      %forming the diagonal matrix in Uh for one block (nTx x num_MPC,num_MPC)
            for ll=1:num_MPC
                Uh2((ll-1)*nTx+1:ll*nTx,ll) = Uh(:,ll);
            end
            Uh = Uh2;
            %Power Loading: because there is a full coupling between signals transmitted at different antennas with far scatterers,
            %the algorithm is reduced to calculate for num_MPC possible beams
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
            codSignal_out2 = zeros(nTx*num_MPC,N_blocks*(long_ch+length_block-1));
            codSignal_out3 = zeros(nTx*num_MPC,N_blocks*length_block+long_ch-1);
%             Dh_temp = (1/Ntp)*sum(1./Dhp(1:Ntp));
%             for rr=1:Ntp
%                 f_2(rr) = (1/Ntp)+(1/(g*10^(Eb_N0_dB(ii)/10)))*(Dh_temp-1/Dhp(rr)); 
%                 if f_2(rr) < 0
%                     f_2(rr) = 0;
%                 end
%             end
            f_2 = [1 0 0 0]';
            for nb=1:N_blocks
                codSignal_out(:,(nb-1)*length_block+1:nb*length_block) = (codSignal(:,(nb-1)*length_block+1:nb*length_block).'*diag(sqrt(f_2))*...
                    Uh').';
            end

            %Si efectivamente todos los retardos son iguales durante la simulación
            for ll=1:num_MPC
                ref_ini = long_ch-delayS_matrix(ll)-1;
                codSignal_out3((ll-1)*nTx+1:ll*nTx,ref_ini+1:ref_ini+N_blocks*length_block) = codSignal_out((ll-1)*nTx+1:ll*nTx,:);
            end
            
            %CHANNEL
            y = zeros(nRx,N_blocks*length_block+long_ch-1);
            y2 = zeros(nRx*num_MPC,N_blocks*length_block+long_ch-1);
            snr = 10^(Eb_N0_dB(ii)/10); % SNR in linear scale
            for nb=1:N_blocks
                [H_Total(nb),h_Total(nb)] = CalculoH(LoS);   
                for ll=1:num_MPC
                    ref_ini = long_ch-delayS_matrix(ll)-1;
                    ref_ini2 = ref_ini+delayS_matrix(ll);
                    y_temp = H_Total(nb).Tap_norm{1}(:,(ll-1)*nTx+1:ll*nTx)*...
                        codSignal_out3((ll-1)*nTx+1:ll*nTx,ref_ini+(nb-1)*length_block+1:ref_ini+nb*length_block);
                    y2((ll-1)*nRx+1:ll*nRx,ref_ini2+(nb-1)*length_block+1:ref_ini2+nb*length_block) = y_temp;
                end
            end
            for ll=1:num_MPC
                y = y+y2((ll-1)*nRx+1:ll*nRx,:);
            end
            y = y(:,long_ch:end);
            
            n = 1/sqrt(2)*(randn(nRx,length_block*N_blocks) + j*randn(nRx,length_block*N_blocks));
            variance = 1/(2*snr); % Variance
            sigma = sqrt(variance); % Standard Deviation
            y = y+sigma*n;   %received signal

            %RECEIVER
            y= y/k_QAM;
            %Combining
            f_2_rx = zeros(size(f_2));
            f_2_rx(find(f_2~=0)) = 1./f_2(find(f_2~=0));

            s_1 = zeros(1,N_blocks); s_2 = zeros(1,N_blocks); s_3 = zeros(1,N_blocks);
            for nb=1:N_blocks
                for ll=1:num_MPC
                    s_1R = trace(real(conj(H_Total(nb).Tap_norm{1}(:,(ll-1)*nTx+1:ll*nTx))*Uh((ll-1)*nTx+1:ll*nTx,ll)*...
                        sqrt(f_2_rx(ll))*conj(X1(ll,:))*y(:,(nb-1)*length_block+1:nb*length_block).'));
                    s_2R = trace(real(conj(H_Total(nb).Tap_norm{1}(:,(ll-1)*nTx+1:ll*nTx))*Uh((ll-1)*nTx+1:ll*nTx,ll)*...
                        sqrt(f_2_rx(ll))*conj(X2(ll,:))*y(:,(nb-1)*length_block+1:nb*length_block).'));
                    s_3R = trace(real(conj(H_Total(nb).Tap_norm{1}(:,(ll-1)*nTx+1:ll*nTx))*Uh((ll-1)*nTx+1:ll*nTx,ll)*...
                        sqrt(f_2_rx(ll))*conj(X3(ll,:))*y(:,(nb-1)*length_block+1:nb*length_block).'));


                    s_1I = trace(real(conj(H_Total(nb).Tap_norm{1}(:,(ll-1)*nTx+1:ll*nTx))*Uh((ll-1)*nTx+1:ll*nTx,ll)*...
                        sqrt(f_2_rx(ll))*conj(j*Y1(ll,:))*y(:,(nb-1)*length_block+1:nb*length_block).'));
                    s_2I = trace(real(conj(H_Total(nb).Tap_norm{1}(:,(ll-1)*nTx+1:ll*nTx))*Uh((ll-1)*nTx+1:ll*nTx,ll)*...
                        sqrt(f_2_rx(ll))*conj(j*Y2(ll,:))*y(:,(nb-1)*length_block+1:nb*length_block).'));
                    s_3I = trace(real(conj(H_Total(nb).Tap_norm{1}(:,(ll-1)*nTx+1:ll*nTx))*Uh((ll-1)*nTx+1:ll*nTx,ll)*...
                        sqrt(f_2_rx(ll))*conj(j*Y3(ll,:))*y(:,(nb-1)*length_block+1:nb*length_block).'));

                    s_1(nb) = s_1(nb)+s_1R+j*s_1I;
                    s_2(nb) = s_2(nb)+s_2R+j*s_2I;
                    s_3(nb) = s_3(nb)+s_3R+j*s_3I;
                end
            end

            %Cálculo de métricas
            X1_dec = zeros(qam,size(s_1,1));            
            X2_dec = zeros(qam,size(s_2,1));
            X3_dec = zeros(qam,size(s_3,1));
            decoded1 = zeros(1,N_blocks); decoded2 = zeros(1,N_blocks); decoded3 = zeros(1,N_blocks);
            % decoded = zeros(3,N_blocks);
            for nb=1:N_blocks
                for k=1:qam
                     %Cómputo de distancias a cada símbolo

                     d1(k,:) = (abs(s_1(nb)-Modulador.Constellation(k))).^2;
                     d2(k,:) = (abs(s_2(nb)-Modulador.Constellation(k))).^2;
                     d3(k,:) = (abs(s_3(nb)-Modulador.Constellation(k))).^2;

                    d12=(real(Modulador.Constellation(k))-real(s_1(nb))).^2+(imag(Modulador.Constellation(k))-imag(s_1(nb))).^2;
                    d22=(real(Modulador.Constellation(k))-real(s_2(nb))).^2+(imag(Modulador.Constellation(k))-imag(s_2(nb))).^2;
                    d32=(real(Modulador.Constellation(k))-real(s_3(nb))).^2+(imag(Modulador.Constellation(k))-imag(s_3(nb))).^2;

                    %Construcción de vector de decisiones
                    X1_dec(k,:)= (abs(Modulador.Constellation(k))^2)*(sum(sum((abs(H_Total(nb).Tap_norm{1})).^2))-1)+d12;    
                    X2_dec(k,:)= (abs(Modulador.Constellation(k))^2)*(sum(sum((abs(H_Total(nb).Tap_norm{1})).^2))-1)+d22;
                    X3_dec(k,:)= (abs(Modulador.Constellation(k))^2)*(sum(sum((abs(H_Total(nb).Tap_norm{1})).^2))-1)+d32;    

                end

                %Decisiones: las que tengan la mínima métrica
                [min_metric1,index1]=min(X1_dec);
                [min_metric2,index2]=min(X2_dec);
                [min_metric3,index3]=min(X3_dec);

                decoded1(nb) = Modulador.Constellation(index1); 
                decoded2(nb) = Modulador.Constellation(index2);
                decoded3(nb) = Modulador.Constellation(index3);    
                decoded(:,nb) = [decoded1(nb); decoded2(nb); decoded3(nb)];
                decodSignal = reshape(decoded,1,size(decoded,2)*3);
            end
            
            %Demodulator
            received = zeros(N,1);
            received(1:3:end) = demodulate(Demodulador,decoded1);           
            received(2:3:end) = demodulate(Demodulador,decoded2);
            received(3:3:end) = demodulate(Demodulador,decoded3);
            
            %Counting Bit Error
            errors = sum(sum(abs(dec2bin(transmitted)-dec2bin(received))));
            BitErr(ii) = BitErr(ii) + errors;
            errors_ser = length(find(transmitted~=received));
            SerErr(ii) = SerErr(ii) + errors_ser;
        end
        if(BitErr(ii)==0)
            maxEb = ii;
            break;
        else
            maxEb = Eb_N0_dB(end);
        end
        BitErr(ii) = BitErr(ii)./(N*log2(qam)*countiter);
        SerErr(ii) = SerErr(ii)./(N*countiter);        
    end
    if(kk==1)
        Ber = BitErr;
        Ser = SerErr;
        tempname = [filename '-' strtrim(filevar(kk+2,:)) '.mat'];
        save(tempname, strtrim(filevar(1,:)),...
            strtrim(filevar(2,:)), strtrim(filevar(kk+2,:)), strtrim(filevar(kk+6,:)));
    elseif(kk==2)
        Ber4 = BitErr;
        Ser4 = SerErr;
        tempname = [filename '-' strtrim(filevar(kk+2,:)) '.mat'];
        save(tempname, strtrim(filevar(1,:)),...
            strtrim(filevar(2,:)), strtrim(filevar(kk+2,:)), strtrim(filevar(kk+6,:)));
    elseif(kk==3)
        Ber16 = BitErr;
        Ser16 = SerErr;
        tempname = [filename '-' strtrim(filevar(kk+2,:)) '.mat'];        
        save(tempname, strtrim(filevar(1,:)),...
            strtrim(filevar(2,:)), strtrim(filevar(kk+2,:)), strtrim(filevar(kk+6,:)));
    else
        Ber64 = BitErr;
        Ser64 = SerErr;
        tempname = [filename '-' strtrim(filevar(kk+2,:)) '.mat'];
        save(tempname, strtrim(filevar(1,:)),...
            strtrim(filevar(2,:)), strtrim(filevar(kk+2,:)), strtrim(filevar(kk+6,:)));
    end
    kk=kk-1;
end
fprintf('\n');
done = 1;
% DONE