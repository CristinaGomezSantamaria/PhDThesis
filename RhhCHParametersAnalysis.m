%Plots and analysis of the results for different channels varying
%num_scatt, radious and pdf
clear all; clc; close all;

MAXITER = 1000;
N_vec = 4;
% P_vec = [1 2 4];
% d_vec = [0.25 0.5];
% scatt_vec = [12 50 100 200]';

%Para el único caso de NLoS
P_vec = 4;
d_vec = 0.5;
scatt_vec = [4 6 12 50 100 200]';
radio_vec = [100 200];
tipo = ['-.'; '--'];
cont = 0;
for nn=1:length(N_vec)
    for pp=1:length(P_vec)
        s_graph = cell(1,4,2);      %OJO CON ESTA INICIALIZACIÓN QUE TOCA DE PARTIDA ASUMIR CONOCIDOS PARÁMETROS PUES ES: 
                                    %CELL{Escenario.num_MS,H.Limite_Taps,length(radio_vec), pero como aún no se conocen ni 
                                    %Escenario ni H porque no se han cargado, toca introducirlos desde aquí
                                    %manualmente conociéndolos previamente
        for d=1:length(d_vec)
            for rd=1:length(radio_vec)
                for sc=1:length(scatt_vec)
                    cd AnalisisNumScattyRadio\'GAUSSCYL';      %PDF: GaussCyl suggested by Molisch, or Gauss by COST259
                    name1 = strcat(num2str(N_vec(nn)),'x',num2str(P_vec(pp)),'NLoSIIdBS',num2str(d_vec(d)),'dMS',num2str(d_vec(d)));
                    cd(name1);        %IF 'both' IS SELECTED THE CHANNEL TAKES INTO ACCOUNT DS & MS, IF 'DS' THEN ONLY DS IS CONSIDERED                
                    name2 = strcat(num2str(N_vec(nn)),'x',num2str(P_vec(pp)),'Physical',num2str(d_vec(d)),num2str(d_vec(d)),...
                        'SUMacrocellBUscatt',num2str(scatt_vec(sc)),'r',num2str(radio_vec(rd)),'SSnoDSsi.mat');  
                                      %IF 'both' WAS SELECTED PREVIOUSLY, HERE IT MUST BE 'SSsiDSsi', IF 'DS' WAS SELECTED PREVIOUSLY THEN HERE IT
                                      %MUST BE 'SSnoDSsi'
                    load(name2);
                    cd ..; cd ..;
                    Rhh_temp = cell(Escenario.num_MS,H_Total(1).Limite_Taps);
                    cont = cont+1;
                    for bs=1:Escenario.num_BS
                        for ms=1:Escenario.N_d{bs}
                            for lt=1:H_Total(1).Limite_Taps
                                Rhh_temp{Escenario.deseados{bs}(ms),lt} = zeros(Escenario.N*Escenario.P,Escenario.N*Escenario.P);
                                for countiter=1:MAXITER
                                    h = reshape(H_Total(countiter).Tap_norm{bs,ms}(:,(lt-1)*Escenario.N+1:lt*Escenario.N).',[],1);
                                    Rhh_temp{Escenario.deseados{bs}(ms),lt} = Rhh_temp{Escenario.deseados{bs}(ms),lt}+h*h';
                                end
                                Rhh_temp{Escenario.deseados{bs}(ms),lt} = Rhh_temp{Escenario.deseados{bs}(ms),lt}/MAXITER;
                                [u,s,v] = svd(Rhh_temp{Escenario.deseados{bs}(ms),lt});
                                s_graph{Escenario.deseados{bs}(ms),lt,rd}(sc,:) = sum(s);
                                Rhh.EigValues{Escenario.deseados{bs}(ms),cont}(lt) = s(1,1);
                                Rhh.EigVectors{Escenario.deseados{bs}(ms),cont}(:,lt) = u(:,1);
                                Rhh.Escenario{Escenario.deseados{bs}(ms),cont} = Escenario;
                                Rhh.LoS{Escenario.deseados{bs}(ms),cont} = LoS;
                                Rhh.num_MPC{Escenario.deseados{bs}(ms),cont} = H_Total(1).Limite_Taps;
                            end
                        end
                    end
                end
            end
            for bs=1:Escenario.num_BS
                for ms=1:Escenario.N_d{bs}
                    figure; H = [];
                    for lt=1:H_Total(1).Limite_Taps
                        subplot(2,2,lt);
                        for rd=1:length(radio_vec)
                                h = plot(kron(scatt_vec,ones(1,size(s_graph{Escenario.deseados{bs}(ms),lt,rd},2))),...
                                    s_graph{Escenario.deseados{bs}(ms),lt,rd},'LineStyle',tipo(rd,:)); 
                                hold on; grid on;
                                title(['System',num2str(N_vec(nn)),'x',num2str(P_vec(pp)),'dBS',num2str(d_vec(d)),'dMS',num2str(d_vec(d)),...
                                    'tap',num2str(lt)]);
                                H = [H; h(1:P_vec(pp))];
                        end
                        %axis([0 200 0 max(Rhh.EigValues{Escenario.deseados{bs}(ms),cont})+0.1])
                        switch pp
                            case 1
                                legend(H,'\sigma_1 r 100','\sigma_1 r 200');
                            case 2
                                legend(H,'\sigma_1 r 100','\sigma_2 r 100','\sigma_1 r 200','\sigma_2 r 200');
                                        %'\sigma_1 r 100','\sigma_2 r 100','\sigma_3 r 100','\sigma_4 r 100',...
                                        %'\sigma_5 r 100','\sigma_6 r 100','\sigma_7 r 100','\sigma_8 r 100',...
                                        %'\sigma_1 r 200','\sigma_2 r 200','\sigma_3 r 200','\sigma_4 r 200',...
                                        %'\sigma_5 r 200','\sigma_6 r 200','\sigma_7 r 200','\sigma_8 r 200');
                            case 3
                                legend(H,'\sigma_1 r 100','\sigma_2 r 100','\sigma_3 r 100','\sigma_4 r 100',...
                                         '\sigma_1 r 200','\sigma_2 r 200','\sigma_3 r 200','\sigma_4 r 200');
%                                          '\sigma_1 r 100','\sigma_2 r 100','\sigma_3 r 100','\sigma_4 r 100',...
%                                          '\sigma_5 r 100','\sigma_6 r 100','\sigma_7 r 100','\sigma_8 r 100',...
%                                          '\sigma_9 r 100','\sigma_10 r 100','\sigma_11 r 100','\sigma_12 r 100',...
%                                          '\sigma_12 r 100','\sigma_14 r 100','\sigma_15 r 100','\sigma_16 r 100',...
%                                          '\sigma_1 r 200','\sigma_2 r 200','\sigma_3 r 200','\sigma_4 r 200',...
%                                          '\sigma_5 r 200','\sigma_6 r 200','\sigma_7 r 200','\sigma_8 r 200',...
%                                          '\sigma_9 r 200','\sigma_10 r 200','\sigma_11 r 200','\sigma_12 r 200',...
%                                          '\sigma_13 r 200','\sigma_14 r 200','\sigma_15 r 200','\sigma_16 r 200');
                        end
                    end
                end
            end
        end
    end
end

name = ['RhhGaussCylNLoSII','.mat'];
save(name,'Rhh');