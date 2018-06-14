%Calculation of Rhh and RTx for Final Scenarios INCLUDING ALL THE
%EIGENVALUES DIFFERENT FROM ZERO, NOT ONLY THE HIGHEST IN EACH TAP
%OPTIONS:
%- with LoS and pure DS (num_scatt=12, radio=100)
%- without LoS and pure DS (num_scatt=12, radio=100)
%- with LoS and DS&SS (num_scatt=50, radio=100)
%- without LoS and DS&SS (num_scatt=50, radio=100)
%CONFIGURE INPUT AND OUTPUT DIRECTORY AND FILE NAMES ACCORDING TO THE
%SCENARIO TO CALCULATE, AND NUMBER OF DISPERSORS AND RADIO

clear all; clc; close all;

MAXITER = 1000;
N_vec = 4; %[2 3];
P_vec = 2; %[1 2 4];
d_vec = 0.5; %[0.25 0.5];
r = 100;
num_disp = 30;
nube = 'LoSsiSSsiDSsi';
for nn=1:length(N_vec)
    for pp=1:length(P_vec)
        for d=1:length(d_vec)
            cd RhhDEFSistemas\'Gauss'\'StrongestEigenvalue'\VariandoNumDisp;     
            name2 = strcat(num2str(N_vec(nn)),'x',num2str(P_vec(pp)),'Physical',num2str(d_vec(d)),num2str(d_vec(d)),...
                'SUMacrocellBUscatt',num2str(num_disp),'r',num2str(r),nube,'.mat');  
                              %it can be LoSsiSSnoDSsi, LoSnoSSnoDSsi, LoSsiSSsiDSsi, LoSnoSSsiDSsi
            load(name2);
            Rtx_temp = cell(Escenario.num_MS,H_Total(1).Limite_Taps);
            Rhh_temp = cell(Escenario.num_MS,H_Total(1).Limite_Taps);
            Rhh.EigVectors = cell(Escenario.num_MS);
            Rtx.EigVectors = cell(Escenario.num_MS);
            for bs=1:Escenario.num_BS
                for ms=1:Escenario.N_d{bs}
                    for lt=1:H_Total(1).Limite_Taps
                        Rhh_temp{Escenario.deseados{bs}(ms),lt} = zeros(Escenario.N*Escenario.P,Escenario.N*Escenario.P);
                        Rtx_temp{Escenario.deseados{bs}(ms),lt} = zeros(Escenario.N,Escenario.N);
                        for countiter=1:MAXITER
                            h = reshape(H_Total(countiter).Tap_norm{bs,ms}(:,(lt-1)*Escenario.N+1:lt*Escenario.N).',[],1);
                            htx = H_Total(countiter).AMVBS_norm{Escenario.deseados{bs}(ms),lt};
                            Rhh_temp{Escenario.deseados{bs}(ms),lt} = Rhh_temp{Escenario.deseados{bs}(ms),lt}+h*h';
                            Rtx_temp{Escenario.deseados{bs}(ms),lt} = Rtx_temp{Escenario.deseados{bs}(ms),lt}+htx*htx';
                        end
                        Rhh_temp{Escenario.deseados{bs}(ms),lt} = Rhh_temp{Escenario.deseados{bs}(ms),lt}/MAXITER;
                        Rtx_temp{Escenario.deseados{bs}(ms),lt} = Rtx_temp{Escenario.deseados{bs}(ms),lt}/MAXITER;
                        [u,s,v] = svd(Rhh_temp{Escenario.deseados{bs}(ms),lt});
                        Rhh.EigValues{Escenario.deseados{bs}(ms)}((lt-1)*N_vec(nn)*P_vec(pp)+1:lt*N_vec(nn)*P_vec(pp)) = diag(s)';
                        Rhh.EigVectors{Escenario.deseados{bs}(ms)}(:,(lt-1)*N_vec(nn)*P_vec(pp)+1:lt*N_vec(nn)*P_vec(pp)) = u;
                        Rhh.Escenario{Escenario.deseados{bs}(ms)} = Escenario;
                        Rhh.LoS{Escenario.deseados{bs}(ms)} = LoS;
                        Rhh.num_MPC{Escenario.deseados{bs}(ms)} = H_Total(1).Limite_Taps;
                        [utx,stx,vtx] = svd(Rtx_temp{Escenario.deseados{bs}(ms),lt});
                        Rtx.EigValues{Escenario.deseados{bs}(ms)}((lt-1)*N_vec(nn)+1:lt*N_vec(nn)) = diag(stx)';
                        Rtx.EigVectors{Escenario.deseados{bs}(ms)}(:,(lt-1)*N_vec(nn)+1:lt*N_vec(nn)) = utx;
                        Rtx.Rtx{Escenario.deseados{bs}(ms),lt} = Rtx_temp{Escenario.deseados{bs}(ms),lt};
                     end
                end
            end            
            name = ['RhhGauss',strcat(num2str(N_vec(nn))),'x',num2str(P_vec(pp)),num2str(d_vec(d)),num2str(d_vec(d)),...
                'num_disp',num2str(num_disp),nube,'.mat'];
            save(name,'Rhh','Rtx');
        end
    end
end