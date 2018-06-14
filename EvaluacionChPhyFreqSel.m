%PARA EVALUACIÓN Y VALIDACIÓN DEL CANAL FÍSICO UTILIZADO PARA EL SISTEMA COMPLETO

close all; clear all; clc;

global Escenario snr
tic;
[Escenario]=Plantilla();                                                    %Configuración escenario, arreglos de antena y canal en general
[LoS]=LineofSight();                                                        %Cálculo de rayos que van por via directa, se hace en este punto pues
                                                                            %solo depende de Escenario, por tanto debe calcularse 1 sola vez
snr_vec = 15;
MAXITER = 1000;                                                              %Number of iterrations
snr = snr_vec;

for countiter=1:MAXITER
    [H_Total(countiter),h_Total(countiter)]=CalculoH(LoS);                                                      % Generates Scatterers, and Calculates Channel Matrix
%     if countiter == 1
%        GraphicsPowerSpectrums; 
%     end
end
nombre=strcat(num2str(Escenario.N),'x',num2str(Escenario.P),Escenario.channel.ID,num2str(Escenario.d_BS),num2str(Escenario.d_MS),...
    Escenario.id,Escenario.channel.Type,'scatt',num2str(Escenario.MS_cluster.num_scattMS),'r',num2str(Escenario.MS_cluster.radio_clusterMS),...
    'LoS',Escenario.LoS,'SS',Escenario.Single_Scatt,'DS',Escenario.Double_Scatt,'.mat');
save(nombre,'Escenario','LoS','H_Total','h_Total');
%OtherGraphics;
toc;                    