%This function allows to plot the Angular (Elevation and Azimuthal) and
%Delay Power Spectrums for an instantaneos Channel Impulse Response.
H = H_Total;
h = h_Total;
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
    for ms=1:Escenario.N_d{bs}
        %PLOT OF ANGLES vs DELAY
        figure;
        subplot(2,2,1);
        h_AoDazim = polar(h.DDCIR{bs,ms}(3,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
        title(['Az AoD vs Delay BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
        subplot(2,2,2);        
        h_AoDelev = polar(h.DDCIR{bs,ms}(4,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
        title(['El AoD vs Delay BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
        subplot(2,2,3);
        h_AoAazim = polar(h.DDCIR{bs,ms}(5,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
        title(['Az AoA vs Delay BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
        subplot(2,2,4);
        h_AoAelev = polar(h.DDCIR{bs,ms}(6,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
        title(['El AoA vs Delay BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
        %PLOT OF ANGLES vs POWER
        figure;
        subplot(2,2,1);
        h_AoDazim = polar(h.DDCIR{bs,ms}(3,:)*pi/180,10*log10(h.DDCIR{bs,ms}(7,:)),'dk');
        title(['Az AoD vs Power BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
        subplot(2,2,2);        
        h_AoDelev = polar(h.DDCIR{bs,ms}(4,:)*pi/180,10*log10(h.DDCIR{bs,ms}(7,:)),'dk');
        title(['El AoD vs Power BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
        subplot(2,2,3);
        h_AoAazim = polar(h.DDCIR{bs,ms}(5,:)*pi/180,10*log10(h.DDCIR{bs,ms}(7,:)),'dk');
        title(['Az AoA vs Power BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
        subplot(2,2,4);
        h_AoAelev = polar(h.DDCIR{bs,ms}(6,:)*pi/180,10*log10(h.DDCIR{bs,ms}(7,:)),'dk');
        title(['El AoA vs Power BS',num2str(bs),'MS',num2str(Escenario.deseados{bs}(ms))])
    end
    end
end
%PLOT PDP
figure;
P_lt_norm = zeros(1,H.Limite_Taps);
P = zeros(Escenario.num_MS,H.Limite_Taps);
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
    for ms=1:Escenario.N_d{bs}
        for lt=1:H.Limite_Taps
            P_lt_norm(lt) = 10*log10((norm(H.Tap_norm{bs,ms}(:,(lt-1)*Escenario.N+1:lt*Escenario.N),'fro'))^2);
        end
        P_lt_norm(find(P_lt_norm==-Inf)) = NaN;
        P(Escenario.deseados{bs}(ms),:) = P_lt_norm;
    end
    end
end
hPDP = bar(P',.5,'group');
% %Para SU
% legend(hPDP,['MS',num2str(1)]);
%Para MU con 7 usuarios:
legend(hPDP,['MS',num2str(1)],['MS',num2str(2)],['MS',num2str(3)],['MS',num2str(4)],['MS',num2str(5)],['MS',num2str(6)],['MS',num2str(7)]);
title('Power Delay Profile');

% %PLOT OF H(tao) vs ANGLES
% for bs=1:Escenario.num_BS
%     if ~isempty(Escenario.deseados{bs})
%     for ms=1:Escenario.N_d{bs}
%         for lt=1:H.Limite_Taps
%             ind_tap_lt = find(h.DDCIR{bs,ms}(1,:)>=((lt-1)*H.Tsymbol) & h.DDCIR{bs,ms}(1,:)<(lt*H.Tsymbol));
%             amp = h.DDCIR{bs,ms}(8,ind_tap_lt).*exp(j*h.DDCIR{bs,ms}(2,ind_tap_lt));
%             k_BS = 2*pi*[cos(h.DDCIR{bs,ms}(3,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(4,ind_tap_lt)*pi/180); ...
%                       sin(h.DDCIR{bs,ms}(3,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(4,ind_tap_lt)*pi/180); ...
%                       sin(h.DDCIR{bs,ms}(4,ind_tap_lt)*pi/180)].';
%             AMV_BS = kron(ones(Escenario.N,1),amp).*exp(-j*k_BS*Escenario.BS_array)';
%             k_MS = 2*pi*[cos(h.DDCIR{bs,ms}(5,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(6,ind_tap_lt)*pi/180); ...
%                       sin(h.DDCIR{bs,ms}(5,ind_tap_lt)*pi/180).*cos(h.DDCIR{bs,ms}(6,ind_tap_lt)*pi/180); ...
%                       sin(h.DDCIR{bs,ms}(6,ind_tap_lt)*pi/180)].';
%             AMV_MS = exp(-j*k_MS*Escenario.MS_array)';
%             AMV_BS_temp = kron(ones(Escenario.P,1),reshape(AMV_BS,1,[]));
%             AMV_MS_temp = kron(AMV_MS,ones(1,Escenario.N));
%             ch_temp = AMV_BS_temp.*AMV_MS_temp;
%             for ii=1:length(ind_tap_lt)
%                 P_lt(ii) = (norm(ch_temp((ii-1)*Escenario.N+1:ii*Escenario.N),'fro'))^2;
%             end
%             P_lt_dB = 10*log10(P_lt);
%             figure;
%             [X,Y] = meshgrid(h.DDCIR{bs,ms}(3,ind_tap_lt),h.DDCIR{bs,ms}(5,ind_tap_lt));
%             plot3(X,Y,P_lt_dB);
%         end
%     end
%     end
% end
% 
% 
