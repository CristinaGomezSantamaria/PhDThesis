%USE WITH ComparisonChPhyFreqSel
%This function automatically plots angular vs delay spectrums for different
%types of clusters configurations in the channel
global cs

H = H_Total(cs);
h = h_Total(cs);
for bs=1:Escenario.num_BS
    if ~isempty(Escenario.deseados{bs})
    for ms=1:Escenario.N_d{bs}
        %PLOT OF ANGLES vs DELAY
        switch cs
            case 1
                figura = figure;
                figure(figura);
                h_AoDazim = polar(h.DDCIR{bs,ms}(3,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk'); hold on;
                figure(figura+1);
                h_AoDelev = polar(h.DDCIR{bs,ms}(4,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk'); hold on;
                figure(figura+2);
                h_AoAazim = polar(h.DDCIR{bs,ms}(5,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk'); hold on;
                figure(figura+3);
                h_AoAelev = polar(h.DDCIR{bs,ms}(6,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk'); hold on;
            case 2
                figura = figure;
                figure(figura);
                h_AoDazim = polar(h.DDCIR{bs,ms}(3,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+1);
                h_AoDelev = polar(h.DDCIR{bs,ms}(4,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+2);
                h_AoAazim = polar(h.DDCIR{bs,ms}(5,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+3);
                h_AoAelev = polar(h.DDCIR{bs,ms}(6,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');                
            case 3
                figura = figure;
                figura = figure;
                figure(figura);
                h_AoDazim = polar(h.DDCIR{bs,ms}(3,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+1);
                h_AoDelev = polar(h.DDCIR{bs,ms}(4,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+2);
                h_AoAazim = polar(h.DDCIR{bs,ms}(5,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+3);
                h_AoAelev = polar(h.DDCIR{bs,ms}(6,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
            case 4
                figura = figure;
                figure(figura);
                h_AoDazim = polar(h.DDCIR{bs,ms}(3,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+1);
                h_AoDelev = polar(h.DDCIR{bs,ms}(4,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+2);
                h_AoAazim = polar(h.DDCIR{bs,ms}(5,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
                figure(figura+3);
                h_AoAelev = polar(h.DDCIR{bs,ms}(6,:)*pi/180,h.DDCIR{bs,ms}(1,:),'dk');
        end
    end
    end
end