%This function plots other results from the Physical Channel which require
%statistical measures over several samples of the Matrix Channel
MAXITER = 1000;
snr_vec = 15;

%CAPACITY
C = zeros(MAXITER,Escenario.num_MS);
for bs=1:Escenario.num_BS
    for ms=1:Escenario.N_d{bs}
        for countiter=1:MAXITER
            H_temp = H_Total(countiter).Tap_norm_awgn{bs,ms};
            C(countiter,Escenario.deseados{bs}(ms)) = abs(log2(det(eye(Escenario.P)+(snr_vec/Escenario.N)*H_temp*H_temp')));
        end
    end
end
x = min(min(C)):max(max(C))/20:max(max(C));
distr_prob = hist(C,x);
if Escenario.num_MS == 1
    distr_prob_acum = cumsum(distr_prob)/sum(distr_prob);
else
    distr_prob_acum = cumsum(distr_prob)./kron(sum(distr_prob),ones(length(x),1));
end
figure;
bar(x,distr_prob_acum,0.5,'grouped')
figure;
plot(x,distr_prob_acum);

%ANGLE SPREAD: GENEAL FOR EACH LINK
%AoD Azimuthal
x_az = 0:15:150;
AoDazAS_temp = zeros(MAXITER,Escenario.num_MS);
for countiter=1:MAXITER
    AoDazAS_temp(countiter,:) = h_Total(countiter).A_S_AoD_az;
end
figure;
pdf_AoDaz = hist(AoDazAS_temp,x_az);
if Escenario.num_MS == 1
    cdf_AoDaz = cumsum(pdf_AoDaz)/sum(pdf_AoDaz);
else
    cdf_AoDaz = cumsum(pdf_AoDaz)./kron(sum(pdf_AoDaz),ones(length(x_az),1));
end
bar(x_az,cdf_AoDaz,0.5,'grouped');
title('CDF Angle Spread for Azimuthal AoD');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
plot(x_az,cdf_AoDaz);
title('CDF Angle Spread for Azimuthal AoD');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
bar(x_az,pdf_AoDaz,0.5,'grouped');
title('PDF Angle Spread for Azimuthal AoD');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');

%AoD Elevation
x_el = 20:5:50;
AoDelAS_temp = zeros(MAXITER,Escenario.num_MS);
for countiter=1:MAXITER
    AoDelAS_temp(countiter,:) = h_Total(countiter).A_S_AoD_el;
end
figure;
pdf_AoDel = hist(AoDelAS_temp,x_el);
if Escenario.num_MS == 1
    cdf_AoDel = cumsum(pdf_AoDel)/sum(pdf_AoDel);
else
    cdf_AoDel = cumsum(pdf_AoDel)./kron(sum(pdf_AoDel),ones(length(x_el),1));
end
bar(x_el,cdf_AoDel,0.5,'grouped')
title('CDF Angle Spread for Elevation AoD');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
plot(x_el,cdf_AoDel);
title('CDF Angle Spread for Elevation AoD');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
bar(x_el,pdf_AoDel,0.5,'grouped');
title('PDF Angle Spread for Elevation AoD');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');

%AoA Azimuthal
x_az = 0:30:180;
AoAazAS_temp = zeros(MAXITER,Escenario.num_MS);
for countiter=1:MAXITER
    AoAazAS_temp(countiter,:) = h_Total(countiter).A_S_AoA_az;
end
figure;
pdf_AoAaz = hist(AoAazAS_temp,x_az);
if Escenario.num_MS == 1
    cdf_AoAaz = cumsum(pdf_AoAaz)/sum(pdf_AoAaz);
else
    cdf_AoAaz = cumsum(pdf_AoAaz)./kron(sum(pdf_AoAaz),ones(length(x_az),1));
end
bar(x_az,cdf_AoAaz,0.5,'grouped')
title('CDF Angle Spread for Azimuthal AoA');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
plot(x_az,cdf_AoAaz);
title('CDF Angle Spread for Azimuthal AoA');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
bar(x_az,pdf_AoAaz,0.5,'grouped');
title('PDF Angle Spread for Azimuthal AoA');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');

%AoA Elevation
x_el = -10:5:50;
AoAelAS_temp = zeros(MAXITER,Escenario.num_MS);
for countiter=1:MAXITER
    AoAelAS_temp(countiter,:) = h_Total(countiter).A_S_AoA_el;
end
figure;
pdf_AoAel = hist(AoAelAS_temp,x_el);
if Escenario.num_MS == 1
    cdf_AoAel = cumsum(pdf_AoAel)/sum(pdf_AoAel);
else
    cdf_AoAel = cumsum(pdf_AoAel)./kron(sum(pdf_AoAel),ones(length(x_el),1));
end
bar(x_el,cdf_AoAel,0.5,'grouped')
title('CDF Angle Spread for Elevation AoA');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
plot(x_el,cdf_AoAel);
title('CDF Angle Spread for Elevation AoA');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
figure;
bar(x_el,pdf_AoAel,0.5,'grouped');
title('PDF Angle Spread for Elevation AoA');
xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');

%ANGLE SPREAD: SPECIFICALLY FOR EACH CLUSTER --> AROUND MS, BS, FC
%AROUND MS
if strcmp(Escenario.MS_cluster.id,'si')
    %AoD Azimuthal
    x_az = -0.01:0.005:0.03;
    AoDazMS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoDazMS(countiter,:) = h_Total(countiter).A_S_AoD_az_ms;
    end
    figure;
    pdf_AoDazMS = hist(AoDazMS,x_az);
    if Escenario.num_MS == 1
        cdf_AoDazMS = cumsum(pdf_AoDazMS)/sum(pdf_AoDazMS);
    else
        cdf_AoDazMS = cumsum(pdf_AoDazMS)./kron(sum(pdf_AoDazMS),ones(length(x_az),1));
    end
    bar(x_az,cdf_AoDazMS,0.5,'grouped');
    title('CDF Angle Spread for Azimuthal AoD, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_az,cdf_AoDazMS);
    title('CDF Angle Spread for Azimuthal AoD, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_az,pdf_AoDazMS,0.5,'grouped');
    title('PDF Angle Spread for Azimuthal AoD, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoD Elevation
    x_el = -0.01:0.005:0.03;
    AoDelMS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoDelMS(countiter,:) = h_Total(countiter).A_S_AoD_el_ms;
    end
    figure;
    pdf_AoDelMS = hist(AoDelMS,x_el);
    if Escenario.num_MS == 1
        cdf_AoDelMS = cumsum(pdf_AoDelMS)/sum(pdf_AoDelMS);
    else
        cdf_AoDelMS = cumsum(pdf_AoDelMS)./kron(sum(pdf_AoDelMS),ones(length(x_el),1));
    end
    bar(x_el,cdf_AoDelMS,0.5,'grouped');
    title('CDF Angle Spread for Elevation AoD, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_el,cdf_AoDelMS);
    title('CDF Angle Spread for Elevation AoD, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_el,pdf_AoDelMS,0.5,'grouped');
    title('PDF Angle Spread for Elevation AoD, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoA Azimuthal
    x_az = -30:5:180;
    AoAazMS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoAazMS(countiter,:) = h_Total(countiter).A_S_AoA_az_ms;
    end
    figure;
    pdf_AoAazMS = hist(AoAazMS,x_az);
    if Escenario.num_MS == 1
        cdf_AoAazMS = cumsum(pdf_AoAazMS)/sum(pdf_AoAazMS);
    else
        cdf_AoAazMS = cumsum(pdf_AoAazMS)./kron(sum(pdf_AoAazMS),ones(length(x_az),1));
    end
    bar(x_az,cdf_AoAazMS,0.5,'grouped');
    title('CDF Angle Spread for Azimuthal AoA, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_az,cdf_AoAazMS);
    title('CDF Angle Spread for Azimuthal AoA, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_az,pdf_AoAazMS,0.5,'grouped');
    title('PDF Angle Spread for Azimuthal AoA, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoA Elevation
    x_el = 0:0.05:0.1;
    AoAelMS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoAelMS(countiter,:) = h_Total(countiter).A_S_AoA_el_ms;
    end
    figure;
    pdf_AoAelMS = hist(AoAelMS,x_el);
    if Escenario.num_MS == 1
        cdf_AoAelMS = cumsum(pdf_AoAelMS)/sum(pdf_AoAelMS);
    else
        cdf_AoAelMS = cumsum(pdf_AoAelMS)./kron(sum(pdf_AoAelMS),ones(length(x_el),1));
    end
    bar(x_el,cdf_AoAelMS,0.5,'grouped');
    title('CDF Angle Spread for Elevation AoA, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_el,cdf_AoAelMS);
    title('CDF Angle Spread for Elevation AoA, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_el,pdf_AoAelMS,0.5,'grouped');
    title('PDF Angle Spread for Elevation AoA, Local MS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
end

%AROUND BS
if strcmp(Escenario.BS_cluster.id,'si')
    %AoD Azimuthal
    x_az = 0:10:150;
    AoDazBS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoDazBS(countiter,:) = h_Total(countiter).A_S_AoD_az_bs;
    end
    figure;
    pdf_AoDazBS = hist(AoDazBS,x_az);
    if Escenario.num_MS == 1
        cdf_AoDazBS = cumsum(pdf_AoDazBS)/sum(pdf_AoDazBS);
    else
        cdf_AoDazBS = cumsum(pdf_AoDazBS)./kron(sum(pdf_AoDazBS),ones(length(x_az),1));
    end
    bar(x_az,cdf_AoDazBS,0.5,'grouped');
    title('CDF Angle Spread for Azimuthal AoD, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_az,cdf_AoDazBS);
    title('CDF Angle Spread for Azimuthal AoD, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_az,pdf_AoDazBS,0.5,'grouped');
    title('PDF Angle Spread for Azimuthal AoD, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoD Elevation
    x_el = -5:15:100;
    AoDelBS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoDelBS(countiter,:) = h_Total(countiter).A_S_AoD_el_bs;
    end
    figure;
    pdf_AoDelBS = hist(AoDelBS,x_el);
    if Escenario.num_MS == 1
        cdf_AoDelBS = cumsum(pdf_AoDelBS)/sum(pdf_AoDelBS);
    else
        cdf_AoDelBS = cumsum(pdf_AoDelBS)./kron(sum(pdf_AoDelBS),ones(length(x_el),1));
    end
    bar(x_el,cdf_AoDelBS,0.5,'grouped');
    title('CDF Angle Spread for Elevation AoD, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_el,cdf_AoDelBS);
    title('CDF Angle Spread for Elevation AoD, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_el,pdf_AoDelBS,0.5,'grouped');
    title('PDF Angle Spread for Elevation AoD, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoA Azimuthal
    x_az = -0.1:0.05:0.1;
    AoAazBS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoAazBS(countiter,:) = h_Total(countiter).A_S_AoA_az_bs;
    end
    figure;
    pdf_AoAazBS = hist(AoAazBS,x_az);
    if Escenario.num_MS == 1
        cdf_AoAazBS = cumsum(pdf_AoAazBS)/sum(pdf_AoAazBS);
    else
        cdf_AoAazBS = cumsum(pdf_AoAazBS)./kron(sum(pdf_AoAazBS),ones(length(x_az),1));
    end
    bar(x_az,cdf_AoAazBS,0.5,'grouped');
    title('CDF Angle Spread for Azimuthal AoA, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_az,cdf_AoAazBS);
    title('CDF Angle Spread for Azimuthal AoA, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_az,pdf_AoAazBS,0.5,'grouped');
    title('PDF Angle Spread for Azimuthal AoA, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoA Elevation
    x_el = 0:0.05:0.2;
    AoAelBS = zeros(MAXITER,Escenario.num_MS);
    for countiter=1:MAXITER
        AoAelBS(countiter,:) = h_Total(countiter).A_S_AoA_el_bs;
    end
    figure;
    pdf_AoAelBS = hist(AoAelBS,x_el);
    if Escenario.num_MS == 1
        cdf_AoAelBS = cumsum(pdf_AoAelBS)/sum(pdf_AoAelBS);
    else
        cdf_AoAelBS = cumsum(pdf_AoAelBS)./kron(sum(pdf_AoAelBS),ones(length(x_el),1));
    end
    bar(x_el,cdf_AoAelBS,0.5,'grouped');
    title('CDF Angle Spread for Elevation AoA, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_el,cdf_AoAelBS);
    title('CDF Angle Spread for Elevation AoA, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_el,pdf_AoAelBS,0.5,'grouped');
    title('PDF Angle Spread for Elevation AoA, Local BS Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
end

%AROUND FC
if strcmp(Escenario.FC_cluster.id,'si')
    %AoD Azimuthal
    x_az = 0:0.5:5;
    AoDazFC = zeros(MAXITER,length(h_Total(1).A_S_AoD_az_fc));
    for countiter=1:MAXITER
        AoDazFC(countiter,:) = h_Total(countiter).A_S_AoD_az_fc;
    end
    figure;
    pdf_AoDazFC = hist(AoDazFC,x_az);
    if Escenario.num_MS == 1
        cdf_AoDazFC = cumsum(pdf_AoDazFC)/sum(pdf_AoDazFC);
    else
        cdf_AoDazFC = cumsum(pdf_AoDazFC)./kron(sum(pdf_AoDazFC),ones(length(x_az),1));
    end
    bar(x_az,cdf_AoDazFC,0.5,'grouped');
    title('CDF Angle Spread for Azimuthal AoD, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_az,cdf_AoDazFC);
    title('CDF Angle Spread for Azimuthal AoD, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_az,pdf_AoDazFC,0.5,'grouped');
    title('PDF Angle Spread for Azimuthal AoD, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoD Elevation
    x_el = 0:0.5:5;
    AoDelFC = zeros(MAXITER,length(h_Total(1).A_S_AoD_az_fc));
    for countiter=1:MAXITER
        AoDelFC(countiter,:) = h_Total(countiter).A_S_AoD_el_fc;
    end
    figure;
    pdf_AoDelFC = hist(AoDelFC,x_el);
    if Escenario.num_MS == 1
        cdf_AoDelFC = cumsum(pdf_AoDelFC)/sum(pdf_AoDelFC);
    else
        cdf_AoDelFC = cumsum(pdf_AoDelFC)./kron(sum(pdf_AoDelFC),ones(length(x_el),1));
    end
    bar(x_el,cdf_AoDelFC,0.5,'grouped');
    title('CDF Angle Spread for Elevation AoD, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_el,cdf_AoDelFC);
    title('CDF Angle Spread for Elevation AoD, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_el,pdf_AoDelFC,0.5,'grouped');
    title('PDF Angle Spread for Elevation AoD, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoA Azimuthal
    x_az = -0.1:0.05:0.1;
    AoAazFC = zeros(MAXITER,length(h_Total(1).A_S_AoD_az_fc));
    for countiter=1:MAXITER
        AoAazFC(countiter,:) = h_Total(countiter).A_S_AoA_az_fc;
    end
    figure;
    pdf_AoAazFC = hist(AoAazFC,x_az);
    if Escenario.num_MS == 1
        cdf_AoAazFC = cumsum(pdf_AoAazFC)/sum(pdf_AoAazFC);
    else
        cdf_AoAazFC = cumsum(pdf_AoAazFC)./kron(sum(pdf_AoAazFC),ones(length(x_az),1));
    end
    bar(x_az,cdf_AoAazFC,0.5,'grouped');
    title('CDF Angle Spread for Azimuthal AoA, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_az,cdf_AoAazFC);
    title('CDF Angle Spread for Azimuthal AoA, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_az,pdf_AoAazFC,0.5,'grouped');
    title('PDF Angle Spread for Azimuthal AoA, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    
    %AoA Elevation
    x_el = -0.02:0.005:0.02;
    AoAelFC = zeros(MAXITER,length(h_Total(1).A_S_AoD_az_fc));
    for countiter=1:MAXITER
        AoAelFC(countiter,:) = h_Total(countiter).A_S_AoA_el_fc;
    end
    figure;
    pdf_AoAelFC = hist(AoAelFC,x_el);
    if Escenario.num_MS == 1
        cdf_AoAelFC = cumsum(pdf_AoAelFC)/sum(pdf_AoAelFC);
    else
        cdf_AoAelFC = cumsum(pdf_AoAelFC)./kron(sum(pdf_AoAelFC),ones(length(x_el),1));
    end
    bar(x_el,cdf_AoAelFC,0.5,'grouped');
    title('CDF Angle Spread for Elevation AoA, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    plot(x_el,cdf_AoAelFC);
    title('CDF Angle Spread for Elevation AoA, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
    figure;
    bar(x_el,pdf_AoAelFC,0.5,'grouped');
    title('PDF Angle Spread for Elevation AoA, FC Cluster');
    xlabel('A_S (°)'); ylabel('Prob A_S<Abscissa');
end


%DELAY SPREAD: GENEAL FOR EACH LINK
x_delay = 0:H_Total(1).Tsymbol/5:H_Total(1).Max_Retardo;
Delay_temp = zeros(MAXITER,Escenario.num_MS);
for countiter=1:MAXITER
    Delay_temp(countiter,:) = H_Total(countiter).D_S_norm2;
end
figure;
pdf_delay = hist(Delay_temp,x_delay);
if Escenario.num_MS == 1
    cdf_delay = cumsum(pdf_delay)/sum(pdf_delay);
else
    cdf_delay = cumsum(pdf_delay)./kron(sum(pdf_delay),ones(length(x_delay),1));
end
bar(x_delay,cdf_delay,0.5,'grouped');
title('CDF Delay Spread');
xlabel('Delay Spread (seg)'); ylabel('Prob A_S<Abscissa');
figure;
plot(x_delay,cdf_delay);
title('CDF Delay Spread');
xlabel('Delay Spread (seg)'); ylabel('Prob A_S<Abscissa');
figure;
bar(x_delay,pdf_delay,0.5,'grouped');
title('PDF Delay Spread');
xlabel('Delay Spread (seg)'); ylabel('Prob A_S<Abscissa');

