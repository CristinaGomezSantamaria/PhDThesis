%Performance of the Variable Random Generation for the position of the
%scatterers with pdf given by Molisch

close all;
Num_Disp = 10000;
k1 = 1/15;
radio = 100;
k2 = 1/(2*radio^2);
count = 0;
count_total = 0;
r_vec = zeros(1,Num_Disp);

ptos_graph = 1000;
x_graph = 1:ptos_graph;
f_x_graph = k1*exp(-k2*x_graph.^2);
Gauss_clasica = (1/sqrt(2*pi*100^2))*exp(-k2*x_graph.^2);

while count < sum(Num_Disp)
    count_total = count_total+1;
    x = randint(1,1,[0 1000]);                                                                                            
    u = k1*rand;
    f_x = k1*exp(-k2*x^2);
    if u <= f_x
        r=x;
        r_vec(count+1) = r;
        count=count+1;
    end
end
[n,r] = hist(r_vec,50);
figure;
h_bar = bar(r,n/Num_Disp); hold on;
h_plot = plot(x_graph,f_x_graph,'k--',x_graph,Gauss_clasica,'r--');
legend([h_bar h_plot(1) h_plot(2)],'PDF R.V.Generated','Theorical PDF by Molisch','Classical Gauss')

