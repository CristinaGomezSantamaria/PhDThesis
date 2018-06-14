%CANAL

function [CH]=channel_gen()

global Escenario bs

switch Escenario.channel
    case 'IIDRayleigh'
        for ms=1:Escenario.N_d{bs}
            CH.hS{ms,1}=zeros(Escenario.P,Escenario.N);
            CH.Rtx=eye(Escenario.N);
            CH.Rrx=eye(Escenario.P);
            CH.R=kron(CH.Rtx,CH.Rrx);
            X=(randn(Escenario.N*Escenario.P,1)+j*randn(Escenario.N*Escenario.P,1))*Escenario.Pot;         
            CH.hS{ms,1}=reshape(CH.R'*X,Escenario.P,Escenario.N);  
        end
    case 'CorrRayleigh'
        for ms=1:Escenario.N_d{bs}
            CH.hS{ms,1}=zeros(Escenario.P,Escenario.N);     
            for pp=1:Escenario.N
                kte(pp,:)=2*pi*(pp-(1:Escenario.N))*Escenario.A_S*Escenario.d_BS;
                CH.Rtx(pp,:)=(1/2*pi)*quadv(@(x)myfunction(x,kte(pp,:)),0,2*pi);
            end
            CH.Rtx=CH.Rtx/max(max(CH.Rtx));
%             figure;
%             plot(kte,CH.Rtx,'-*b');
%             hold on; grid on;
%             J=besselj(0,kte);
%             plot(kte,J,'--r');
%             CH.Rtx=eye(Escenario.N)+fliplr(eye(Escenario.N))*Escenario.ro;
%             útil para insertar Rtx manualmente
            CH.Rrx=eye(Escenario.P);
            CH.R=kron(CH.Rtx,CH.Rrx);
            X=(randn(Escenario.N*Escenario.P,1)+j*randn(Escenario.N*Escenario.P,1))*Escenario.Pot;         
            CH.hS{ms,1}=reshape(CH.R'*X,Escenario.P,Escenario.N);  
        end
end