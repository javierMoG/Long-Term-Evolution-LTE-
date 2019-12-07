%Parámetros de simulación
tic
%Número de frames simulados
nFr=10;
%Espaciación entre valores del SNR
espSnr=.5;
%Modulaciones a simular
modulacion=["16QAM";"64QAM";"QPSK"];
%Vector de SNR con los valores del SNR simulado
SNRIn=[];
%Vector de throughput con los valores del throughput simulado (en Mbps)
trgh=[];
%Vector de Bit Error Rate con los valores de Bit Error Rate simulado
BER=[];
for i=1:length(modulacion)
    [SNRIn(i,:),BER(i,:),trgh(i,:)]=PHY(modulacion(i),nFr,espSnr,false);
end
toc

%Gráfica de resultados de Throughput
figure(1)
plot(SNRIn(1,:),trgh(1,:),'*-.')
grid on;
hold on;
plot(SNRIn(2,:),trgh(2,:),'*-.')
plot(SNRIn(3,:),trgh(3,:),'*-.')
xlabel('SNR (dB)');
ylabel('Throughput (Mbps)');
legend(modulacion,'Location','NorthWest');
%Gráfica de resultados de BER
figure(2)
plot(SNRIn(1,:), BER(1,:));
grid on;
hold on;
plot(SNRIn(2,:), BER(2,:));
plot(SNRIn(3,:), BER(3,:));
xlabel('SNR (dB)');
ylabel('BER');
legend(modulacion);
