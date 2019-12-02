%Parámetros de simulación
nFr=300;
espSnr=.5;
modulacion=["16QAM";"64QAM";"QPSK"];
SNRIn=[];
trgh=[];
BER=[];
for i=1:length(modulacion)
    [SNRIn(i,:),BER(i,:),trgh(i,:)]=PHY(modulacion(i),nFr,espSnr,false);
end
figure(1)
plot(SNRIn(1,:),trgh(1,:),'*-.')
grid on;
hold on;
plot(SNRIn(2,:),trgh(2,:),'*-.')
plot(SNRIn(3,:),trgh(3,:),'*-.')
%plot(SNRIn, 1e-6*simThroughput/(NFrames*10e-3),'*-.');
xlabel('SNR (dB)');
ylabel('Throughput (Mbps)');
legend(modulacion,'Location','NorthWest');
figure(2)
plot(SNRIn(1,:), BER(1,:));
grid on;
hold on;
plot(SNRIn(2,:), BER(2,:));
plot(SNRIn(3,:), BER(3,:));
xlabel('SNR (dB)');
ylabel('BER');
legend(modulacion);
