Long Term Evolution. Proyecto Final. 
===========================================================================================================================
===========================================================================================================================
1. Instrucciones de configuraci�n
Descargar Matlab r2019b.
Descargar la paqueter�a LTE-Toolbox. 
Clonar el siguiente repositorio en Documents/MATLAB
	https://github.com/javierMoG/Long-Term-Evolution-LTE-.git
===========================================================================================================================
===========================================================================================================================
2. Instrucciones de operaci�n
La simulaci�n principal del modelo se hace mediante el archivo 
	simulaci�n.m
Se configuran s�lo los siguientes par�metros
	nF: N�mero de frames totales a simular. 
	espSnr: Espaciaci�n de los valores del SNR a simular. La simulaci�n corre de 0 a 20 Db con una espaciaci�n de espSNR
	modulacion: Modulaci�nes a simular. S�lo pueden ser 'QPSK', '16QAM' y/o '64QAM'. 

Para correr la simulaci�n en Matlab es necesario situarse en Documents/Matlab/Long-Term-Evolution-LTE- y escribir en la
ventana de comandos el siguiente comando
	>>simulacion
===========================================================================================================================
===========================================================================================================================
3. Manifiesto
Software.
	Matlab r2019b
`	Paqueteria LTE-Toolbox

Archivos dentro de la carpeta.
	hDiaplayENBParameterSummary.m
	hHARQScheduling.m
	hNewHARQProcess.m
	hPlotRVSequence.m
	PHY.m
	simulacion.m
	BER.fig
	Thr.fig
===========================================================================================================================
===========================================================================================================================
4. Informaci�n sobre los desarrolladores
	Javier Montiel Gonz�lez. 	javiermg19978@gmail.com
	Alexis Calvillo Madrid.  	acalvillm@gmail.com
