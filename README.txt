Long Term Evolution. Proyecto Final. 
===========================================================================================================================
===========================================================================================================================
1. Instrucciones de configuración
Descargar Matlab r2019b.
Descargar la paquetería LTE-Toolbox. 
Clonar el siguiente repositorio en Documents/MATLAB
	https://github.com/javierMoG/Long-Term-Evolution-LTE-.git
===========================================================================================================================
===========================================================================================================================
2. Instrucciones de operación
La simulación principal del modelo se hace mediante el archivo 
	simulación.m
Se configuran sólo los siguientes parámetros
	nF: Número de frames totales a simular. 
	espSnr: Espaciación de los valores del SNR a simular. La simulación corre de 0 a 20 Db con una espaciación de espSNR
	modulacion: Modulaciónes a simular. Sólo pueden ser 'QPSK', '16QAM' y/o '64QAM'. 

Para correr la simulación en Matlab es necesario situarse en Documents/Matlab/Long-Term-Evolution-LTE- y escribir en la
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
4. Información sobre los desarrolladores
	Javier Montiel González. 	javiermg19978@gmail.com
	Alexis Calvillo Madrid.  	acalvillm@gmail.com
