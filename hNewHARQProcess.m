% hNewHARQProcess Initialize HARQ processes
%
%   HARQPROCESSES = hNewHARQProcess(ENBUE) es una estructura que contiene
%   los parametros incializados de los procesos HARQ (Hybrid Automatic 
%   Repeat Request) dado una configuracion ENBUE de enlace descendente, 
%   ascendente o enlance lateral
%
%   HARQPROCESSES = hNewHARQProcess(ENB,PDSCH) inicializa los procesos HARQinitializes HARQ processes
%   para una estructura ENB y configuracion de estructura PDSCH dados. 

%   Copyright 2015-2018 The MathWorks, Inc.

function harqProcesses = hNewHARQProcess(enbue,varargin)
    
    if(nargin>1)
        enbue.PDSCH = varargin{1};
    end
    % Identifica la estructura que contiene la entidad de parametros HARQ
    % (PDSCH en enlace descente, PUSCH en enlace ascendente, UE en enlace
    % lateral)
    if isfield(enbue,'PDSCH')
        harqstruct = enbue.PDSCH;
    elseif isfield(enbue,'PUSCH')
        harqstruct = enbue.PUSCH;
    else
        harqstruct = enbue; % La estructura de entrada es 'UE'
    end
    
    % Si la Modulacion es un arreglo de celdad multiples palabras codigo
    % son usadas
    if (iscell(harqstruct.Modulation)) || isstring(harqstruct.Modulation)
        ncw = length(harqstruct.Modulation);
    else
        ncw = 1;
    end
    harqProcess.data = {};                     % Inicializamos la informacion
    harqProcess.iniConfig = harqstruct;        % Inicializamos la configuracion
    harqProcess.iniConfig.RVIdx = ones(1,ncw); % Add RVIdx to config Agregamos RVIdx a la configuracion
    harqProcess.blkerr = zeros(1,ncw);         % Inicializamos el bloque de errores
    harqProcess.ncw = ncw;                     % Establecemos el numero de palabras codigo
    
    % Inicializamos la configuracion de la estructira transmitir
    harqProcess.txConfig = harqProcess.iniConfig;
    
    % Creamos una estructura vacia para cada palabra codigo
    harqProcess.decState(1:ncw) = deal(struct());
    
    % Creamos los procesos HARQ como los indica NHARQProcesses
    harqProcesses(1:harqstruct.NHARQProcesses) = harqProcess;
end