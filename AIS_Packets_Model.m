 clear; clc;
% Исходные параметры 

% Ненастраевымые параметры
% Частота дискретизации, Гц
Fs = 249600;
% Длительность одного пакета АИС в битах (26.67 мс)
PacketLen_bit = 256; 
% Символьная скорость передачи, б/c
SymbolRate = 9600; 
% Количество отсчётов на один символ (бит, т.к ФМ-2). Формула round(Fs / SymbolRate)
sps = 26; 
% Коэффициент фильтра в GMSK модуляторе
filterBT = 1.0;
% Количество отсчётов в символьных интервалах в GMSK модуляторе
span = 1;

% Принятый сигнал АИС
% Длина одного АИС пакета в отсчётах 
PacketLen_smpls = PacketLen_bit*sps;  
% Число пакетов АИС в записи, шт
NumPackets = 10;  
% Длина записи в отсчётах
RecordLen_smps = NumPackets*PacketLen_smpls;

% Вектор под запись демодулированных сигналов
RecordAisSignal = zeros(RecordLen_smps*100, 1); 
%RecordAisSignal = zeros(RecordLen_smps - PacketLen_smpls + 1 + NumPackets*PacketLen_smpls, 1);

% Длительность паузы между сообщениями в записи, в длительностях сообщения АИС
PauseLen = 0; 

% Массив со сдвигами по частоте
deltaF = zeros(NumPackets, 1); 
% Массив с индексами начала пакетов
StartPositions = zeros(NumPackets, 1); 
% Массив с индексами конца пакетов
EndPositions = zeros(NumPackets, 1); 
% Массив с количеством коллизий для каждого пакета
CollisionsCounter = zeros(NumPackets, 1); 


% Принятый сигнал
RecievedAisSignal = 2*randi ([0 1], 1, RecordLen_smps) - 1;


modGMSK = comm.GMSKModulator('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitInput',true,'PulseLength',...
    span);

demodGMSK = comm.GMSKDemodulator ('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitOutput',true,'PulseLength',...
    span, 'TracebackDepth', 16);




for i = 1 : NumPackets
    % Формирование принятого пакета
    % Информационные биты
    BitPacket(i,:) = randi ([0 1], PacketLen_bit, 1);
    % Формирование комплексной огибающей АИС-сигнала (GMSK модуляция)
    %PacketLPE = GMSKmod(BitPacket(i,:), span, sps, filterBT);
    PacketLPE = modGMSK(BitPacket(i,:)');
    
    % Допплеровский сдвиг (пока только сдвиг частоты)
    % Нижняя граница доплеровского сдвига, Гц
    dF_Bottom = -4000;
    % Верхняя граница доплеровского сдвига, Гц 
    dF_Upper = 4000; 
    % Случайная величина частотного сдвига
    deltaFreq(i) = randi([dF_Bottom, dF_Upper]); 
    % Сигнал (пакет), сдвинутый по частоте
    FreqShiftedAisPacketLPE = PacketLPE .* exp(1i*2*pi*deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';

    % Случайный поворот фазы
    % Phi = -pi/2 + pi*rand; % от -pi/2 до pi/2
    % AisLPE = AisLPE * exp(1i*Phi);
    
    % Случайное ослабление сигнала
    % AisLPE_deltaF = randi([1000, 30000]) .* AisLPE_deltaF;


    % Коллизии
    % Маяк наличия коллизии
    %CollisionsOn = randi([0 1]);


    % Добавление текущего сообщения в общую запись на случайное место
    %if CollisionsOn
        StartPositions(i) = randi([1, RecordLen_smps*100 - PacketLen_smpls + 1]);
    % else
    %     if i == 1
    %         StartPositions(i) = 1;
    %     else
    %         StartPositions(i) = StartPositions(i-1) + PacketLen_smpls * (1 + PauseLen);
    %     end
    % end

    EndPositions(i) = StartPositions(i) + PacketLen_smpls - 1;
    % Запишем текущую КО на места (start:end) с учётом уже разположеных там значений
    RecordAisSignal(StartPositions(i) : EndPositions(i)) = (RecordAisSignal(StartPositions(i) : EndPositions(i)) + PacketLPE);

end


% Нормировка результата
    RecordAisSignal = RecordAisSignal / max([abs(real(RecordAisSignal)); ...
        abs(imag(RecordAisSignal))]);

% Фиксация коллизий
for j = 1 : NumPackets
    % Зафиксировали один пакет j и сравниваем с ним все оставшиеся пакеты k 
    for k = 1 : NumPackets
        if k ~= j
            if StartPositions(j) >= StartPositions(k) && StartPositions(j) <= EndPositions(k) % либо наложились, либо j-й позже k-го
                CollisionsCounter(j) = CollisionsCounter(j) + 1;
            end
            if EndPositions(j) >= StartPositions(k) && EndPositions(j) < EndPositions(k) % либо наложилось, лтбо k-й позже j-го
                CollisionsCounter(j) = CollisionsCounter(j) + 1;
            end
        end
    end
end


demodData = zeros (10, 256);
%%%%% Демодуляция
for i = 1:NumPackets
    SignalToDemod = RecordAisSignal(StartPositions(i):EndPositions(i));
    
    %demodData(i,:) = GMSKDemod(SignalToDemod, sps, filterBT, span);
    %demodGMSK (i,:) = demodGMSK(SignalToDemod);
    demodData =  demodGMSK(SignalToDemod);
    OutBitPacket (i, :) = demodData';
    [errorCounter(i,:), ratio(i,:)] = biterr (BitPacket(1,:), OutBitPacket(i,:));

    
end 



% % Зонные фильтры (дб внутри цикла демодуляции)
% Fs = 250000;  % Sampling Frequency
% Fpass = 4800;            % Passband Frequency
% Fstop = 9000;            % Stopband Frequency
% Dpass = 0.057501127785;  % Passband Ripple
% Dstop = 0.0001;          % Stopband Attenuation
% dens  = 20;              % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% ImpulseResponse  = firpm(N, Fo, Ao, W, {dens});
% 
% 
% bandWidth = 1/Params.symbolInterval; % Hz;
% 
% t = (0:length(SignalToDemod)-1)/Fs;
% fc1 = -3840.1; fc2 = 0; fc3 = 3840.1;
% 
% narrow_band1_leftpart = conv(SignalToDemod.*exp(-1i*2*pi*fc1.*t)',flip(ImpulseResponse));
% z1 = narrow_band1_leftpart.*exp(1i*2*pi*fc1.*(0:length(narrow_band1_leftpart)-1)/Fs)';
% 
% narrow_band2_centered = conv(SignalToDemod,flip(ImpulseResponse));
% z2 = narrow_band2_centered;
% 
% narrow_band3_rightpart = conv(SignalToDemod.*exp(-1i*2*pi*fc3.*t)',flip(ImpulseResponse));
% z3 = narrow_band3_rightpart.*exp(1i*2*pi*fc3.*(0:length(narrow_band3_rightpart)-1)/Fs)';
% 
% ZonalBandSignals = [z1, z2, z3];

    






function s = GMSKmod(input, span, sps, filterBT)
gmskModulator; 
% N = length(input);
% h = gaussdesign(filterBT, span, sps);
% h = h/std(h);
% 
% % NRZI кодирование; когда приходит "0" меняется, "1" - не меняется
% SymNrzi = zeros(length(input), 1);
% SymNrzi(1) = 2 * input(1) - 1;
% for i = 2 : length(input)
%     if (input(i) == 0)
%         if (SymNrzi(i - 1) == 1)
%             SymNrzi(i) = -1;
%         else
%             SymNrzi(i) = 1;
%         end
%     else
%         SymNrzi(i) = SymNrzi(i - 1);
%     end
% end
% 
% symbols = repmat(SymNrzi, 1, sps);
% symbols = reshape(symbols.', sps*N, 1); % эти строки это repelem(SymNrzi, samplesPerSymbol);
% 
% base = filter(h, 1, symbols);
% base = base/max(abs(base))/sps;
% integrated = zeros(length(base), 1);
% integrated(1) = base(1);
%     for i = 2 : length(base)
%         integrated(i) = integrated(i-1) + base(i);
%     end
% 
% wd = 0.5*pi; % девиация?
% inphase = cos(wd*integrated);
% quadrature = sin(wd*integrated);
% s = inphase + 1i*quadrature;
end

function outBits  = GMSKDemod(signal, sps, filterBT, span)

trainingSeqNRZI = logical([1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1]);
modGMSK = comm.GMSKModulator('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitInput',true,'PulseLength',...
    span,'SymbolPrehistory', [1 -1]);

prbdet = comm.PreambleDetector('Input','Symbol','Detections','First');
prbdet.Preamble = (2*[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]-1).';
prbdet.Threshold = 20;


syncSyms = modGMSK(trainingSeqNRZI');
syncAngles = unwrap(angle(syncSyms));
rcvAngles = unwrap(angle(signal));
lenCorrWin = 2*length(syncAngles);

pos = 32;
    % for i = 1 : 1 : length(rcvAngles)-2*length(syncAngles)       
    %     if (length(rcvAngles) > lenCorrWin)
    %         [acor,lag] = xcorr(syncAngles,rcvAngles(i : i + lenCorrWin - 1));
    %         [~,I] = max(abs(acor));
    %         lagDiff = lag(I);
    %         Index = lagDiff;
    %         pos = -Index+1;
    %     else
    %         pos = 1;
    %     end        
    % 
    %     % samplePhase = mod(pos, sps) + floor(sps/2);   
    %     % rxDownsample = rcvAngles(samplePhase : sps: end);        
    %     % demodBits = zeros(size(rxDownsample));
    %     % idx = find(abs(diff(rxDownsample)) > pi/4);
    %     % demodBits(idx) = 1;
    %     % [indx,~] = prbdet(2*demodBits - 1);
    %     % 
    %     % if(isempty(indx))            
    %     %     pos = 1;
    %     % else            
    %     %     break;
    %     % end
    % end 

rcvAngles = unwrap(angle(signal));
samplePhase = mod(pos, sps) + floor(sps/2);
rxDownsample = rcvAngles(samplePhase: sps : end);
outBits = zeros(size(rxDownsample));
idx = find(abs(diff(rxDownsample)) > pi/4);
outBits(idx) = 1;

   if length (outBits) ~= 256
       outBits(256) = 5;
   end
end

