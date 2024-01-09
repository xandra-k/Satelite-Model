 clear; clc;
%% Исходные параметры (Setup) 

% Ненастраевымые параметры
% Частота дискретизации, Гц
    Fs = 249600;
% Длительность одного пакета АИС в битах (26.67 мс)
    PacketLen_bit = 256; 
% Символьная скорость передачи, б/c
    SymbolRate = 9600; 
% Символьный интервал
    SymbolInterval = 1/SymbolRate;  
% Количество отсчётов на один символ (бит, т.к ФМ-2). Формула round(Fs / SymbolRate)
    sps = round(Fs / SymbolRate);
    %sps = 26; 
% Коэффициент фильтра в GMSK модуляторе
    filterBT = 1.0;
% Количество отсчётов в символьных интервалах в GMSK модуляторе
    span = 1;

% Параметры модели    
% Индикатор сдвига Доплера
    DopplerShift = 0;
% Тип демодуляции
    Demodulator = 'Zonal'; % 'Zonal' 'Simple'

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

% Модулятор, демодулятор, счётчик ошибок
modGMSK = comm.GMSKModulator('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitInput',true,'PulseLength',...
    span);

demodGMSK = comm.GMSKDemodulator ('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitOutput',true,'PulseLength',...
    span, 'TracebackDepth', 16);

errorRate = comm.ErrorRate('ReceiveDelay', ...
                            demodGMSK.TracebackDepth);

% Параметры зонных фильтров
    Fs = Fs;                 % Sampling Frequency
    Fpass = 4800;            % Passband Frequency
    Fstop = 9000;            % Stopband Frequency
    Dpass = 0.057501127785;  % Passband Ripple
    Dstop = 0.0001;          % Stopband Attenuation
    dens  = 20;              % Density Factor
    [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
    ZonalFilterIR  = firpm(N, Fo, Ao, W, {dens}); % Импульсная характеристика
% Полоса фильтра, Гц
    bandWidth = 1/SymbolInterval; 

%% Модель 
for i = 1 : NumPackets
    % Формирование принятого пакета
    % Информационные биты
        BitPacket(i,:) = randi ([0 1], PacketLen_bit, 1);
    % Формирование комплексной огибающей АИС-сигнала (GMSK модуляция)
    %PacketLPE = GMSKmod(BitPacket(i,:), span, sps, filterBT);
        PacketLPE(i,:) = modGMSK(BitPacket(i,:)');
    
    % Допплеровский сдвиг (пока только сдвиг частоты)
    % Нижняя граница доплеровского сдвига, Гц
        dF_Bottom = -4000;
    % Верхняя граница доплеровского сдвига, Гц 
        dF_Upper = 4000; 
    % Случайная величина частотного сдвига
        deltaFreq(i) = randi([dF_Bottom, dF_Upper]); 

    if DopplerShift 
        % Сигнал (пакет), сдвинутый по частоте
        PacketLPE(i,:) = PacketLPE(i,:) .* exp(1i*2*pi*deltaFreq(i)/Fs*(0:PacketLen_smpls-1));
    else
        PacketLPE(i,:) = PacketLPE(i,:);
    end

    % Случайный поворот фазы
    % Phi = -pi/2 + pi*rand; % от -pi/2 до pi/2
    % PacketLPE(i,:) = PacketLPE(i,:) * exp(1i*Phi);
    
    % Случайное ослабление сигнала
    % PacketLPE(i,:) = randi([1000, 30000]) .* PacketLPE(i,:);


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
        RecordAisSignal(StartPositions(i) : EndPositions(i)) = (RecordAisSignal(StartPositions(i) : EndPositions(i)) + PacketLPE(i,:)');

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
%disp (['Общее число коллизий в каждом пакете: ', num2str(CollisionsCounter')])

DemodData = zeros (NumPackets, PacketLen_bit);
%%% ДЕМОДУЛЯЦИЯ %%%
if isequal(Demodulator, 'Simple')
    for i = 1:NumPackets
        % Часть сигнала, предназнач. для демодуляции
            SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
        % GMSK-демодулция i-го пакета
            DemodData =  demodGMSK(SignalToDemod(i,:)');
        % Запись в массив выходных бит
            OutBitPacket (i, :) = DemodData';
        %[errorCounter(i,:), ratio(i,:)] = biterr (BitPacket(i,:), OutBitPacket(i,:));
        
        % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
            y = errorRate(BitPacket(i,:)', OutBitPacket(i,:)');
        % Сохраним BER
            BER(i) = y(1);        
    end 
elseif isequal(Demodulator, 'Zonal') 
    for i = 1:NumPackets
        % Часть сигнала, предназнач. для демодуляции
            SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
        % ....
        t = (0:length(SignalToDemod)-1)/Fs;
        % Центральные частоты поддиапазонов (пока не параметризированы ~0.8Rb)
            fc1 = -3840.1; fc2 = 0; fc3 = 3840.1;
        % Выделяем три поддиапазона демодулятора (временная область, exp - сдвиг на fc)
            narrow_band1_leftpart = filter (flip(ZonalFilterIR), 1, SignalToDemod.*exp(-1i*2*pi*fc1.*t));
            z1 = narrow_band1_leftpart.*exp(1i*2*pi*fc1.*(0:length(narrow_band1_leftpart)-1)/Fs);
        
            narrow_band2_centered = filter (flip(ZonalFilterIR), 1, SignalToDemod);
            z2 = narrow_band2_centered;
            
            narrow_band3_rightpart = filter (flip(ZonalFilterIR), 1, SignalToDemod.*exp(-1i*2*pi*fc3.*t));
            z3 = narrow_band3_rightpart.*exp(1i*2*pi*fc3.*(0:length(narrow_band3_rightpart)-1)/Fs);

        % Массив сигналов соответвующих поддиапазонов
            ZonalBandSignals = [z1', z2', z3'];
        % Цикл по зонным демодуляторам
            for ZonalDemodNum = 1:3
                % Сигнал текущего поддиапазона
                    SubbandSig = ZonalBandSignals(:, ZonalDemodNum);
                % GMSK-демодулция i-го пакета, ZonalDemodNum-го демодулятора
                    DemodData =  demodGMSK(SubbandSig);
                % Запись в массив выходных бит
                    OutBitPacket (i, :) = DemodData';
                % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
                    y = errorRate(BitPacket(i,:)', OutBitPacket(i,:)');
                % Сохраним BER
                    BER(ZonalDemodNum, i) = y(1); 

            end       
    end 
end