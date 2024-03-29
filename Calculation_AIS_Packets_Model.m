clear; clc; close all;
%% Исходные параметры (Setup) 
tic
load ("DopDist.mat")
% Ненастраевымые параметры
% Частота дискретизации, Гц
    Fs = 249600; % /2
    %Fs = 9600*8;
% Длительность одного пакета АИС в битах (26.67 мс)
    lenKoef = 1;
    PacketLen_bit = 256*lenKoef; 
% Символьная скорость передачи, б/c
    SymbolRate = 9600; 
% Символьный интервал
    SymbolInterval = 1/SymbolRate;  
% Количество отсчётов на один символ (бит, т.к ФМ-2). Формула round(Fs / SymbolRate)
    sps = round(Fs / SymbolRate);
    %sps = 26; 
% Коэффициент фильтра в GMSK модуляторе
    filterBT = 0.3;
% Количество отсчётов в символьных интервалах в GMSK модуляторе
    span = 3;
% Мощность передатчика 
    Ptx = 2; % Вт
    Ptx_dB = 10*log10(Ptx); % дБ
    Ptx_dBm = 10*log10(Ptx) + 30; % дБм
% Длина волны
    f = 161.971*10^6;
    lamda = 1/f * 3 * 10^8;

% Параметры модели    
% Индикатор сдвига Доплера
    DopplerShift = 1;
% Тип демодуляции
    % Demodulator = 'Zonal'; % 'Zonal' 'Simple'
% Канал АБГШ
    Channel = 1;
% Фиксироаванный ОСШ 
    SNR = 15; 

% Принятый сигнал АИС
% Длина одного АИС пакета в отсчётах 
    PacketLen_smpls = PacketLen_bit*sps;  
% Число пакетов АИС в записи, шт
    NumPackets = 2;  
% Длина записи в отсчётах
    RecordLen_smps = NumPackets*PacketLen_smpls;
% Расширяющий коэффициент записи (практически влияет на плотность коллизий)
    koef = 1;

% Вектор под запись демодулированных сигналов
    RecordAisSignal = zeros(RecordLen_smps*koef, 1); 
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


% % Принятый сигнал
%     bitsAisSignal = 2*randi ([0 1], 1, RecordLen_smps) - 1;

% Модулятор, демодулятор, счётчик ошибок
modGMSK = comm.GMSKModulator('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitInput',true,'PulseLength',...
    span);

demodGMSK = comm.GMSKDemodulator('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitOutput',true,'PulseLength',...
    span, 'TracebackDepth', 16);

errorRate = comm.ErrorRate('ReceiveDelay', ...
                            demodGMSK.TracebackDepth);

    %% МОДЕЛЬ 

% Допплеровский сдвиг (пока только сдвиг частоты)
% Вектор всевозможных доплеровских сдвигов (-3.8150e3 3.8150e3)
    dopFreqsVect = DopDist (1,:);
    %dopFreqsMatrix = zeros (length(dopFreqsVect), length(dopFreqsVect));
% Вектор расстояний от передатчика до приёмника (от судна до спутника)
    distanceVect = DopDist (2,:);

cntTX = 1;    
% Цикл по частотным сдвигам
for cnt1 = 1:3
    for cnt2 = 1:3 %1:length(dopFreqsVect)
    % Текущий допплеровский сдвиг
        deltaFreq = [dopFreqsVect(cnt1) dopFreqsVect(cnt2)]; 
        distance = [distanceVect(cnt1) distanceVect(cnt2)];

        % Цикл по временной коллизии
        for packet_overlap = 0:0.1:0.3 %0:0.1:1
        % Обнуляем вектор для записи    
            RecordAisSignal = zeros(RecordLen_smps*koef, 1); 

        % Передаваемые биты
            BitPackets = randi ([0 1], PacketLen_bit, 2);

        % GMSK модуляция двух пакетов
            PacketLPE1 = modGMSK(BitPackets(:, 1));
            PacketLPE2 = modGMSK(BitPackets(:, 2));
            PacketLPE = [PacketLPE1.'; PacketLPE2.'];

        % Затухание в свободном пространстве (Free-space path loss)
            Ls_dB1 = 10*log10((lamda/(4*pi*distance(1)))^2);  % в дБ
            Ls_dB2 = 10*log10((lamda/(4*pi*distance(2)))^2);  % в дБ
    
        % Мощность принятого сигнала (c учётом мощности передатчика)
            Prx_dBm(1,:) = Ptx_dBm + Ls_dB1;
            Prx_dBm(2,:) = Ptx_dBm + Ls_dB2;

        % Разница принятых мощностей в дБ
            %delta = Prx1_dBm - Prx2_dBm;
            delta = Prx_dBm(1,:) - Prx_dBm(2,:);

        % Энергия 1го сигнала
            E1 = sum(abs(PacketLPE(1,:).^2)); % разы
            E1_db = 10*log10(E1); % дБ

        % Энергия 2го сигнала
            E2_db = E1_db - delta;
            E2 =  10^(E2_db/10); 

        % Амплитудный множитель
            x = sqrt(E2/E1);

        % Перерасчёт пакетов с учётом разной энергии в зависимости от
        % расстояния м/у TX и RX
            PacketLPE(1,:) = PacketLPE(1,:);
            PacketLPE(2,:) = x*PacketLPE(2,:);

        % Сигнал (пакет), сдвинутый по частоте
            PacketLPE(1,:) = PacketLPE(1,:) .* exp(1i*2*pi*deltaFreq(1)/Fs*(0:PacketLen_smpls-1));
            PacketLPE(2,:) = PacketLPE(2,:) .* exp(1i*2*pi*deltaFreq(2)/Fs*(0:PacketLen_smpls-1));

        % Цикл по расстановке пакетов
        % Стартовая позиция    
            StartPositions = [1 packet_overlap*PacketLen_smpls+1]; % не пересекаются на packet_overlap процентов
            for i = 1 : NumPackets
                    EndPositions(i) = StartPositions(i) + PacketLen_smpls - 1;
                % Запишем текущую КО на места (start:end) с учётом уже разположеных там значений
                    RecordAisSignal(StartPositions(i) : EndPositions(i)) = (RecordAisSignal(StartPositions(i) : EndPositions(i)) + PacketLPE(i,:).');
            end

        % Нормировка результата
            RecordAisSignal = RecordAisSignal / max([abs(real(RecordAisSignal)); ...
                abs(imag(RecordAisSignal))]);           

        % Запись в структуру    
            TotalRecord.BitPacket(cntTX, :) = BitPackets(:, 1);
            TotalRecord.Signal(cntTX, :) =  RecordAisSignal;
            TotalRecord.freq1(cntTX,:) = dopFreqsVect(cnt1);
            TotalRecord.freq2(cntTX,:) = dopFreqsVect(cnt2);
            TotalRecord.timeCollision(cntTX,:) = packet_overlap;
            TotalRecord.StartPositions(cntTX,:) = StartPositions;
            TotalRecord.EndPositions(cntTX,:) = EndPositions;

            
        % Cчётчик
            cntTX = cntTX + 1;
        end
    end
end
    
    
    %% КАНАЛ АБГШ

    % Добавляем шумовые отсчёты к сигналу
        TotalRecord.Signal = awgn(TotalRecord.Signal, SNR, "measured"); 
    
    %% ДЕМОДУЛЯЦИЯ 

for cntRX = 1:cntTX-1
    % Исходные биты для проверки
        BitPacket = TotalRecord.BitPacket(cntRX, :);
    % Текущая запись для демодуляции
        RecordAisSignal = TotalRecord.Signal(cntRX, :).';
    % Частотный сдвиг
        deltaFreq = TotalRecord.freq1(cntRX,:); 
    % Коллизии
        timeCollision = TotalRecord.timeCollision(cntRX,:);
    % Начало и конец 1го и 2го пакетов
        StartPositions = TotalRecord.StartPositions(cntRX,:);
        EndPositions = TotalRecord.EndPositions(cntRX,:);
    % Сигнал, предназначенный для демодуляции
        SignalToDemod = RecordAisSignal(StartPositions(1):EndPositions(1)).* exp(1i*2*pi*-deltaFreq/Fs*(0:PacketLen_smpls-1)).';
    
    % Классический демодулятор
        DemodData1 =  demodGMSK(SignalToDemod);
        OutBitPacket1  = DemodData1; 
        y = biterr(BitPacket(1:end-16), OutBitPacket1(17:end).');
        BER1(cntRX) = y/(256*lenKoef - 16); 

    % Зонный демодулятор
        t = (0:length(SignalToDemod)-1)/Fs;
        narrow_band2_centered = lowpass(SignalToDemod, 700, Fs, ImpulseResponse="fir");
        SubbandSig = [narrow_band2_centered];        
        DemodData2 =  demodGMSK(SubbandSig);
        OutBitPacket2  = DemodData2.';
        y = biterr(BitPacket(1:end-16), OutBitPacket2(17:end));
        BER2(cntRX) = y/(256*lenKoef - 16); 


    % Cчётчик
        cntRX = cntRX + 1;        
end
toc