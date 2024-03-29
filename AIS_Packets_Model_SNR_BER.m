clear; clc; close all;
%% Исходные параметры (Setup) 

% Ненастраевымые параметры
% Частота дискретизации, Гц
    %Fs = 249600; % /2
    Fs = 9600*6;
% Длительность одного пакета АИС в битах (26.67 мс)
    PacketLen_bit = 256*100; 
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

% Параметры модели    
% Индикатор сдвига Доплера
    DopplerShift = 1;
% Тип демодуляции
    Demodulator = 'Zonal'; % 'Zonal' 'Simple'
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

SNRs = 20:2:20;
%% Модель 
for n = 1:length (SNRs) % цикл по осш
    RecordAisSignal = zeros(RecordLen_smps*koef, 1);    
    for i = 1 : NumPackets
        % Формирование принятого пакета
        % Информационные биты
            BitPacket(i,:) = randi ([0 1], PacketLen_bit, 1);
        % Формирование комплексной огибающей АИС-сигнала (GMSK модуляция)
        %PacketLPE = GMSKmod(BitPacket(i,:), span, sps, filterBT);
            PacketLPE(i,:) = modGMSK(BitPacket(i,:)');
        
        % Допплеровский сдвиг (пока только сдвиг частоты)
        % Нижняя граница доплеровского сдвига, Гц
            dF_Bottom = -3000;
        % Верхняя граница доплеровского сдвига, Гц 
            dF_Upper = 3000; 
        % Случайная величина частотного сдвига
            %deltaFreq(i) = randi([dF_Bottom, dF_Upper]); 
            deltaFreq = [-2200, +2200];
    
        if DopplerShift 
            % Сигнал (пакет), сдвинутый по частоте
            PacketLPE(i,:) = PacketLPE(i,:) .* exp(1i*2*pi*deltaFreq(i)/Fs*(0:PacketLen_smpls-1));
            % Ослабление сигнала в зависимости от частотного сдвига
            %PacketLPE(i,:) = 1/abs(deltaFreq(i)) * PacketLPE(i,:);
            if i == 1
                PacketLPE(i,:) = PacketLPE(i,:);
            else
                PacketLPE(i,:) = PacketLPE(i,:);
            end

        else
            PacketLPE(i,:) = PacketLPE(i,:);
        end

       % Carrier_to_Interference = 10*log10(PacketLPE(1,:)/sum(PacketLPE(2:end,:))); % надо посчитать в дб 
    
        % Случайный поворот фазы
        % Phi = -pi/2 + pi*rand; % от -pi/2 до pi/2
        % PacketLPE(i,:) = PacketLPE(i,:) * exp(1i*Phi);
        
        % Случайное ослабление сигнала
        % PacketLPE(i,:) = randi([1000, 30000]) .* PacketLPE(i,:);
    
    
        % Добавление текущего сообщения в общую запись на случайное место
        %if CollisionsOn
            %StartPositions(i) = randi([1, RecordLen_smps*koef - PacketLen_smpls + 1]);
            StartPositions = [1 1]; % пересекаются на 100%
            %StartPositions = [1 0.5*PacketLen_smpls+1]; % пересекаются на 50%
            % StartPositions = [1, PacketLen_smpls+1]; % не пересекаются
        % else
        %     if i == 1
        %         StartPositions(i) = 1;
        %     else
        %         StartPositions(i) = StartPositions(i-1) + PacketLen_smpls * (1 + PauseLen);
        %     end
        % end
    
            EndPositions(i) = StartPositions(i) + PacketLen_smpls - 1;
        % Запишем текущую КО на места (start:end) с учётом уже разположеных там значений
            RecordAisSignal(StartPositions(i) : EndPositions(i)) = (RecordAisSignal(StartPositions(i) : EndPositions(i)) + PacketLPE(i,:).');
    
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
    
    %% КАНАЛ АБГШ
        if Channel 
            % Добавляем шумовые отсчёты к сигналу
            RecordAisSignal = awgn(RecordAisSignal, SNRs(n), "measured"); 
        else
            RecordAisSignal = RecordAisSignal;
        end
    
    %% ДЕМОДУЛЯЦИЯ 
    %if isequal(Demodulator, 'Simple')
    SignalToDemod2 = [];
    SignalToDemod1 = [];
        for i = 1:NumPackets
            % Частотная синхронизация
                if DopplerShift 
                    % Часть сигнала, предназнач. для демодуляции
                    SignalToDemod1(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i)).* exp(1i*2*pi*-deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';
                else
                    SignalToDemod1(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
                end
    
            % GMSK демодуляция
                DemodData1 =  demodGMSK(SignalToDemod1(i,:).');
            % Запись в массив выходных бит
                OutBitPacket1 (i, :) = DemodData1;
            % [errorCounter(i,:), ratio(i,:)] = biterr (BitPacket(i,:), OutBitPacket(i,:));
            
            % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
                %y = errorRate(BitPacket(i,:)', OutBitPacket1(i,:)');
                y = biterr(BitPacket(i,1:end-16), OutBitPacket1(i,17:end));
            % Сохраним BER
                BER1(i, n) = y/(256*100 - 16);        
        end 
    
    SignalToDemod = [];
    %elseif isequal(Demodulator, 'Zonal') 
        for i = 1:NumPackets
            % Частотная синхронизация
                if DopplerShift 
                    % Часть сигнала, предназнач. для демодуляции
                    SignalToDemod2(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i)).* exp(1i*2*pi*-deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';
                else
                    SignalToDemod2(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
                end
    
            % ....
                t = (0:length(SignalToDemod2)-1)/Fs;
            % Центральные частоты поддиапазонов (пока не параметризированы )
                fc1 = -3840; fc2 = 0; fc3 = 3840;
            % Выделяем три поддиапазона демодулятора (временная область, exp - сдвиг на fc)
                
                %narrow_band1_leftpart = filter (flip(ZonalFilterIR), 1, SignalToDemod2(i,:).*exp(1i*2*pi*fc1.*t));
                % narrow_band1_leftpart = lowpass(SignalToDemod2(i,:).*exp(1i*2*pi*fc1.*t), 4800, Fs, ImpulseResponse="fir");
                % z1 = narrow_band1_leftpart.*exp(1i*2*pi*-fc1.*t);
                % z1 = z1./ max([abs(real(z1)); abs(imag(z1))]);
    
                %narrow_band2_centered = filter (flip(ZonalFilterIR), 1, SignalToDemod(i,:));
                narrow_band2_centered = lowpass(SignalToDemod2(i,:), 200, Fs, ImpulseResponse="fir");
                z2 = narrow_band2_centered;
                %z2 = z2 / max([abs(real(z2)); abs(imag(z2))]); % нормировка
    
                
                % narrow_band3_rightpart = filter (flip(ZonalFilterIR), 1, SignalToDemod2(i,:).*exp(1i*2*pi*fc3.*t));
                % narrow_band3_rightpart = lowpass(SignalToDemo2d2(i,:).*exp(1i*2*pi*fc3.*t), 4800, Fs, ImpulseResponse="fir"); 
                % z3 = narrow_band3_rightpart.*exp(1i*2*pi*-fc3.*t);
                % z3 = z3./ max([abs(real(z3)); abs(imag(z3))]);
    
    
            % Массив сигналов соответвующих поддиапазонов
                %ZonalBandSignals = [z1', z2', z3'];
                ZonalBandSignals = [z2.'];     
    
            % Цикл по зонным демодуляторам
                for ZonalDemodNum = 1:1 % для z2 полностью инвертированны
                    % Сигнал текущего поддиапазона
                        SubbandSig = ZonalBandSignals(:, ZonalDemodNum);
                    % GMSK-демодулция i-го пакета, ZonalDemodNum-го демодулятора
                        DemodData2 =  demodGMSK(SubbandSig);
                    % Запись в массив выходных бит
                        OutBitPacket2 (i, :) = DemodData2.';
                    % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
                        y = biterr(BitPacket(i,1:end-16), OutBitPacket2(i,17:end));
                    % Сохраним BER
                        %BER2(ZonalDemodNum, i) = y/240; 
                        BER2(i, n) = y/(256*100-16); 
                end       
        end 
    %end
end % цикл по SNR

%% График
figure()
semilogy(SNRs, BER1(1,:))
hold on
semilogy(SNRs, BER2(1,:))
grid on
legend ('Simple GMSK demod', 'Zonal demod')