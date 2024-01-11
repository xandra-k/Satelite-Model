clear; clc;
%% Исходные параметры (Setup) 

% Ненастраевымые параметры
% Частота дискретизации, Гц
    Fs = 249600; % /2
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
    filterBT = 0.3;
% Количество отсчётов в символьных интервалах в GMSK модуляторе
    span = 3;

% Параметры модели    
% Индикатор сдвига Доплера
    DopplerShift = 1;
% Тип демодуляции
    Demodulator = 'Zonal'; % 'Zonal' 'Simple'

% Принятый сигнал АИС
% Длина одного АИС пакета в отсчётах 
    PacketLen_smpls = PacketLen_bit*sps;  
% Число пакетов АИС в записи, шт
    NumPackets = 100;  
% Длина записи в отсчётах
    RecordLen_smps = NumPackets*PacketLen_smpls;

% Вектор под запись демодулированных сигналов
    RecordAisSignal = zeros(RecordLen_smps*10, 1); 
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

demodGMSK = comm.GMSKDemodulator ('BandwidthTimeProduct', filterBT, ...
    'SamplesPerSymbol', sps, 'BitOutput',true,'PulseLength',...
    span, 'TracebackDepth', 16);

errorRate = comm.ErrorRate('ReceiveDelay', ...
                            demodGMSK.TracebackDepth);

% Параметры зонных фильтров
    Fs = Fs;                 % Sampling Frequency
    Fpass = 4800;            % Passband Frequency
    Fstop = 30000;            % Stopband Frequency
    Dpass = 0.057501127785;  % Passband 0.
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
        StartPositions(i) = randi([1, RecordLen_smps*10 - PacketLen_smpls + 1]);
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

% 
% %%Тут тестирую BER после фильтрации
% for i = 1:NumPackets
%     if DopplerShift 
%         % Сигнал (пакет), сдвинутый по частоте
%         SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i)).* exp(1i*2*pi*-deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';
%     else
%         SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
%     end
%     DemodData1 =  demodGMSK(SignalToDemod(i,:).');
%     %filtered2(i,:) = filter (flip(ZonalFilterIR), 1, SignalToDemod(i,:));
%     filtered(i,:) = lowpass(SignalToDemod(i,:), 4800, Fs, ImpulseResponse="fir"); % мб 4000 (найти границу)
%     SignalToDemod2(i,:) = filtered(i,:);
% 
%     %pwelch(filtered2(i,:), [], [], [], Fs, 'centered')
% 
%     DemodData2 =  demodGMSK(SignalToDemod2(i,:).');
%     %DemodData2 = circshift(DemodData2, -2);
%     % Запись в массив выходных бит
%     OutBitPacket1 (i, :) = DemodData1;
%     OutBitPacket2 (i, :) = DemodData2;
%     y1 = biterr(BitPacket(i,1:end-16), OutBitPacket1(i,17:end));
%     y2 = biterr(BitPacket(i,1:end-16), OutBitPacket2(i,17:end));
%     BER1(i) = y1/240;
%     BER2(i) = y2/240;
% end 
% ber1 = sum(BER1)/NumPackets
% ber2 = sum(BER2)/NumPackets

    % 
    % % pwelch(SignalToDemod(i,:), [], [], [], Fs, 'centered')
    % % hold on
    % pwelch(filtered(i,:), [], [], [], Fs, 'centered')
    % hold on




%%% ДЕМОДУЛЯЦИЯ %%%
if isequal(Demodulator, 'Simple')
    for i = 1:NumPackets
        % Частотная синхронизация
            if DopplerShift 
                % Сигнал (пакет), сдвинутый по частоте
                SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i)).* exp(1i*2*pi*-deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';
            else
                SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
            end

        % Часть сигнала, предназнач. для демодуляции
           % SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i)).* exp(1i*2*pi*-deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';
        % GMSK-демодулция i-го пакета
            DemodData =  demodGMSK(SignalToDemod(i,:).');
        % Запись в массив выходных бит
            OutBitPacket (i, :) = DemodData;
        % [errorCounter(i,:), ratio(i,:)] = biterr (BitPacket(i,:), OutBitPacket(i,:));
        
        % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
            y = errorRate(BitPacket(i,:)', OutBitPacket(i,:)');
            y = biterr(BitPacket(i,1:end-16), OutBitPacket(i,17:end));
        % Сохраним BER
            BER(i) = y(1)/240;        
    end 

elseif isequal(Demodulator, 'Zonal') 
    for i = 1:NumPackets
        % Частотная синхронизация
            if DopplerShift 
                % Сигнал (пакет), сдвинутый по частоте
                SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i)).* exp(1i*2*pi*-deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';
            else
                SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
            end
        % ....
            t = (0:length(SignalToDemod)-1)/Fs;
        % Центральные частоты поддиапазонов (пока не параметризированы )
            fc1 = -3840; fc2 = 0; fc3 = 3840;
        % Выделяем три поддиапазона демодулятора (временная область, exp - сдвиг на fc)
            
            %narrow_band1_leftpart = filter (flip(ZonalFilterIR), 1, SignalToDemod(i,:).*exp(1i*2*pi*fc1.*t));
            narrow_band1_leftpart = lowpass(SignalToDemod(i,:).*exp(1i*2*pi*fc1.*t), 4800, Fs, ImpulseResponse="fir");
            z1 = narrow_band1_leftpart.*exp(1i*2*pi*0-fc1.*t);
            z1 = z1./ max([abs(real(z1)); abs(imag(z1))]);

            %narrow_band2_centered = filter (flip(ZonalFilterIR), 1, SignalToDemod(i,:));
            narrow_band2_centered = lowpass(SignalToDemod(i,:), 4800, Fs, ImpulseResponse="fir");
            z2 = narrow_band2_centered;
            z2 = z2./ max([abs(real(z2)); abs(imag(z2))]);

            
            % narrow_band3_rightpart = filter (flip(ZonalFilterIR), 1, SignalToDemod(i,:).*exp(1i*2*pi*fc3.*t));
            narrow_band3_rightpart = lowpass(SignalToDemod(i,:).*exp(1i*2*pi*fc3.*t), 4800, Fs, ImpulseResponse="fir"); 
            z3 = narrow_band3_rightpart.*exp(1i*2*pi*0-fc3.*t);
            z3 = z3./ max([abs(real(z3)); abs(imag(z3))]);


        % Массив сигналов соответвующих поддиапазонов
            ZonalBandSignals = [z1', z2', z3'];


        % Цикл по зонным демодуляторам
            for ZonalDemodNum = 1:3 % для z2 полностью инвертированны
                % Сигнал текущего поддиапазона
                    SubbandSig = ZonalBandSignals(:, ZonalDemodNum);
                % GMSK-демодулция i-го пакета, ZonalDemodNum-го демодулятора
                    DemodData =  demodGMSK(SubbandSig);
                % Запись в массив выходных бит
                    OutBitPacket (i, :) = ((-1)*DemodData' +1);
                % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
                    y = biterr(BitPacket(i,1:end-16), OutBitPacket(i,17:end));
                % Сохраним BER
                    BER(ZonalDemodNum, i) = y(1)/240; 

            end       
    end 
end