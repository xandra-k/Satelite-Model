% Регулятор коллизий
%StartPositions = [1 1]; % пересекаются на 100%
StartPositions = [1 1024]; % пересекаются на 50%
%StartPositions = [1 2049]; % не пересекаются

% Количество пакетов
NumPackets = 2;

% Массив допплеровских сдвигов
deltaFreq = [4000 -4000];

DopplerShift = 1;
% квадрат расстояния в зенитте 500 км, на горизонте 3000
% maxD = 1/(3000*1e3)^2;
% minD = 1/(500*1e3)^2;


 RecordAisSignal = zeros(RecordLen_smps*koef, 1);

for i = 1 : NumPackets
    % Информационные биты
        BitPacket(i,:) = randi ([0 1], PacketLen_bit, 1);
    % Формирование комплексной огибающей АИС-сигнала (GMSK модуляция)
        PacketLPE(i,:) = modGMSK(BitPacket(i,:)');

    if DopplerShift 
        % Сигнал (пакет), сдвинутый по частоте
        PacketLPE(i,:) = PacketLPE(i,:) .* exp(1i*2*pi*deltaFreq(i)/Fs*(0:PacketLen_smpls-1));
        % Ослабление сигнала в зависимости от частотного сдвига
        PacketLPE(i,:) = 1/abs(deltaFreq(i)) * PacketLPE(i,:);
    else
        PacketLPE(i,:) = PacketLPE(i,:);
    end
    

    % Добавление текущего сообщения в общую запись на случайное место
        EndPositions(i) = StartPositions(i) + PacketLen_smpls - 1;

    % Запишем текущую КО на места (start:end) с учётом уже разположеных там значений
        RecordAisSignal(StartPositions(i) : EndPositions(i)) = (RecordAisSignal(StartPositions(i) : EndPositions(i)) + PacketLPE(i,:).');
end

    RecordAisSignal = RecordAisSignal / max([abs(real(RecordAisSignal)); ...
        abs(imag(RecordAisSignal))]);

    CollisionsCounter = zeros(1,NumPackets);
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


BER1 = [];
BER2 = [];

%% ДЕМОДУЛЯЦИЯ 
%if isequal(Demodulator, 'Simple')
    for i = 1:NumPackets
        % Частотная синхронизация
            if DopplerShift 
                % Часть сигнала, предназнач. для демодуляции
                SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i)).* exp(1i*2*pi*-deltaFreq(i)/Fs*(0:PacketLen_smpls-1)).';
            else
                SignalToDemod(i,:) = RecordAisSignal(StartPositions(i):EndPositions(i));
            end

        % GMSK демодуляция
            DemodData =  demodGMSK(SignalToDemod(i,:).');
        % Запись в массив выходных бит
            OutBitPacket (i, :) = DemodData;
        % [errorCounter(i,:), ratio(i,:)] = biterr (BitPacket(i,:), OutBitPacket(i,:));
        
        % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
            %y = errorRate(BitPacket(i,:)', OutBitPacket(i,:)');
            y = biterr(BitPacket(i,1:end-16), OutBitPacket(i,17:end));
        % Сохраним BER
            BER1(i) = y/240;        
    end 

SignalToDemod = [];
OutBitPacket = [];
%elseif isequal(Demodulator, 'Zonal') 
    for i = 1:NumPackets
        % Частотная синхронизация
            if DopplerShift 
                % Часть сигнала, предназнач. для демодуляции
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
            % narrow_band1_leftpart = lowpass(SignalToDemod(i,:).*exp(1i*2*pi*fc1.*t), 4800, Fs, ImpulseResponse="fir");
            % z1 = narrow_band1_leftpart.*exp(1i*2*pi*-fc1.*t);
            % z1 = z1./ max([abs(real(z1)); abs(imag(z1))]);

            %narrow_band2_centered = filter (flip(ZonalFilterIR), 1, SignalToDemod(i,:));
            narrow_band2_centered = lowpass(SignalToDemod(i,:), 3000, Fs, ImpulseResponse="fir");
            z2 = narrow_band2_centered;
            %z2 = z2./ max([abs(real(z2)); abs(imag(z2))]);

            
            % narrow_band3_rightpart = filter (flip(ZonalFilterIR), 1, SignalToDemod(i,:).*exp(1i*2*pi*fc3.*t));
            % narrow_band3_rightpart = lowpass(SignalToDemod(i,:).*exp(1i*2*pi*fc3.*t), 4800, Fs, ImpulseResponse="fir"); 
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
                    DemodData =  demodGMSK(SubbandSig);
                % Запись в массив выходных бит
                    OutBitPacket (i, :) = DemodData.';
                % Вычисление ошибок с учётом задержки TraceBack Витерби Декодера
                    y = biterr(BitPacket(i,1:end-16), OutBitPacket(i,17:end));
                % Сохраним BER
                    BER2(ZonalDemodNum, i) = y/240; 

            end       
    end 