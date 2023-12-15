clc; clear; close all;

%% Добавление LEO-спутника
% Создаём спутниковый сценарий
startTime = datetime(2023,5,5,0,0,0);
stopTime = startTime + hours(6);
sampleTime = 10; % секунд
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Добавляем спутник из tle файла
tleFile = "cubesat.tle";
sat = satellite(sc, tleFile, Name = "Cubesat");
%% Одна базовая станция
targetGs = groundStation(sc, Name="Location to Catch", ...
    Latitude=-16.943443,Longitude=177.835693);    

%% Поле видимости антенны к спутнику
% Common approach is to fix the sensor on a gimbal and orient the sensor 
% by manuevering the gimbal, rather than the spacecraft body itself
g = gimbal(sat);

% Управление направлением спутника/сенсора по отношению к цели 
% pointAt(g, gs); % при помощи gimbal'a
% pointAt(sat, "nadir"); % всегда надир

% Добавляем сенсор к спутнику
SensorAngle = 100;
camSensor = conicalSensor(g, "MaxViewAngle", SensorAngle, Name = "Satellie Sensor");

%% Создаём связь между камерой и базовой станцией
% 
% % БС в поле видимости спутника, но связь с камерой не установлена - красный
% acFalse = access(sat, gs);
% acFalse.LineColor = 'red';
% 
% % Связь между камерой и БС установлена - зелёный
% acTrue = access(camSensor, gs);
% acTrue.LineColor = 'green';

% accessStatusTrue = accessStatus(acTrue);
% accessStatusFalse = accessStatus(acFalse);
% 
% % Определим отсчёты времени в которые: 
% % БС в поле видимости спутника, но связь с камерой не установлена
% indFoV = find (accessStatusFalse ==1);
% % Связь между камерой и БС установлена
% indConnected = find (accessStatusTrue ==1);
% % Индексы массива indFoV характеризующие одновременное попадание в зону видимости и установку связи 
% indConnected_FoV = find(ismember(indFoV.',indConnected.','rows'));


%% Визуализация сценария 
v = satelliteScenarioViewer(sc,"ShowDetails", true);
% Включить визуализацию поля видимости спутника
fov = fieldOfView(camSensor);
%play(sc)

%% Добавляем сетку базовых-станций (кораблей)
% target location: Latitude=-16.943443,Longitude=177.835693

latlim = [-30 -1];
lonlim = [162 192];
spacingInLatLon = 1; % degrees

numShips = 10;

groundStationsLocationslats = randi([latlim], 1, numShips);
groundStationsLocationslons = randi([lonlim], 1, numShips);


gslat = groundStationsLocationslats;
gslon = groundStationsLocationslons;

groundStations = groundStation(sc, gslat, gslon, Name="Ship" + string(1:numShips)');
%% Создаём связь между камерой и базовой станцией
% БС в поле видимости спутника, но связь с камерой не установлена - красный
acStationsFalse = access(sat, [groundStations targetGs]);
acStationsFalse.LineColor = 'red';

% Связь между камерой и БС установлена - зелёный
acStationsTrue = access(camSensor, [groundStations targetGs]);
acStationsTrue.LineColor = 'green';

%% Определяем отсчёты времени, когда связь м/у спутником и объектом устанволена
intvlsTrue = accessIntervals(acStationsTrue)
intvlsFalse = accessIntervals(acStationsFalse)

%% Add Transmitter to LEO  Satellite

antennaf = 1090e6; % in Hz
satelliteAntenna = arrayConfig("Size",[1 1]); % Add code comment    
satelliteTransmitter = transmitter( ...
    sat, ...
    Antenna=satelliteAntenna, ...
    Power = 10*log10(125),...
    MountingAngles= [0,0,0]);
pattern(satelliteTransmitter,Size=50000); 

% Play the scenario.
play(sc);

%% Add Receiver to ground object
GroundStationsAntenna = arrayConfig("Size",[1 1]); % Create an isotropic antenna element
GroundStationsReceiver = receiver(targetGs, ...
    Antenna=GroundStationsAntenna);
pattern(GroundStationsReceiver,antennaf, Size=50000);
%% Определяем расстояние до спутнка
a = [ intvlsTrue.StartTime(3) intvlsTrue.EndTime(3)]
for i = 1:length(a)
    [targetGs_data(i).azimuth, targetGs_data(i).elevation, targetGs_data(i).range, targetGs_data(i).timeOut] = aer(sat,targetGs, a(i));
    ship1(i).Distance = distance (targetGs.Latitude, targetGs.targetGs.Longitude,...
        satelliteTransmitter);
end


