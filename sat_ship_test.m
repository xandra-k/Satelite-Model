clc; clear; close all;

%% Константы 
EarthRadius = physconst('EarthRadius');

%% Добавление LEO-спутника
% Создаём спутниковый сценарий
startTime = datetime(2024,3,24,0,16,56);
stopTime = startTime + hours(6);
sampleTime = 10; % секунд
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Добавляем спутник из tle файла
tleFile = "PolytechUniverse1.tle";
sat = satellite(sc, tleFile, Name = "Cubesat");

% Определим широту, долготу, высоту спутника в начальный момент 
lla = states(sat,startTime,"CoordinateFrame","geographic");
%% Одна базовая станция
targetGs = groundStation(sc, Name="Location to Catch", ...
    Latitude=lla(1),Longitude=lla(2));    

%% Поле видимости антенны к спутнику
% Common approach is to fix the sensor on a gimbal and orient the sensor 
% by manuevering the gimbal, rather than the spacecraft body itself
g = gimbal(sat);

% Управление направлением спутника/сенсора по отношению к цели 
% pointAt(g, gs); % при помощи gimbal'a
% pointAt(sat, "nadir"); % всегда надир

% Добавляем сенсор к спутнику
SensorAngle = 90;
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
% target (-37.3397603503485, -14.5454154598156)
latlim = [-38.864670 -36.571163]; 
lonlim = [-15.296714 -14.480875];
spacingInLatLon = 0.01; % degrees

proj = projcrs(4087); %??
spacingInXY = deg2km(spacingInLatLon)*1000; % meters
[xlimits,ylimits] = projfwd(proj,latlim,lonlim);
R = maprefpostings(xlimits,ylimits,spacingInXY,spacingInXY);
[X,Y] = worldGrid(R);
[gridlat,gridlon] = projinv(proj,X,Y);
AtlanticOcean = readgeotable("C:\Users\AKuznecova\Desktop\Github\Satelite-Model\AtlanticOcean\AtlanticOcean.shp");
NorthAtlanticOcean = AtlanticOcean(AtlanticOcean.name == "North Atlantic Ocean",:);
T = geotable2table(AtlanticOcean,["Latitude","Longitude"]);
[landlat,landlon] = polyjoin(T.Latitude,T.Longitude);
[landlat,landlon] = maptrimp(landlat,landlon,latlim,lonlim);

bufwidth = 1;
[landlatb,landlonb] = bufferm(landlat,landlon,bufwidth,"outPlusInterior");
australiab = geopolyshape(landlatb,landlonb);
gridpts = geopointshape(gridlat,gridlon);
inregion = isinterior(australiab,gridpts);

gslat = gridlat(inregion);
gslon = gridlon(inregion);
gs = groundStation(sc,gslat,gslon);
fq = 162e6; % Hz
txpower = 20; % dBW
antennaType = "Monopole"; % "Monopole";

effectiveLine = sqrt((500e3 + EarthRadius)^2 - EarthRadius^2); % касательная к земле по углу видимости 
halfBeamWidth = acosd(effectiveLine/EarthRadius);


% numShips = 10;
% 
% groundStationsLocationslats = randi([latlim], 1, numShips);
% groundStationsLocationslons = randi([lonlim], 1, numShips);
% 
% 
% gslat = groundStationsLocationslats;
% gslon = groundStationsLocationslons;
% 
% groundStations = groundStation(sc, gslat, gslon, Name="Ship" + string(1:numShips)');
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
antennaType = "Monopole";
if antennaType == "Monopole"
    lambda = physconst('lightspeed')/fq; % meters
    monopoleHeight = lambda/2;
    monopoleWidth = 1e-2;
    mpl = monopole('Height', monopoleHeight, 'Width', monopoleWidth); 
    tx = transmitter(iridiumSatellites, ...
        Frequency=fq, ...
        Power=txpower, ...
        Antenna=mpl); 
end

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


