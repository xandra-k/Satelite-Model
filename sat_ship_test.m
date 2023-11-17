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
gs = groundStation(sc, Name="Location to Catch", ...
    Latitude=-16.943443,Longitude=177.835693);    

%% Поле видимости антенны к спутнику
% Common approach is to fix the sensor on a gimbal and orient the sensor 
% by manuevering the gimbal, rather than the spacecraft body itself
g = gimbal(sat);

% Управление направлением спутника/сенсора по отношению к цели 
% pointAt(g, gs); % при помощи gimbal'a
% pointAt(sat, "nadir"); % всегда надир

% Добавляем сенсор к спутнику
SensorAngle = 180;
camSensor = conicalSensor(g, "MaxViewAngle", SensorAngle, Name = "Satellie Sensor");

%% Создаём связь между камерой и базовой станцией

% БС в поле видимости спутника, но связь с камерой не установлена - красный
acFalse = access(sat, gs);
acFalse.LineColor = 'red';

% Связь между камерой и БС установлена - зелёный
acTrue = access(camSensor, gs);
acTrue.LineColor = 'green';

s = accessStatus(acFalse);

%% Визуализация сценария 
v = satelliteScenarioViewer(sc,"ShowDetails", true);
fov = fieldOfView(camSensor)
%play(sc)
%% Определяем отсчёты времени, когда связь м/у спутником и объектом устанволена
% intvlsTrue = accessIntervals(acTrue)
% intvlsFalse = accessIntervals(acFalse)

%% Добавляем cсетку базовых-станций (кораблей)
% Grid of Ground Stations Covering Fiji region
% target location: Latitude=-16.943443,Longitude=177.835693

latlim = [-30 -1];
lonlim = [162 192];
spacingInLatLon = 1; % degrees

groundStationsLocationslats = randi([latlim], 1, 1);
groundStationsLocationslons = randi([lonlim], 1, 1);


% % Create a projected coordinate reference system (CRS) that is appropriate for Fiji
% proj = projcrs(3139-42); 
% 
% % Calculate the projected X-Y map coordinates for the data grid. 
% spacingInXY = deg2km(spacingInLatLon)*1000; % meters
% [xlim,ylim] = projfwd(proj,latlim,lonlim);
% R = maprefpostings(xlim,ylim,spacingInXY,spacingInXY);
% [X,Y] = worldGrid(R);
% 
% % Transform the grid locations to latitude-longitude coordinates.
% [gridlat,gridlon] = projinv(proj,X,Y);
% landareas = readgeotable("landareas.shp");
% surface = landareas(landareas.Name == "Fiji",:);
% 
% T = geotable2table(surface,["Latitude","Longitude"]);
% [landlat,landlon] = polyjoin(T.Latitude,T.Longitude);
% 
% bufwidth = 1;
% [landlatb,landlonb] = bufferm(landlat,landlon,bufwidth,"outPlusInterior");
% australiab = geopolyshape(landlatb,landlonb);
% 
% gridpts = geopointshape(gridlat,gridlon);
% inregion = isinterior(australiab,gridpts);
% 
% gslat = gridlat(inregion);
% gslon = gridlon(inregion);

gslat = groundStationsLocationslats;
gslon = groundStationsLocationslons;

groundStations = groundStation(sc, gslat, gslon);

acStationsTrue = access(camSensor, groundStations);
acStationsTrue.LineColor = 'green';
acStationsFalse = access(sat, groundStations);
acStationsFalse.LineColor = 'red';

% intvlsTrue = accessIntervals(acStationsTrue)
% intvlsFalse = accessIntervals(acStationsFalse)

%% Add Transmitter to LEO  Satellite
interferenceFreq = 2.99e9;                              % In Hz
rng("default");
txInterferingSat = transmitter(interferingSat, ...
    Frequency = interferenceFreq, ...                   % In Hz
    Power = 10+10*rand(1,numel(interferingSat)));       % In dBW
gaussianAntenna(txInterferingSat, ...
    DishDiameter = 0.2);                                % In m            % In m

play(sc)

