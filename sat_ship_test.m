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
gs = groundStation(sc,Name="Location to Catch", ...
    Latitude=-16.943443,Longitude=177.835693);    

%% Поле видимости антенны к спутнику
% Common approach is to fix the sensor on a gimbal and orient the sensor 
% by manuevering the gimbal, rather than the spacecraft body itself
g = gimbal(sat);

% Направляем спутник на цель при помощи gimbal'a
%pointAt(g,gs);

% Добавляем сенсор к спутнику
camSensor = conicalSensor(g, "MaxViewAngle", 60);

%% Создаём связь между камерой и базовой станцией
% Связь между камерой и БС установлена - зелёный
%acTrue = access(camSensor, gs);
%acTrue.LineColor = 'green';
% БС в поле видимости спутника, но связь с камерой не установлена - красный
acFalse = access(sat, gs);
acFalse.LineColor = 'red';

%% Визуализация сценария 
% 
v = satelliteScenarioViewer(sc,"ShowDetails", true);
fov = fieldOfView(camSensor)
show(sat);

%% Определяем отсчёты времени, когда связь м/у спутником и объектом устанволена
intvlsTrue = accessIntervals(acTrue)
intvlsFalse = accessIntervals(acFalse)




%% Добавляем корабли
