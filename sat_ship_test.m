clc; clear; close all;
%% Константы 
EarthRadius = physconst('EarthRadius');

%% Добавление LEO-спутника
% Создаём спутниковый сценарий
startTime = datetime(2024,3,24,0,15,45);
stopTime = startTime + hours(6);
sampleTime = 10; % секунд
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,ShowDetails=false);

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

% Добавляем сенсор к спутнику
SensorAngle = 179;
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
% Включить визуализацию поля видимости спутника
fov = fieldOfView(camSensor);
%play(sc)

%% Добавляем сетку базовых-станций (кораблей)
% target (-37.3397603503485, -14.5454154598156)

latlim = [-50 0]; 
lonlim = [-50 10];

spacingInLatLon = 2; % degrees

proj = projcrs(4087); %?? % 4087
spacingInXY = deg2km(spacingInLatLon)*1000; % meters
[xlimits,ylimits] = projfwd(proj,latlim,lonlim);
R = maprefpostings(xlimits,ylimits,spacingInXY,spacingInXY);
[X,Y] = worldGrid(R);

[gridlat,gridlon] = projinv(proj,X,Y);
AtlanticOcean = readgeotable("C:\Users\AKuznecova\Desktop\Github\Satelite-Model\AtlanticOcean\AtlanticOcean.shp");
SouthAtlanticOcean = AtlanticOcean(AtlanticOcean.name == "South Atlantic Ocean",:);
T = geotable2table(SouthAtlanticOcean,["latitude", "longitude"]);
[landlat,landlon] = polyjoin(T.latitude, T.longitude);
[landlat,landlon] = maptrimp(landlat,landlon,latlim,lonlim);

bufwidth = 0.5;
[landlatb,landlonb] = bufferm(landlat,landlon,bufwidth,"outPlusInterior");
australiab = geopolyshape(landlatb,landlonb);
gridpts = geopointshape(gridlat,gridlon);
inregion = isinterior(australiab,gridpts);

gslat = gridlat(inregion);
gslon = gridlon(inregion);
gs = groundStation(sc,gslat,gslon);

%% Антенна передатчика (спутник)
fq = 162e6; % Hz
txpower = 20; % dBW
antennaType = "Monopole"; % "Monopole";

effectiveLine = sqrt((500e3 + EarthRadius)^2 - EarthRadius^2); % касательная к земле по углу видимости 
halfBeamWidth = acosd(effectiveLine/EarthRadius);

if antennaType == "Monopole"
    lambda = physconst('lightspeed')/fq; % meters
    monopoleHeight = lambda/2;
    monopoleWidth = 1e-2;
    mpl = monopole('Height', monopoleHeight, 'Width', monopoleWidth); 
    tx = transmitter(sat, ...
        Frequency=fq, ...
        Power=txpower, ...
        Antenna=mpl); 
end


%% Антенна приёмника (корабли)
isotropic = arrayConfig(Size=[1 1]);
rx = receiver(gs,Antenna=isotropic);
pattern(tx,Size=500000);

%% 
[doppler_shifts, lat0, lon0] = satcoverage(gridpts,sc,startTime,inregion,halfBeamWidth);

figure
worldmap(latlim,lonlim)
mlabel south

colormap turbo
geoshow(gridlat,gridlon,(doppler_shifts), DisplayType="surface")
geoshow(SouthAtlanticOcean,FaceColor="none")

cBar = contourcbar;
title(cBar,"Hz");
title("Doppler Shift at " + string(startTime) + " UTC")


function [coveragedata, lat0, lon0] = satcoverage(gridpts,sc,timeIn,inregion,halfBeamWidth)

    % Get satellites and ground station receivers
    sats = sc.Satellites;
    rxs = [sc.GroundStations.Receivers];

    % Compute the latitude, longitude, and altitude of all satellites at the input time
    lla = states(sats,timeIn,"CoordinateFrame","geographic");

    lat0 = lla(1,1,1);
    lon0 = lla(2,1,1);

    % Initialize coverage data
    coveragedata = NaN(size(gridpts));

    for satind = 1
        % Create a geopolyshape for the satellite field-of-view
        fov = fieldOfViewShape(lla(:,1,satind), halfBeamWidth);

        % Find grid and rx locations which are within the field-of-view
        gridInFOV = isinterior(fov,gridpts);
        rxInFOV = gridInFOV(inregion);

        % Compute sigstrength at grid locations using temporary link objects
        gridsigstrength = NaN(size(gridpts));
        if any(rxInFOV)
            tx = sats(satind).Transmitters;
            lnks = link(tx,rxs(rxInFOV));
            % rxsigstrength = sigstrength(lnks,timeIn)+30; % Convert from dBW to dBm
            rxsigstrength = dopplershift(sats(satind),sc.GroundStations((rxInFOV)), timeIn, Frequency=162e6);
            gridsigstrength(inregion & gridInFOV) = rxsigstrength;
            delete(lnks)
        end

        % Update coverage data with maximum signal strength found
        coveragedata = max(gridsigstrength,coveragedata);
    end
end

function satFOV = fieldOfViewShape(lla,halfViewAngle)

    % Find the Earth central angle using the beam view angle
    if isreal(acosd(sind(halfViewAngle)*(lla(3)+earthRadius)/earthRadius))
        % Case when Earth FOV is bigger than antenna FOV 
        % The Earth central angle, whose value serves as a limiting condition to 
        % ﬁnd all the points on the Earth surface that fall into the instantaneousfootprint 
        % of the antenna:
    %     earthCentralAngle = 90-acosd(sind(halfViewAngle)*(lla(3)+earthRadius)/earthRadius)-halfViewAngle;
    % else
        % Case when antenna FOV is bigger than Earth FOV 
        %earthCentralAngle = 90-halfViewAngle;
        earthCentralAngle = 20; % формулу от высоты
    end

    % Create a buffer zone centered at the position of the satellite with a buffer of width equaling the Earth central angle
    [latBuff,lonBuff] = bufferm(lla(1),lla(2),earthCentralAngle,"outPlusInterior",100);

    % Handle the buffer zone crossing over -180/180 degrees
    if sum(abs(lonBuff) == 180) > 2
        crossVal = find(abs(lonBuff)==180) + 1;
        lonBuff(crossVal(2):end) = lonBuff(crossVal(2):end) - 360 *sign(lonBuff(crossVal(2)));
    elseif sum(abs(lonBuff) == 180) == 2
        lonBuff = [lonBuff; lonBuff(end); lonBuff(1); lonBuff(1)];
        if latBuff(1) > 0
            latBuff = [latBuff; 90; 90; latBuff(1)];
        else
            latBuff = [latBuff; -90; -90; latBuff(1)];
        end
    end

    % Create geopolyshape from the resulting latitude and longitude buffer zone values
    satFOV = geopolyshape(latBuff,lonBuff);
end



% %% Создаём связь между камерой и базовой станцией
% % БС в поле видимости спутника, но связь с камерой не установлена - красный
%  acStationsFalse = access(sat, [groundStations targetGs]);
% acStationsFalse.LineColor = 'red';
% 
% % Связь между камерой и БС установлена - зелёный
% acStationsTrue = access(camSensor, [groundStations targetGs]);
% acStationsTrue.LineColor = 'green';
% 
% %% Определяем отсчёты времени, когда связь м/у спутником и объектом устанволена
% intvlsTrue = accessIntervals(acStationsTrue)
% intvlsFalse = accessIntervals(acStationsFalse)
% 
% %% Add Transmitter to LEO  Satellite
% antennaType = "Monopole";
% if antennaType == "Monopole"
%     lambda = physconst('lightspeed')/fq; % meters
%     monopoleHeight = lambda/2;
%     monopoleWidth = 1e-2;
%     mpl = monopole('Height', monopoleHeight, 'Width', monopoleWidth); 
%     tx = transmitter(iridiumSatellites, ...
%         Frequency=fq, ...
%         Power=txpower, ...
%         Antenna=mpl); 
% end
% 
% % Play the scenario.
% play(sc);
% 
% %% Add Receiver to ground object
% GroundStationsAntenna = arrayConfig("Size",[1 1]); % Create an isotropic antenna element
% GroundStationsReceiver = receiver(targetGs, ...
%     Antenna=GroundStationsAntenna);
% pattern(GroundStationsReceiver,antennaf, Size=50000);
% %% Определяем расстояние до спутнка
% a = [ intvlsTrue.StartTime(3) intvlsTrue.EndTime(3)]
% for i = 1:length(a)
%     [targetGs_data(i).azimuth, targetGs_data(i).elevation, targetGs_data(i).range, targetGs_data(i).timeOut] = aer(sat,targetGs, a(i));
%     ship1(i).Distance = distance (targetGs.Latitude, targetGs.targetGs.Longitude,...
%         satelliteTransmitter);
% end
% 
% 
