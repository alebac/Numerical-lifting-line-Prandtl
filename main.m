clear all
close all
clc

%% Input data
stationsWing = 1000;
fourierModes = 150;

rootChord = 0.2;
taperRatio = 1;
wingspan = 0.8;
airspeed = [30 0 0];
density = 1.225;                % kg/m^3

alphaFly = 4; 
twist = 0;

alpha0lift = 0;
alphaRad = alphaFly * pi/180;
alpha0liftRad = alpha0lift * pi/180;
clalpha = 5.73; %2*pi;

%% Propellers

hubPositions = 0.4;
diameters = 0.15;
thrust = 5;
rpm = -25400;                                   % rounds per minute
propellerRotationSpeed = rpm/60 * 2*pi;         % radians per second

spinnerRadius = 0.1 * max(diameters);

%% PRANDTL %%

% Wing twist
clear alpha;
alphaTwist = abs((0:stationsWing) / stationsWing - 1/2)*2 * twist;
alpha = alphaFly + alphaTwist;

% Wing discretization
stations = 0:stationsWing;
dy = wingspan/stationsWing;
y = stations'*dy - wingspan/2;
y(1) = y(1) + 0.0001 * wingspan;        % this is to make sure that a chord = 0 
y(end) = y(end) - 0.0001 * wingspan;    % at the tip doesn't create a mess
theta = acos(-2*y/wingspan);

%% Chord
% chord = interp1(a(2,:), a(1,:), abs(y3));
% chord = rootChord * sqrt(1-(2*y3/wingspan).^2); 
chord = ones(stationsWing+1,1) * rootChord;
% chord = rootChord - 2*abs(y3)/wingspan * (1-taperRatio)*rootChord;

%% Surface
surface = trapz(y,chord);

%% Propellers
propData = propellersPosition(hubPositions, diameters, stationsWing, wingspan, ...
    propellerRotationSpeed, airspeed, density, thrust);
inducedVelocity = induced_velocity_propeller(propData, airspeed, ...
    spinnerRadius);

%% Solving the system
M = zeros(stationsWing+1,fourierModes);

for m = 1:stationsWing+1
    for n = 1:fourierModes
        M(m,n) = sin(n*theta(m)) * ( sin(theta(m)) ...
            + n * (clalpha*chord(m)/(4*wingspan)) );
    end
end

% Propellers induction
Vtot = zeros(stationsWing+1,3);
Vtot(:,1) = ones(stationsWing+1,1)*norm(airspeed);
Vtot = Vtot + inducedVelocity;

normVtot = zeros(stationsWing+1,1);
for i = 1:stationsWing+1
    normVtot(i) = norm(Vtot(i,:));
end

alphaP = angle_of_attack( Vtot );

b = zeros(stationsWing+1,1);
for m = 1:stationsWing+1
        b(m) = (clalpha*chord(m)/(4*wingspan)) * ...
            normVtot(m)/norm(airspeed) * ((alpha(m)-alpha0lift) * pi/180 ...
            - Vtot(m,3)/normVtot(m)) * sin(theta(m));
end

% b = ones(k+1,1) * (alpha-alpha0lift) * pi/180;
A = M\b;

G = zeros(stationsWing+1,fourierModes);
for m = 1:stationsWing+1 
    for n = 1:fourierModes
        G(m,n) = 2*wingspan*norm(airspeed) * sin( n* theta(m) );
    end
end
gamma = G*A;

Z = zeros(stationsWing+1,fourierModes);
for m = 1:stationsWing+1
    for n = 1:fourierModes
        Z(m,n) = n * sin(n* theta(m)) / sin(theta(m));
    end
end
alphai = Z*A;

%% Post processing
AR = wingspan^2 / surface;
CL = A(1) * pi * AR;
CDi = pi * AR * A(1)^2 * ( 1 + sum( (2:fourierModes)' .* (A(2:end)/A(1)).^2 ) );

AR_check = CL^2/(pi*CDi);
e1 = CL^2/(pi*CDi*AR);
e2 = 1 / ( 1 + sum( (2:fourierModes)' .* (A(2:end)/A(1)).^2 ) );

%% Propellers induction
downwash = norm(airspeed) * sin(alphai);
Vtot(:,3) = downwash;

alphai = angle_of_attack( Vtot );

clprandtl = gamma ./ (1/2 * norm(airspeed) .* chord);

%% Plots
plots