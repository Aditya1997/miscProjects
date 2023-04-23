clear all, clc

% Articles
% https://thesportjournal.org/article/normative-profiles-for-serve-speed-for-the-training-of-the-serve-and-reception-in-volleyball/
% https://www.mdpi.com/2076-3417/9/19/4007
% https://www.mdpi.com/2073-4441/14/17/2593 - Drag/lift spinning sphere
% https://www.engineersedge.com/calculators/magnus_effect_calculator_15766.htm
% Note, Reynolds Re ranges from 2e5 (v = 10 m/s) to 6e5 (v = 28 m/s)

%% Arrays and Variables

x = zeros(1000,1);
y = zeros(1000,1);
xdrag = zeros(1000,1);
ydrag = zeros(1000,1);
xdragmagnus = zeros(1000,1);
ydragmagnus = zeros(1000,1);
Vx = zeros(1000,1);
Vy = zeros(1000,1);
Vxdrag = zeros(1000,1);
Vydrag = zeros(1000,1);
Vxdragmagnus = zeros(1000,1);
Vydragmagnus = zeros(1000,1);
xnet = zeros(300,1);
ynet = zeros(300,1);
xfloor = zeros(300,1);
yfloor = zeros(300,1);

yfinished = 0;
ydfinished = 0;
ydmfinished = 0;

% Serve variables
angledeg = 10; % angle in degrees
angle = angledeg*pi/180; % angle for calcs
Vi = 17; % velocity in meters/s
x(1) = -.5;
y(1) = 2.70;

% angledeg = 10; % angle in degrees
% angle = angledeg*pi/180; % angle for calcs
% Vi = 22; % velocity in meters/s
% x(1) = -.5;
% y(1) = 2.70;

% Ball starts
ydrag(1) = y(1);
ydragmagnus(1) = y(1);
xdrag(1) = x(1);
xdragmagnus(1) = x(1);
Vix = Vi*cos(angle);
Viy = Vi*sin(angle);
Vx(1) = Vix;
Vy(1) = Viy;
Vxdrag(1) = Vix;
Vydrag(1) = Viy;
Vxdragmagnus(1) = Vix;
Vydragmagnus(1) = Viy;

% Making net and floor
% Net
xnet(1) = 2.43;
ynet(1) = 0;

% Floor
xfloor(1) = 0;
yfloor(1) = 0;

for j = 2:300
    xnet(j) = 18.288/2;
    ynet(j) = ynet(j-1)+(2.4384/300);
    xfloor(j) = xfloor(j-1)+(18.288/300);
    yfloor(j) = 0;
end

%% Creating arrays based on accel/velocity/position
timesteps = 500;
for i = 2:timesteps
    timestepsize = 3/timesteps;
    t = i*timestepsize;
    % Drag variables, https://www.mdpi.com/2076-3417/9/19/4007
    radius = .105;
    A = pi*radius^2;
    rho = 1.293; % density of air
    %Cd = .3; % Ranges from ~.45 for v=10 to .2 for v >= 18
    Cd = .2; % Ranges from ~.45 for v=10 to .2 for v >= 18
    m = .28; % mass in kg
    % Drag force calculation
    if Vy(i-1) > 0
        Dy = -Cd*rho*((Vydrag(i-1))^2)*A/2;
    else
        Dy = Cd*rho*((Vydrag(i-1))^2)*A/2;
        if Dy > 9.81*m
            Dy = 9.81*m;
        end
    end
    Dx = -Cd*rho*((Vxdrag(i-1))^2)*A/2;

    % Calculating velocity and position with drag
    Vxdrag(i) = Vxdrag(i-1) + (Dx/m)*timestepsize;
    Vydrag(i) = Vydrag(i-1) + (-9.81+Dy/m)*timestepsize;
    xdrag(i) = xdrag(i-1) + Vxdrag(i)*timestepsize;
    ydrag(i) = ydrag(i-1) + Vydrag(i)*timestepsize;

    % Calculating effects of Magnus force
    velapprox = Vi;
    omega = 30; % radians per s
    S = radius*omega/velapprox;
    disp(S);
    %CL = 0; % Float serve, no spin
    %CL = .1; % Lift coefficient for Magnus, approx based on S
    CL = .2; % Lift coefficient for Magnus, approx based on S
    Vdragmagnustotal = sqrt((Vydragmagnus(i-1)^2)+(Vxdragmagnus(i-1)^2)); % magnitude of ball velocity
    Vdragmagnusangle = atan((Vydragmagnus(i-1))/(Vxdragmagnus(i-1))); % angle of ball vel traj
    disp(rad2deg(Vdragmagnusangle));
    % FMagnus = .5*Cl*rho*A*v^2
    FMagnus = .5*CL*rho*A*(Vdragmagnustotal)^2;
    if Vdragmagnusangle > 0
        FMagnusx = FMagnus*sin(Vdragmagnusangle);
        FMagnusy = FMagnus*-cos(Vdragmagnusangle);
    else
        FMagnusx = FMagnus*-sin(Vdragmagnusangle);
        FMagnusy = FMagnus*-cos(Vdragmagnusangle);
    end
    disp(FMagnusx);
    disp(FMagnusy);
    DMx = Dx + FMagnusx; % Drag force + Magnus force
    DMy = Dy + FMagnusy; % Drag force + Magnus force
    % Calculating velocity and position with drag + Magnus force
    Vxdragmagnus(i) = Vxdragmagnus(i-1) + (DMx/m)*timestepsize;
    Vydragmagnus(i) = Vydragmagnus(i-1) + (-9.81+DMy/m)*timestepsize;
    xdragmagnus(i) = xdragmagnus(i-1) + Vxdragmagnus(i)*timestepsize;
    ydragmagnus(i) = ydragmagnus(i-1) + Vydragmagnus(i)*timestepsize;

    % No drag
    Dx = 0;
    Dy = 0;
    % m = 1;
    Vx(i) = Vx(1) + (Dx/m)*t;
    Vy(i) = Vy(1) + (-9.81+Dy/m)*t;
    x(i) = x(1) + Vix*t + .5*(Dx/m)*t^2;
    y(i) = y(1) + Viy*t + .5*(-9.81+Dy/m)*(t^2);
end

%% Ending arrays upon floor
for yvals=1:1000
    if ydragmagnus(yvals) <= 0 && ydmfinished == 0
        endvaldm = yvals;
        ydmfinished = 1;
    end
    if ydrag(yvals) <= 0 && ydfinished == 0
        endvald = yvals;
        ydfinished = 1;
    end
    if y(yvals) <= 0 && yfinished == 0
        endval = yvals;
        yfinished = 1;
    end
end
xdragmagnusbeforefloor = xdragmagnus(1:endvaldm);
ydragmagnusbeforefloor = ydragmagnus(1:endvaldm);
xdragbeforefloor = xdrag(1:endvald);
ydragbeforefloor = ydrag(1:endvald);
xbeforefloor = x(1:endval);
ybeforefloor = y(1:endval);

%% Plotting arrays
plot(xnet,ynet);
hold on
plot(xfloor,yfloor);
hold on
plot(xbeforefloor,ybeforefloor);
hold on
plot(xdragbeforefloor,ydragbeforefloor);
hold on
plot(xdragmagnusbeforefloor,ydragmagnusbeforefloor);
axis equal
legend('','','no effects','drag only','magnus+drag')
hold off;
