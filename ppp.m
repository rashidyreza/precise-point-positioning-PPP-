% Reza Rashidi
% email: rashidyreza76@gmail.com
% PPP
clc;clear;close all;
% ==============================================================
% [rinexobs] = read_rinex_obs('bcov0010.21o');
load('matlab.mat')
% ==============================================================
algorithm_type = 1;
Filetype = 'sp3';% Filetype = 'sp3';
readnewfile = 'no';
usethisfile = rinexobs;
r0rcv = 'r0rinex';
[rinexobs] = T_emission_collins_klobuchar_XYZsat(algorithm_type,Filetype,readnewfile,usethisfile,r0rcv);
% ========================3=dtsat sp3=========================
[SP3_XYZ] = ReadSP3('igs21385.sp3');
sp3time = unique(SP3_XYZ.data(:,SP3_XYZ.col.gps_seconds));
TOW = rinexobs.data(:,rinexobs.col.TOW); 
Find = [];
for j = 1:length(sp3time)
    i = sp3time(j);
    Find = [Find ;find(TOW==i)];
end
rinexobs.data = rinexobs.data(Find,:);
% =======================SUNBern============================
PRN_Time=SUNBern(SP3_XYZ);
TOW = rinexobs.data(:,rinexobs.col.TOW); 
PRN = rinexobs.data(:,rinexobs.col.PRN); 
PTn = [[1:size(TOW)]',TOW,PRN];
Find = [];
for i = 1: size(PRN_Time,1)
if find(TOW == PRN_Time(i,2) )
    ff = find(TOW == PRN_Time(i,2) );
    fff = PTn(ff,:);
    if find(fff(:,3) == PRN_Time(i,1) )
        FFF = find(fff(:,3) == PRN_Time(i,1) );
        Find = [Find;fff(FFF,1)];
    end
end
end
rinexobs.data(Find,:) = [];
%===============================================================
[ant_PCV] = read_antenna_PCV('igs14.atx', {'SEPCHOKE_B3E6   SPKE'});
options.system.gps = 1;options.system.glo =0;options.system.gal =0;options.system.bds =0;
[atx] = r_antx('igs14.atx',  {'SEPCHOKE_B3E6   SPKE'} ,options);
%===============================================================
Find = find( rinexobs.data(:,rinexobs.col.L1)== 0 );
rinexobs.data(Find,:) =  [];
Find = find( rinexobs.data(:,rinexobs.col.L2)== 0 );
rinexobs.data(Find,:) =  [];
Find = find( rinexobs.data(:,rinexobs.col.C1)== 0 );
rinexobs.data(Find,:) =  [];
Find = find( rinexobs.data(:,rinexobs.col.P2)== 0 );
rinexobs.data(Find,:) =  [];
%====================cut of angle================================
Find = find( rinexobs.data(:,rinexobs.col.Elevation)< 5 );
rinexobs.data(Find,:) =  [];
%%======================epoch====================================
TOW = rinexobs.data(:,rinexobs.col.TOW); 
epoch = unique(TOW);
epoch = epoch(1+1);
m = find(TOW ==epoch ); m = m(1)-1;
rinexobs.data = rinexobs.data(1:m,:);
%=================================================
[windup,rsun]=windup_ppp(rinexobs);
%================================================
for i = 1:m
PCV_receiver= PCV_interp(ant_PCV, 90-rinexobs.data(i,rinexobs.col.Elevation),wrapTo360( rinexobs.data(i,rinexobs.col.Azimuth)), 1, 1);
PCO_satellite = sat_apc( rinexobs.data(i,18:20),rinexobs.r0',rsun(i,:), atx.sat.neu , rinexobs.data(i,rinexobs.col.PRN), 1);
rinexobs.data(i,18:20) = rinexobs.data(i,18:20) + PCO_satellite';
end


%=================================================
L1 = rinexobs.data(:,rinexobs.col.L1);
L2 = rinexobs.data(:,rinexobs.col.L2);
freqs1=1575.42e6;   
freqs2=1227.60e6;
Phic_obs = (L1*freqs1^2 - L2*freqs2^2)/(freqs1^2 - freqs2^2);
C1 = rinexobs.data(:,rinexobs.col.C1);
P2 = rinexobs.data(:,rinexobs.col.P2);
Rc_obs = (C1*freqs1^2 - P2*freqs2^2)/(freqs1^2 - freqs2^2);
c = 299792458;
Xsat = rinexobs.data(:,rinexobs.col.Xsat);
Ysat = rinexobs.data(:,rinexobs.col.Ysat);
Zsat = rinexobs.data(:,rinexobs.col.Zsat);
Tr0  = rinexobs.data(:,rinexobs.col.Tr0);
Mwet  = rinexobs.data(:,rinexobs.col.Mwet);
dtsat = rinexobs.data(:,rinexobs.col.dtsat);
TOW = rinexobs.data(:,rinexobs.col.TOW); 
PRN = rinexobs.data(:,rinexobs.col.PRN);
TOW = rinexobs.data(:,rinexobs.col.TOW); 
delta_ro_rel = rinexobs.data(:,rinexobs.col.delta_ro_rel);
yphi = Phic_obs+c*dtsat - Tr0-delta_ro_rel-windup;
yR = Rc_obs+c*dtsat - Tr0-delta_ro_rel;
y = [yphi;yR];
%-------------------------------------------------------------------------------
rcv0 = rinexobs.r0;
tr0 = 1;
cdtrcv = 1;
Bi0 = ones(m,1);
X0 = [rcv0;cdtrcv;tr0;Bi0];
%===============================================================================
Q = inv([0.003^2*eye(m) zeros(m);zeros(m) 0.3^2*eye(m)]);
dxhat = 1;
while norm(dxhat) > 10^-5
dX = -Xsat+rcv0(1);
dY = -Ysat+rcv0(2);
dZ = -Zsat+rcv0(3);
ro0 = sqrt( dX.^2 +dY.^2 + dZ.^2 );
y0 = ro0+cdtrcv + Mwet*tr0;
dy = y-[y0+eye(m)*Bi0;y0];
A = [dX./ro0 , dY./ro0 , dZ./ro0 , ones(m,1),Mwet,eye(m);...
    dX./ro0 , dY./ro0 , dZ./ro0 , ones(m,1),Mwet,zeros(m,m)];
dxhat = inv(A'*Q*A)*A'*Q*dy;
X0 = X0 + dxhat;
rcv0 = X0(1:3);cdtrcv = X0(4);tr0 = X0(5);Bi0 = X0(6:end)
end
Qx = inv(A'*A);
rinexobs.r0 - rcv0