% Reza Rashidi
% email: rashidyreza76@gmail.com

%{
inputs
                    rinex observation File
                    NAVIGATION FILE or sp3 FILE
                    Filetype: 'sp3' or 'navigation'
                    algorithm_type for calculate T_emission
                           1: for A Pseudorange-Based Algorithm
                           2: for A Purely Geometric Algorithm
                    r0rcv for second Algorithm
                          r0rcv = [X;Y;Z]
                           or
                 r0rcv = 'r0rinex' to use coordinates in rinexobs header;
                 readnewfile = 'yes' or 'no'
                 if no I will use old file of rinex observation File
%}

%{
outinputs
                 Time of emission
                 dt
                 [X,Y,Z] sat in Time of emission
                 dtsat

%}
%-------------------------read files------------------------------
function [rinexobs,satnav] = T_emission_collins_klobuchar_XYZsat(algorithm_type,Filetype,readnewfile,usethisfile,r0rcv)
readnewfile = sum(double(readnewfile));
switch Filetype
    case 'sp3'
        sp3fname = uigetfile('*.sp3*');     % sp3 FILE
        [SP3_XYZ] = ReadSP3(sp3fname);
       alfa=[1 1 1 1];
       beta = alfa;
    case 'navigation'
        Filename = uigetfile('*.21n*');     % NAVIGATION FILE
        [nav,Date,alfa,beta]=ReadNavigation(Filename,38);
end
if readnewfile == sum(double('yes'))
    filename = uigetfile('*.21o*');             % OBSERVATION FILE
    [rinexobs] = read_rinex_obs(filename);
    rec_xyz = rinexobs.r0;
else
    rinexobs = usethisfile;
    rec_xyz = rinexobs.r0;
end
%============================================================
if r0rcv== 'r0rinex'
    r0rcv_xyz = rec_xyz;
else
    r0rcv_xyz = r0rcv;
end
% GNSS Type G for GPS
Find = find( rinexobs.data(:,rinexobs.col.GNSStype)== double('G') );
rinexobs.data =  rinexobs.data(Find,:);
global GM we c
GM = 3.986004418e14;
we = 7.2921151467e-5;
c = 2.99792458e8;
dt=0;
i=4;
while find(dt==0)>0
    dt = rinexobs.data(:,i)/c;
    Find = find( dt==0 );
    i=i+1;
    if i==10
        rinexobs.data(Find,:)=[];
        dt(Find,:)=[];
        break
    end
end
switch algorithm_type
    % A Pseudorange-Based Algorithm
    case 1
        switch Filetype
            case 'sp3'
                PRN = unique(SP3_XYZ.data(:,SP3_XYZ.col.PRN))';
                for i = PRN
                    dtsat = 0;
                    while 1
                        Find = find( SP3_XYZ.data(:,SP3_XYZ.col.PRN) == i );
                        pointY = SP3_XYZ.data(Find,SP3_XYZ.col.dtsat)*1e-6;
                        pointX = SP3_XYZ.data(Find,SP3_XYZ.col.gps_seconds);
                        Find = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                        Temission = rinexobs.data(Find,rinexobs.col.TOW) - dt(Find) - dtsat;
                        dtsat0 = lagrange_interpolation(Temission,pointX,pointY,10);
                        Norm = norm ( dtsat - dtsat0);
                        dtsat = dtsat0;
                        if Norm<1e-10
                            rinexobs.data(Find,rinexobs.col.dtsat) = dtsat;
                            rinexobs.data(Find,rinexobs.col.dt_sat_rcv) = dt(Find);
                            % Temission
                            Temission = rinexobs.data(Find,rinexobs.col.TOW) - dt(Find) - dtsat;
                            rinexobs.data(Find,rinexobs.col.Temission) = Temission;
                            break;
                        end
                    end
                    Findsp3 = find( SP3_XYZ.data(:,SP3_XYZ.col.PRN) == i );
                    Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                    pointX = SP3_XYZ.data(Findsp3,SP3_XYZ.col.gps_seconds);
                    pointYx = SP3_XYZ.data(Findsp3,SP3_XYZ.col.X);
                    pointYy = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Y);
                    pointYz = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Z);
                    XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                    YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                    ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                    rsat = [XX , YY , ZZ ];
                    rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat) = rsat;
                end
            case 'navigation'
                PRN = unique(nav.data(:,nav.col.PRN))';
                Prn = nav.data(:,nav.col.PRN);
                for i = PRN
                    Find = find(Prn == i); Find = Find(1);
                    TGD = nav.data(Find,nav.col.TGD)*299792458;
                    rinexobs.data( find(rinexobs.data(:,rinexobs.col.PRN)==i),rinexobs.col.TGD) = TGD;
                end
                PRN = unique(nav.data(:,nav.col.PRN))';
                for i = PRN
                    dtsat = 0;
                    while 1
                        Findnav = find( nav.data(:,nav.col.PRN) == i );
                        Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                        Temission = rinexobs.data(Findobs,rinexobs.col.TOW) - dt(Findobs) - dtsat;
                        TOE = nav.data(Findnav,nav.col.TOE);
                        Nearest = [];
                        for j=1:length(Temission)
                            MIN = abs(TOE-Temission(j));
                            nnn = find(MIN==min(MIN));
                            Nearest(j,1) = Findnav(nnn);
                        end
                        a0 = nav.data(Nearest,nav.col.SV_Clock_Bias);
                        a1 = nav.data(Nearest,nav.col.SV_Clock_Drift);
                        a2 = nav.data(Nearest,nav.col.SV_Clock_Drift_Rate);
                        toe = nav.data(Nearest,nav.col.TOE);
                        [AAAA,drel] =  NAV2ECEF(nav,rinexobs,i);
                        dtsat0 = a0 + a1.*(Temission - toe) + a2.*(Temission - toe).^2 + drel;
                        Norm = norm ( dtsat - dtsat0);
                        dtsat = dtsat0;
                        if Norm<1e-10
                            rinexobs.data(Findobs,rinexobs.col.dtsat) = dtsat;
                            rinexobs.data(Findobs,rinexobs.col.dt_sat_rcv) = dt(Findobs);
                            % Temission
                            Temission = rinexobs.data(Findobs,rinexobs.col.TOW) - dt(Findobs) - dtsat;
                            rinexobs.data(Findobs,rinexobs.col.Temission) = Temission;
                            break;
                        end
                    end
                    [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                end
        end
    case 2
        switch Filetype
            case 'sp3'
                rinexobs.data(:,rinexobs.col.Temission) = rinexobs.data(:,rinexobs.col.TOW);
                PRN = unique(SP3_XYZ.data(:,SP3_XYZ.col.PRN))';
                for i = PRN
                    Findsp3 = find( SP3_XYZ.data(:,SP3_XYZ.col.PRN) == i );
                    Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                    pointX = SP3_XYZ.data(Findsp3,SP3_XYZ.col.gps_seconds);
                    pointYx = SP3_XYZ.data(Findsp3,SP3_XYZ.col.X);
                    pointYy = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Y);
                    pointYz = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Z);
                    pointY = SP3_XYZ.data(Findsp3,SP3_XYZ.col.dtsat)*1e-6;
                    while 1
                        Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                        XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                        YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                        ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                        rsat = [XX , YY , ZZ ];
                        rsat_r0rcv = (rsat - ones(length(Findobs),1)*r0rcv_xyz');
                        dt = sqrt(sum(rsat_r0rcv.^2,2))/c;
                        rinexobs.data(Findobs,rinexobs.col.Temission) = rinexobs.data(Findobs,rinexobs.col.TOW) - dt;
                        Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                        XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                        YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                        ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                        rsat0 = [XX , YY , ZZ ];
                        Norm = norm ( rsat0 - rsat);
                        if Norm<1e-3
                            rinexobs.data(Findobs,rinexobs.col.dt_sat_rcv) = dt;
                            rinexobs.data(Findobs,rinexobs.col.Temission) = (rinexobs.data(Findobs,rinexobs.col.TOW) - dt -...
                                rinexobs.data(Findobs,rinexobs.col.dtrcv));
                            Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                            XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                            YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                            ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                            rsat = [XX , YY , ZZ ];
                            rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat) = rsat;
                            break;
                        end
                    end
                    dtsat = lagrange_interpolation(Temission,pointX,pointY,10);
                    rinexobs.data(Findobs,rinexobs.col.dtsat) = dtsat;
                end
            case 'navigation'
                
                PRN = unique(nav.data(:,nav.col.PRN))';
                Prn = nav.data(:,nav.col.PRN);
                for i = PRN
                    Find = find(Prn == i); Find = Find(1);
                    TGD = nav.data(Find,nav.col.TGD)*299792458;
                    rinexobs.data( find(rinexobs.data(:,rinexobs.col.PRN)==i),rinexobs.col.TGD) = TGD;
                end
                
                rinexobs.data(:,rinexobs.col.Temission) = rinexobs.data(:,rinexobs.col.TOW);
                PRN = unique(nav.data(:,nav.col.PRN))';
                for i = PRN
                    Findnav = find( nav.data(:,nav.col.PRN) == i );
                    Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                    while 1
                        [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                        rsat = rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat);
                        rsat_r0rcv = (rsat - ones(length(Findobs),1)*r0rcv_xyz');
                        dt = sqrt(sum(rsat_r0rcv.^2,2))/c;
                        rinexobs.data(Findobs,rinexobs.col.Temission) = rinexobs.data(Findobs,rinexobs.col.TOW) - dt;
                        [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                        rsat0 = rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat);
                        Norm = norm ( rsat0 - rsat);
                        if Norm<1e-3
                            rinexobs.data(Findobs,rinexobs.col.dt_sat_rcv) = dt;
                            rinexobs.data(Findobs,rinexobs.col.Temission) = (rinexobs.data(Findobs,rinexobs.col.TOW) - dt -...
                                rinexobs.data(Findobs,rinexobs.col.dtrcv));
                            [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                            break;
                        end
                    end
                    TOE = nav.data(Findnav,nav.col.TOE);
                    Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                    Nearest = [];
                    for j=1:length(Temission)
                        MIN = abs(TOE-Temission(j));
                        nnn = find(MIN==min(MIN));
                        Nearest(j,1) = Findnav(nnn);
                    end
                    a0 = nav.data(Nearest,nav.col.SV_Clock_Bias);
                    a1 = nav.data(Nearest,nav.col.SV_Clock_Drift);
                    a2 = nav.data(Nearest,nav.col.SV_Clock_Drift_Rate);
                    toe = nav.data(Nearest,nav.col.TOE);
                    [AAAA,drel] =  NAV2ECEF(nav,rinexobs,i);
                    dtsat = a0 + a1.*(Temission - toe) + a2.*(Temission - toe).^2 + drel;
                    rinexobs.data(Findobs,rinexobs.col.dtsat) = dtsat;
                end
        end
end
dt = rinexobs.data(:,rinexobs.col.dt_sat_rcv);
dw = we*dt;
for i = 1:length(dw)
    rsat = rinexobs.data(i,rinexobs.col.Xsat:rinexobs.col.Zsat)';
    rsat = Rotation(3,dw(i))*rsat;
    rinexobs.data(i,rinexobs.col.Xsat:rinexobs.col.Zsat) = rsat';
end
% compute Elevation & Azimuth
Res = cart_Geo(r0rcv_xyz,2);
rinexobs.Geodetic = Res;
phi = Res(1);
lambda = Res(2);
h = Res(3);
r0rcv_xyz = r0rcv_xyz(:);
K = size(rinexobs.data,1);
r_ct =  rinexobs.data(:,rinexobs.col.Xsat:rinexobs.col.Zsat)';
dr = r_ct - r0rcv_xyz*ones(1,K);
clear CC;
for k=1:K;
    CC(k,1)=norm(dr(:,k));
end
ro = dr./((CC*ones(1,3))');
eEe = [-sin(lambda),cos(lambda),0]';
nNn=[-cos(lambda)*sin(phi) ,-sin(lambda)*sin(phi) , cos(phi)]';
uUu= [cos(lambda)*cos(phi), sin(lambda)*cos(phi), sin(phi)]';
dU=[eEe nNn uUu]*dr;
Elevationa = asin( sum(ro.*(uUu*ones(1,k))) )*180/pi;
Az = atan2( sum(ro.*(eEe*ones(1,k))),sum(ro.*(nNn*ones(1,k))) )*180/pi;
for ii = 1:K
    d_ion(ii,1) = 0;
    Dataaa(ii,:) = GPSt2Date(rinexobs.data(ii,1),rinexobs.data(ii,2));
    JD(ii,1) = cal2gpstime(Dataaa(ii,:));
end
rinexobs.data(:,rinexobs.col.JD) = JD;
rinexobs.data(:,rinexobs.col.I1_klobuchar) = d_ion;
rinexobs.data(:,rinexobs.col.Azimuth:rinexobs.col.Elevation) = [Az',Elevationa'];
Data = GPSt2Date(rinexobs.data(100,1),rinexobs.data(100,2));
Data = Data(1:3);
DOYy = DOY(Data);
for k = 1:K
    [ Tr0,Mwet ] = collins( phi,h,DOYy,Elevationa(k)*pi/180 );
    rinexobs.data(k,rinexobs.col.Tr0) = Tr0;
    rinexobs.data(k,rinexobs.col.Mwet) = Mwet;
end
% delta ro_rel
rrcv0 = r0rcv_xyz*ones(1,K);
rsat0 = r_ct;
for k=1:K;
    rsat_rcv(k,1)=norm(dr(:,k));
    rrcv(k,1)=norm(rrcv0(:,k));
    rsat(k,1)=norm(rsat0(:,k));
end
delta_ro_rel = 2*GM*log( (rsat + rrcv + rsat_rcv)./(rsat + rrcv - rsat_rcv) )./(c^2);
delta_ro_rel = delta_ro_rel(:);
rinexobs.data(:,rinexobs.col.delta_ro_rel) = delta_ro_rel;

switch Filetype
    case 'sp3'
        satnav = SP3_XYZ;
    case 'navigation'
        satnav = nav;
end
[rinexobs]=Hatch_filter_GNSS(rinexobs);

end
function [R]=Rotation(i,j)
%Rotation matrix
if i==3
    R = [cos(j) sin(j) 0 ; -sin(j) cos(j) 0 ; 0 0 1];
else if i==2
        R = [cos(j) 0 -sin(j) ; 0 1 0 ; sin(j) 0 cos(j)];
    else if i==1
            R = [1 0 0 ; 0 cos(j) sin(j) ; 0 -sin(j) cos(j)];
        end
    end
end
end
function result = cart_Geo(XXX,type)
a = 6378137;
b = 6356752.3142;
e2=(a^2-b^2)/(a^2);
switch type
    case 1
        %%
        %phi lambda h to kartezian
        phi = XXX(1);
        Landa = XXX(2);
        h = XXX(3);
        N=a^2/sqrt(a^2*cos(phi)^2+b^2*sin(phi)^2);
        X=(N+h)*cos(phi)*cos(Landa);
        Y=(N+h)*cos(phi)*sin(Landa);
        Z=(N*(b^2/a^2)+h)*sin(phi);
        result = [X Y Z];
    case 2
        %--------------------------------------
        %kartezian to phi lambda h
        X = XXX(1);
        Y = XXX(2);
        Z = XXX(3);
        p=sqrt(X^2+Y^2);
        phi0=atan((Z/p)/(1-e2));
        N0=a/(1-e2*sin(phi0)^2)^.5;
        h0=p/cos(phi0)-N0;
        Landa=wrapTo2Pi(atan2(Y,X));
        repitition=0;delta_phi=1;delta_h=1;
        while abs(delta_phi)>10^-12 && abs(delta_h)>10^-5
            N=a/(1-e2*sin(phi0)^2)^.5;
            h=(p/cos(phi0))-N;
            phi=atan((Z/p)*((1-((e2*N)/(N+h)))^-1));
            delta_phi=phi-phi0;delta_h=h-h0;phi0=phi;h0=h;
            %             repitition=repitition+1
        end
        result = [phi Landa h];
end
end
function [ Tr0,Mwet ] = collins( phi0,H,D,elev )
phi = phi0*180/pi;
elev = elev*180/pi;
Trdz=2.3*exp(-0.11*H*10^(-3));
Trz0w=0.1;
[ad bd cd aw bw cw]=interpp(phi0,D);
at=2.53e-5;bt=5.49e-3;ct=1.14e-3;
m=(1+(ad/1+(bd/1+cd)))/(sind(elev)+(ad/(sind(elev)+(bd/(sind(elev)+cd)))));
dm=((1/sind(elev))-(1+(at/1+(bt/1+ct)))/(sind(elev)+(at/(sind(elev)+(bt/(sind(elev)+ct))))))*H;
Mdry=m+dm; 
Mwet=(1+(aw/1+(bw/1+cw)))/(sind(elev)+(aw/(sind(elev)+(bw/(sind(elev)+cw)))));
Tr0=Trdz*Mdry+Trz0w*Mwet;
end
function out=DOY(r2)
y=r2(1);m=r2(2);d=r2(3);
k=[31 28 31 30 31 30 31 31 30 31 30 31];
y=2000+y;
R1=abs(2012-y);
R2=(R1/4)-floor(R1/4);
if R2==0
    k(2)=29;
end
if m>1
    doy=sum(k(1:m-1))+d;
elseif m==1
    doy=d;
end
out=doy;
end
function Data = GPSt2Date(gps_week,gps_sec)
JD = 7*( gps_week + (gps_sec/(86400*7)) ) + 2444244.5;
[ Data ] = JD2Date(JD);
end
function [ Data ] = JD2Date(JD)
a = floor(JD + 0.5);
b = a + 1537;
c = floor( (b - 122.1)/365.25 );
d = floor(365.25*c);
e = floor( (b-d)/30.6001 );
D = b - d - floor(30.6001*e) +  (JD + 0.5) - floor(JD + 0.5);
M = e - 1 - 12*floor(e/14);
Y = c - 4715 - floor( (7 + M)/10 );
h = degrees2dms((D - floor(D))*24 );
D = floor(D);
H = h(1);
Min = h(2);
S = ceil(h(3));
S = h(3);
Data = [Y M D H Min S];
end
function d_ion = KlobucharAlgorithm(Az,El,phi0,lambda0,gps_sec,alfa,beta)
global c Re hgps phiP lambdaP
c=299792458;
Re=6378000;
hgps=350000;
phiP=deg2rad(78.3);
lambdaP=deg2rad(291);
AI=0;
PI=0;
%step1 % Earth-centred angle
psi=( (pi/2)- El - asin(Re*cos(El)/(Re+hgps)) );
%step2 % latitude of the IPP
phiIP=asin( sin(phi0)*cos(psi) + cos(phi0)*sin(psi)*cos(Az) );
%step3 % longitude of the IPP
lambdaIP=lambda0 + psi*sin(Az)/cos(phiIP);
% step4: Find the geomagnetic latitude of the IPP
phim = asin( sin(phiIP)*sin(phiP) + cos(phiIP)*cos(phiP)*cos(lambdaIP-lambdaP) );
% step5: Find the local time at the IPP
t1=43200*(lambdaIP/pi)+gps_sec;
t1 = mod(t1,86400);
% step6: Compute the amplitude of ionospheric delay
for i=1:4
    AI=AI+alfa(i,1)*(phim/pi)^(i-1);
end
if AI<0
    AI=0;
end
% step7: Compute the period of ionospheric delay
for i=1:4
    PI=PI+beta(i,1)*(phim/pi)^(i-1);
end
if PI<72000
    PI=72000;
end
%step 8: Compute the phase of ionospheric delay
XI=2*pi*(t1-50400)/PI;
% step9: Compute the slant factor (ionospheric mapping function)
F=(1-(Re*cos(El)/(Re+hgps))^2)^(-0.5);
% step10: Compute the ionospheric time delay
if abs(XI)<pi/2
    I1 =( 5e-9 + AI*cos(XI) )*F;
elseif abs(XI)>=pi/2
    I1 = 5e-9*F;
end
% Ik = ((f1/fk)^2)*I1;
d_ion=I1*c;
freqs1=1575.42e6;   
freqs2=1227.60e6;
alpha1=1/((freqs1/freqs2)^2-1);
d_ion = d_ion*alpha1;
end
function [rinexx]=Hatch_filter_GNSS(rinexx)
%--------------------------------------
N = 360;
freqs1=1575.42e6;   
freqs2=1227.60e6;
alpha1=1/((freqs1/freqs2)^2-1);
 TOW=rinexx.data(:,rinexx(1).col.TOW);
 interval=unique(TOW);
 interval=interval(2)-interval(1);
PRN=rinexx.data(:,rinexx.col.PRN);
rinex = rinexx;
k=1;
for i=[unique(PRN)]'
    rinex(k).PRN=rinexx(1).data(find(PRN==i),:);
    k=k+1;
end
%------------------------------------------------
for i=1:k-1
Arc=rinex(i).PRN(:,rinex(1).col.TOW)./interval;
Arc=[Arc;0]-[0;Arc];
Arc(1)=[];
Arc=find(abs(Arc)>1);
Arc=[0;Arc];
Arc=[Arc+1,Arc];
Arc=[Arc(1:end-1,1),Arc(2:end,2)];
rinex(i).Arc=Arc;
end
for Sat = [unique(PRN)]'
sat = find(unique(PRN)==Sat);
Arc = rinex(sat).Arc;
C1=rinex(sat).PRN(:,rinex(1).col.C1);
P2=rinex(sat).PRN(:,rinex(1).col.P2);
L1=rinex(sat).PRN(:,rinex(1).col.L1);
L2=rinex(sat).PRN(:,rinex(1).col.L2);
Time=rinex(sat).PRN(:,rinex(1).col.TOW);
L1_D=L1+2*alpha1*(L1-L2);
L1_LC=((freqs1^2)*L1-(freqs2^2)*L2)/((freqs1^2)-(freqs2^2));
C1_PC=((freqs1^2)*C1-(freqs2^2)*P2)/((freqs1^2)-(freqs2^2));
for iii = 1:size(Arc,1)
%single frequency smoothed code
 Rs_hat(Arc(iii,1):Arc(iii,2)) = HatchFilter( C1(Arc(iii,1):Arc(iii,2)),L1(Arc(iii,1):Arc(iii,2)),N);
% divergence-free carrier smoothed code
 Rd_hat(Arc(iii,1):Arc(iii,2)) = HatchFilter( C1(Arc(iii,1):Arc(iii,2)),L1_D(Arc(iii,1):Arc(iii,2)),N);
%  Ionosphere Free smoother
RI_hat(Arc(iii,1):Arc(iii,2)) = HatchFilter( C1_PC(Arc(iii,1):Arc(iii,2)),L1_LC(Arc(iii,1):Arc(iii,2)),N);
end
Rs_hat = Rs_hat';
Rd_hat = Rd_hat';
RI_hat = RI_hat';
rinexx.data(find(PRN==Sat),rinexx.col.Hatch_DF) = Rd_hat;
rinexx.data(find(PRN==Sat),rinexx.col.Hatch_IO_free) = RI_hat;
Rs_hat = [];
Rd_hat = [];
RI_hat = [];
end
end
function [ RH ] = HatchFilter( C1,L1,N)
RH = C1(1,:);
    for k=2:length(C1)
        if k<N
            n = k;
        else
            n = N;
        end
        if (C1(k-1)==0)
            RH(k) = C1(k);
        else
            RH(k) = (1/n)*C1(k) + ((n-1)/n)*(RH(k-1)+L1(k)-L1(k-1));
        end
    end
end
function [T,P,e,beta,lambda] = intpl(phii,D)
% phi = abs(rad2deg(phi));
phi=abs(phii);
T0 = [299.65, 294.15, 283.15, 272.15, 263.65]';
P0 = [1013.25, 1017.25, 1015.75, 1011.75, 1013.00]';
e0 = [26.31, 21.79, 11.66, 6.78, 4.11]';
beta0 = [6.3, 6.05, 5.58, 5.39, 4.53]'*10^(-3);
lambda0 = [2.77, 3.15, 2.57, 1.81, 1.55]';

dT = [0 , 7, 11, 15, 14.5]';
dP = [0 , -3.75, -2.25, -1.75, -0.5]';
de = [0 , 8.85, 7.24, 5.36, 3.39]';
dbeta = [0 , 0.25, 0.32, 0.81, 0.62]'*10^(-3);
dlambda = [0 , 0.33, 0.46, 0.74, 0.30]';

id = phi/15;
id_floor = floor(phi/15);
if mod(phi,15)==0
    T0_f = T0(id);
    P0_f = P0(id);
    e0_f = e0(id);
    beta0_f = beta0(id);
    lambda0_f = lambda0(id);
    
    dT_f = dT(id);
    dP_f = dP(id);
    de_f = de(id);
    dbeta_f = dbeta(id);
    dlambda_f = dlambda(id);
    
elseif mod(phi,15)~=0 && id_floor==0
    id = 1;
    T0_f = T0(id);
    P0_f = P0(id);
    e0_f = e0(id);
    beta0_f = beta0(id);
    lambda0_f = lambda0(id);
    
    dT_f = dT(id);
    dP_f = dP(id);
    de_f = de(id);
    dbeta_f = dbeta(id);
    dlambda_f = dlambda(id);

elseif mode(phi/15)~=0 && id_floor==5
    id = 5;
    T0_f = T0(id);
    P0_f = P0(id);
    e0_f = e0(id);
    beta0_f = beta0(id);
    lambda0_f = lambda0(id);
    
    dT_f = dT(id);
    dP_f = dP(id);
    de_f = de(id);
    dbeta_f = dbeta(id);
    dlambda_f = dlambda(id);

else
    
    x = [15*id_floor, 15*(id_floor+1)];
    y_T = [T0(id_floor), T0(id_floor+1)];
    T0_f = interp1(x,y_T,phi);
    y_P = [P0(id_floor), P0(id_floor+1)];
    P0_f = interp1(x,y_P,phi);
    y_e = [e0(id_floor), e0(id_floor+1)];
    e0_f = interp1(x,y_e,phi);
    y_beta = [beta0(id_floor), beta0(id_floor+1)];
    beta0_f = interp1(x,y_beta,phi);
    y_lambda = [lambda0(id_floor), lambda0(id_floor+1)];
    lambda0_f = interp1(x,y_lambda,phi);
    
    y_dT = [dT(id_floor), dT(id_floor+1)];
    dT_f = interp1(x,y_dT,phi);
    y_dP = [dP(id_floor), dP(id_floor+1)];
    dP_f = interp1(x,y_dP,phi);
    y_de = [de(id_floor), de(id_floor+1)];
    de_f = interp1(x,y_de,phi);
    y_dbeta = [dbeta(id_floor), dbeta(id_floor+1)];
    dbeta_f = interp1(x,y_dbeta,phi);
    y_dlambda = [dlambda(id_floor), dlambda(id_floor+1)];
    dlambda_f = interp1(x,y_dlambda,phi);

end

if phii<0
    D_min = 211;            
else
    D_min =28;              
end

T = T0_f - dT_f*cos(2*pi*(D-D_min)/365.25);
P = P0_f - dP_f*cos(2*pi*(D-D_min)/365.25);
e = e0_f- de_f*cos(2*pi*(D-D_min)/365.25);
beta = beta0_f- dbeta_f*cos(2*pi*(D-D_min)/365.25);
lambda = lambda0_f- dlambda_f*cos(2*pi*(D-D_min)/365.25);

end
function [ad bd cd aw bw cw]=interpp(phii,DOY)
x=[15 30 45 60 75];
%average dry
y1=[1.2769934e-3 1.2683230e-3 1.2465397e-3 1.2196049e-3 1.2045996e-3];
y2=[2.9153695e-3 2.9152299e-3 2.9288445e-3 2.9022565e-3 2.9024912e-3];
y3=[62.610505e-3 62.837393e-3 63.721774e-3 63.824265e-3 64.258455e-3];
%amplitude dry
y4=[0 1.2709626e-5 2.6523662e-5 3.4000452e-5 4.1202191e-5];
y5=[0 2.1414979e-5 3.0160779e-5 7.2562722e-5 11.723375e-5];
y6=[0 9.128400e-5 4.3497037e-5 84.795348e-5 170.37206e-5];
% wet
y7=[5.8021898e-4 5.6794847e-4 5.8118019e-4 5.9727542e-4 6.1641693e-4];
y8=[1.4275268e-3 1.5138625e-3 1.4572752e-3 1.5007428e-3 1.7599082e-3];
y9=[4.3472961e-2 4.6729510e-2 4.3908931e-2 4.4626982e-2 5.4736038e-2];
if phii<15
    ad=y1(1);
    bd=y2(1);
    cd=y3(1);
elseif phii>75
    ad=y1(5)-y4(5)*cos(2*pi*((DOY-28)/365.25));
    bd=y2(5)-y5(5)*cos(2*pi*((DOY-28)/365.25));
    cd=y3(5)-y6(5)*cos(2*pi*((DOY-28)/365.25));
else
ad=interp1(x,y1,phii)-interp1(x,y4,phii)*cos(2*pi*((DOY-28)/365.25));
bd=interp1(x,y2,phii)-interp1(x,y5,phii)*cos(2*pi*((DOY-28)/365.25));
cd=interp1(x,y3,phii)-interp1(x,y6,phii)*cos(2*pi*((DOY-28)/365.25));
end
if phii<15
    aw=y7(1);
    bw=y8(1);
    cw=y9(1);
elseif phii>75
    aw=y7(5);
    bw=y8(5);
    cw=y9(5);
else
    aw=interp1(x,y7,phii);
    bw=interp1(x,y8,phii);
    cw=interp1(x,y9,phii);
end
end
function [JD ] = cal2gpstime(varargin)

% Unpack
if nargin == 1
    cal_time = varargin{1};
    year = cal_time(:,1);
    month = cal_time(:,2);
    day = cal_time(:,3);
    hour = cal_time(:,4);
    min = cal_time(:,5);
    sec = cal_time(:,6);
    clear cal_time
else
    year = varargin{1};
    month = varargin{2};
    day = varargin{3};
    hour = varargin{4};
    min = varargin{5};
    sec = varargin{6};
end

% Seconds in one week
secs_per_week = 604800;

% Converts the two digit year to a four digit year.
% Two digit year represents a year in the range 1980-2079.
if (year >= 80 & year <= 99)
    year = 1900 + year;
end

if (year >= 0 & year <= 79)
    year = 2000 + year;
end

% Calculates the 'm' term used below from the given calendar month.
if (month <= 2)
    y = year - 1;
    m = month + 12;
end

if (month > 2)
    y = year;
    m = month;
end

% Computes the Julian date corresponding to the given calendar date.
JD = floor( (365.25 * y) ) + floor( (30.6001 * (m+1)) ) + ...
    day + ( (hour + min / 60 + sec / 3600) / 24 ) + 1720981.5;

% Computes the GPS week corresponding to the given calendar date.
gps_week = floor( (JD - 2444244.5) / 7 );

% Computes the GPS seconds corresponding to the given calendar date.
gps_seconds=round(((((JD-2444244.5)/7)-gps_week)*secs_per_week)/0.5)*0.5;
end