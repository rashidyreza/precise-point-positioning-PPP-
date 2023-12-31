
function [windup,Xctsun]=windup_ppp(rinex)
%% Compute Sun position              ***(KamPar)***
Leapsec = 16;

MJD    = rinex.data(:,rinex.col.JD); %% epo: Y M D h m s
Xctsun  = zeros(length(MJD),3);
Rsun    = zeros(length(MJD),1);
xctsat = rinex.data(:,rinex.col.Xsat:rinex.col.Zsat);
for ii=1:length(MJD)
    [X,xctsun,R,L,B] = SUN_bern(MJD(ii) - (Leapsec - (350.7512/1000) )/86400);
    Xctsun(ii,:) = xctsun';
    Rsun(ii)  = R;
end
A = Xctsun - xctsat;
B = sum(A.*A,2);
rcv = rinex.r0./norm(rinex.r0);
i = [rcv(1),0 0];
j = [0, rcv(2), 0];
k = [0 0 rcv(3)];
for ii=1:length(MJD)
E_sun_sat(ii,:) = A(ii,:)/B(ii);
U(ii,:) = rinex.r0' - xctsat(ii,:);
U(ii,:) = U(ii,:)/norm(U(ii,:));
K(ii,:) = - xctsat(ii,:)./norm(xctsat(ii,:));
J(ii,:) = cross(K(ii,:),E_sun_sat(ii,:));
I(ii,:) = cross(J(ii,:),K(ii,:));
Dp(ii,:) = [[I(ii,:)]' - [U(ii,:)]'*(U(ii,:)*[[I(ii,:)]']) - [cross(U(ii,:),J(ii,:))]']';
D(ii,:) = [[i]' - [U(ii,:)]'*(U(ii,:)*[[i]']) + [cross(U(ii,:),j)]']';
phi(ii,1) = sign(cross(Dp(ii,:),D(ii,:))*[U(ii,:)]')*acos( (Dp(ii,:)*[D(ii,:)]')/(norm(Dp(ii,:))*norm(D(ii,:)) ) );
end
freqs1=1575.42e6;   
freqs2=1227.60e6;
c = 299792458;
lambdaN = c/(freqs1+freqs2);
windup = lambdaN*phi;


end

function [X,Xct,R,L,B]=SUN_bern(XMJD)
%
%
% MATLAB version of bernese SUN.for
% SEPTEMBER 2004 Djamal Asgari
% NAME       :  SUN
%
%      CALL SUN(XMJD,X,R,L,B)
%
% PURPOSE    :  COMPUTATION OF POSITION OF THE SUN AT TIME XMJD
%               (MODIFIED JULIAN DATE). THIS SR WAS WRITTEN USING
%               SIMON NEWCOMB'S "TABLES OF THE SUN"
%
%               PRECISION  :        MAXIMUM        MEAN
%                                       DIFFERENCES
%                              L      .08"         .03"
%                              B      .03"         .005"
%                              R     3.E-7        .5E-7
%
% PARAMETERS :
%         IN :  XMJD   : EPOCH IN MODIFIED JULIAN DATE       R*8
%        OUT :  X(K),K=1,2,3 : RECTANGULAR COORDINATES OF    R*8
%                        THE SUN IN EQUATORIAL SYSTEM 1950.0
%               R      : DISTANCE EARTH-SUN (IN AU           R*8
%               L , B  : ECLIPTICAL LONGITUDE, LATITUDE IN   R*8
%                        SYSTEM OF EPOCH XMJD
%
% SR CALLED  :  PRENEW
%
% REMARKS    :  ---
%
% AUTHOR     :  G.BEUTLER
%
% VERSION    :  3.4  (JAN 93)
%
% CREATED    :  87/12/11 12:02        LAST MODIFIED :  90/01/11 11:45
%COPYRIGHT  :  ASTRONOMICAL INSTITUTE
%      1987      UNIVERSITY OF BERNE
%                    SWITZERLAND
%

%        REAL*8 XMJD,T,L,B,R,X(3),Y(3),PID,GD,EKL,TH,PRAEZ(3,3)
%        INTEGER*4 JME(4),IME(4),S1ME(4),K1ME(4),S2ME(4),K2ME(4),
%     1            JV(39),IV(39),S1V(39),K1V(39),S2V(39),K2V(39),
%     2            JM(45),IM(45),S1M(45),K1M(45),S2M(45),K2M(45),
%     3            JJ(21),IJ(21),S1J(21),K1J(21),S2J(21),K2J(21),
%     4            JS(11),IS(11),S1S(11),K1S(11),S2S(11),K2S(11)
%        INTEGER*4 JBV(22),IBV(22),SBV(22),KBV(22),
%     1            JBM(3),IBM(3),SBM(3),KBM(3),
%     2            JBJ(7),IBJ(7),SBJ(7),KBJ(7),
%     3            JBS(2),IBS(2),SBS(2),KBS(2)
JME=-ones(1,4);
IME=[1,2,3,4];
S1ME=[13,5,15,23];
K1ME=[243,225,357,326];
S2ME=[28,6,18,5];
K2ME=[335,130,267,239];
JV=[-ones(1,4),-2*-ones(1,5),-3*ones(1,5),-4*ones(1,5)  ,-5*ones(1,4) ,-6*ones(1,4)   ,-7*ones(1,4),-8*ones(1,5)  , -9*ones(1,2),-10];
IV=[0,1,2,3,0,1,2,3,4,2,3,4,5,6,3,4,5,6,7,5,6,7,8,...
    6,7,8,9,7,8,9,10,8,9,12,13,14,9,10,10];
S1V=[75,4838,74,9,3,116,5526,2497,44,13,666,1559,...
    1024,17,3,210,144,152,6,84,37,123,154,38,14,...
    10,14,20,6,3,0,11,0,42,0,32,6,0,3];
K1V=[2962,2991,2076,2490,1620,1488,1483,3159,3123,...
    1760,1777,3453,3182,3150,1980,2062,1954,3438,...
    3220,2356,2218,1953,3596,2641,2530,2300,...
    120,2940,2790,2880,0,3220,0,2592,0,488,...
    3510,0,180];
S2V=[94,2359,69,16,4,160,6842,869,52,21,1045,1497,194,...
    19,6,376,196,94,6,163,59,141,26,80,25,14,12,42,12,...
    4,4,24,6,44,12,33,13,4,8];
K2V=[2050,2091,3485,3300,900,584,583,2267,388,900,...
    876,2552,495,430,900,1163,1052,2548,590,1454,...
    1322,1054,2700,1743,1640,1350,2840,2035,1940,...
    1660,1350,2340,2180,1697,2220,1387,2610,2560,...
    2930];
JM=[ ones(1,3),2*ones(1,4),3*ones(1,4),4*ones(1,4),5*ones(1,4),6*ones(1,4),7*ones(1,3),8*ones(1,4),9*ones(1,3),...
    10*ones(1,3),11, 11,12,13, 13,15, 15,17 ,17];
IM=[2,1,0,3,2,1,0,4,3,2,1,4,3,2,1,5,4,3,2,6,5,4,3,...
    6,5,4,7,6,5,4,7,6,5,7,6,5,7,6,7,8,7,...
    9,8,10,9];
S1M=[6,273,48,41,2043,1770,28,4,129,425,8,34,500,585,...
    9,7,85,204,3,0,20,154,101,6,49,106,3,10,52,21,4,...
    28,62,5,19,5,17,44,6,13,45,21,0,4,26];
K1M=[2180,2177,2603,3460,3439,2004,1480,2840,2942,...
    3389,70,710,1052,3341,3250,1720,546,1008,180,...
    0,1860,2274,963,3010,1765,2227,720,3070,3489,...
    2152,570,2980,3460,680,1110,3380,590,1059,2320,...
    1840,2278,3090,0,2430,1130];
S2M=[8,150,28,52,2057,151,31,6,168,215,6,49,478,105,...
    10,12,107,89,3,5,30,139,27,10,60,38,5,15,45,8,...
    6,34,17,8,15,0,20,9,5,15,5,22,6,4,0];
K2M=[1300,1277,3470,2554,2538,2950,2343,1800,2035,...
    2490,900,3397,152,659,530,900,3246,110,1080,...
    2170,957,1373,1880,2090,862,1329,3490,2170,2597,...
    3100,3290,2081,2570,3370,230,0,3300,210,1430,...
    940,1430,2200,2610,1530,0];
JJ=[ones(1,5), 2*ones(1,4),3*ones(1,4),4*ones(1,4),5*ones(1,4)];
IJ=[-3,-2,-1,0,1,-3,-2,-1,0,-4,-3,-2,-1,-4,-3,-2,-1, -5,-4,-3,-2];
S1J=[3,163,7208,2600,73,69,2731,1610,73,5,164,556,210, 16,44,80,23,0,5,7,9];
K1J=[1980,1986,1795,2632,2763,808,871,1095,2526,1580, 1705,827,985,2590,1682,777,930,0,2590,1640,710];
S2J=[5,208,7067,244,80,103,4026,1459,8,9,281,803,174,29, 74,113,17,3,10,12,14];
K2J=[1120,1120,895,3386,65,3505,3571,195,2630,690,812,3526,86,1700,799,3477,30,2520,1690,760, 3430];
JS=[ones(1,4), 2*ones(1,4),3 3,4];
IS=[-2,-1,0,1,-3,-2,-1,0,-2,-1,-2];
S1S=[11,419,320,8,0,108,112,17,21,17,3];
K1S=[1050,1006,2695,2700,0,2906,2936,2770,2890, 2910,2880];
S2S=[15,429,8,8,3,162,112,0,32,17,4];
K2S=[110,106,3530,0,1980,2006,2031,0,2001,2010, 1940];
JBV=[-1*ones(1,4),-2*ones(1,4),-3*ones(1,5),-4,-4,-4,-5,-5,-6,-6,-6,-8];
IBV=[0,1,2,3,1,2,3,4,2,3,4,5,6,3,5,6,6,7, 5,7,8,12];
SBV=[29,5,92,7,23,12,67,14,14,8,210,7,4,6,31, 12,9,19,6,4,4,10];
KBV=[1450,3230,937,2620,1730,1490,1230,1110,2010,...
    1870,1518,1530,2960,2320,18,1800,270,180,2880,...
    570,570,610];
JBM=[2,2,4];
IBM=[-2,0,-3];
SBM=[8,8,7];
KBM=[900,3460,1880];
JBJ=[ones(1,4),2,3,3];
IBJ=[-2,-1,0,1,-1,-2,-1];
SBJ=[7,17,16,23,166,6,18];
KBJ=[1800,2730,1800,2680,2655,1710,2670];
JBS=[1 1];
IBS=[-1,1];
SBS=[6,6];
KBS=[2600,2800];

T=(XMJD-15019.5D0)/36525.D0;
TT=T;
L=279.D0+41.D0/60+48.04D0/3600;
L=L+(129602768.13D0.*T+1.089D0*T.^2)/3600;
L=mod(L/180*pi,2*pi);
GD=358.D0+28.D0/60+33.D0/3600;
GD=GD+(129596579.1D0*T-.54D0*T.^2-.012D0*T.^3)/3600;
GE=mod(GD/180*pi,2*pi);
TH=(XMJD+3242.297D0)/365.25D0;
GME=mod((248.07D0+1494.7235D0*TH)/180*pi,2*pi);
GV =mod((63.07037D0+22518.442986D0*T)/180*pi,2*pi);
GJ =mod((221.64742D0+32964.466939D0*T)/180*pi,2*pi);
GS =mod((193.13230D0+34777.259042D0*T)/180*pi,2*pi);
GM =mod((165.94905D0+16859.069667D0*T)/180*pi,2*pi);
% TT
% ones(length(TT))
XL1=6.4*sin((231.19+20.2*TT)./180*pi)...
    +(1.882-.016*TT)*sin((57.24+150.27*TT)/180*pi)...
    +.266*sin((31.8+119.0*TT)/180*pi) ...
    +.202*sin((315.6+893.3*TT)/180*pi);
L=L+XL1/3600*pi/180;
GE=GE+XL1/3600*pi/180;
GV=GV-XL1/3600/180*pi;
GJ=GJ+XL1/3600/180*pi;
GS=GS+XL1/3600/180*pi;
W1=299.1+(GV-GE)/pi*180;
GD=63.07037D0+22518.442986D0*T;
W2=90.+mod(GD,360.D0);
C=(6910.057-17.24*TT-.052*TT^2)*sin(GE)...
    +(72.338-.361*TT)*sin(2*GE)...
    +(1.054-.001*TT)*sin(3*GE)+.018*sin(4*GE);
R=3057-15*T+cos(GE)*(-727412.D0+1814*T+5*T^2)...
    +cos(2*GE)*(-9138+46*T)+cos(3*GE)*(-145+T)...
    +cos(4*GE)*(-2);
XLME=0;
RME=0;

for K=1:4
    W=-(JME(K)*GME+IME(K)*GE);
    W1=K1ME(K)/180.*pi;
    W2=K2ME(K)/180.*pi;
    XLME=XLME+S1ME(K)*cos(W+W1);
    RME=RME+S2ME(K)*cos(W+W2);
end

XLV=0;
RV=0;
for  K=1:39
    W=-(JV(K)*GV+IV(K)*GE);
    W1=K1V(K)/1800.*pi;
    W2=K2V(K)/1800.*pi;
    if(K==2)W1=299.1017/180*pi ;end
    if(K==7)W1=148.3133/180*pi ;end
    if(K==8)W1=315.9433/180*pi ;end
    if(K==12)W1=345.2533/180*pi ;end
    if(K==11)W1=177.71/180*pi ;end
    if(K==13)W1=318.12/180*pi ;end
    if(K==2)W2=209.08/180*pi ;end
    if(K==7)W2=58.3183/180*pi ;end
    if(K==11)W2=87.57/180*pi ;end
    if(K==12)W2=255.25/180*pi ;end
    if(K==16)W2=116.28/180*pi ;end
    WH=-JV(K)*330.9017/180*pi-(IV(K)+JV(K))*GE-JV(K)*(GV+pi);
    XLV=XLV+S1V(K)*cos(WH+W1) ;
    RV=RV+S2V(K)*cos(WH+W2);
end

XLM=0;
RM=0;
for K=1:45
    W=(-JM(K)*GM+IM(K)*GE);
    W1=K1M(K)/1800.*pi;
    W2=K2M(K)/1800.*pi;
    if(K==5)W1=343.8883/180*pi ;end
    if(K==6)W1=200.4017/180*pi ;end
    if(K==10)W1=338.88/180*pi ;end
    if(K==13)W1=105.18/180*pi ;end
    if(K==14)W1=334.05/180*pi ;end
    if(K==5)W2=253.8283/180*pi ;end
    if(K==13)W2=15.17/180*pi ;end
    WH=-JM(K)*127.0633/180*pi+(IM(K)-JM(K))*GE+JM(K)*GM;
    XLM=XLM+S1M(K)*cos(WH+W1);
    RM=RM+S2M(K)*cos(WH+W2);
end;

XLJ=0;
RJ=0;
for K=1:21
    W=-(JJ(K)*GJ+IJ(K)*GE);
    W1=K1J(K)/1800.*pi;
    W2=K2J(K)/1800.*pi;
    if(K==3)W1=179.5317/180*pi ;end
    if(K==4)W1=263.2167/180*pi ;end
    if(K==7)W1=87.145/180*pi ;end
    if(K==8)W1=109.4933/180*pi ;end
    if(K==12)W1=82.65/180*pi ;end
    if(K==3)W2=89.545/180*pi ;end
    if(K==7)W2=357.1083/180*pi ;end
    if(K==8)W2=19.4667/180*pi ;end
    if(K==12)W2=352.56/180*pi ;end
    WH=-JJ(K)*88.4450/180*pi-(JJ(K)+IJ(K))*GE+JJ(K)*GJ;
    XLJ=XLJ+S1J(K)*cos(WH+W1);
    RJ=RJ+S2J(K)*cos(WH+W2);
end

XLS=0;
RS=0;
for K=1:11
    W=-(JS(K)*GS+IS(K)*GE);
    W1=K1S(K)/1800.*pi;
    W2=K2S(K)/1800.*pi;
    if(K==2)W1=100.58/180*pi ;end
    if(K==3)W1=269.46/180*pi ;end
    WH=-JS(K)*10.2417/180*pi-(IS(K)+JS(K))*GE+JS(K)*GS;
    XLS=XLS+S1S(K)*cos(WH+W1);
    RS=RS+S2S(K)*cos(WH+W2);
end

XBV=0;
for K=1:22
    W=(KBV(K)/10.-JBV(K)*330.9017)/180*pi -(IBV(K)+JBV(K))*GE-JBV(K)*(GV+pi);
    XBV=XBV+SBV(K)*cos(W);
end

XBM=0;
for K=1:3
    W=(KBM(K)/10.-JBM(K)*127.0633)/180*pi...
        -(IBM(K)+JBM(K))*GE+JBM(K)*GM;
    XBM=XBM+SBM(K)*cos(W);
end

XBJ=0;
for K=1:7
    W=(KBJ(K)/10.-JBJ(K)*88.445)/180*pi...
        -(IBJ(K)+JBJ(K))*GE+JBJ(K)*GJ;
    XBJ=XBJ+SBJ(K)*cos(W);
end

XBS=0;
for K=1:2
    W=(KBS(K)/10.-JBS(K)*10.2417)/180*pi...
        -(IBS(K)+JBS(K))*GE+JBS(K)*GS;
    XBS=XBS+SBS(K)*cos(W);
end

GD=296.104608D0+477198.849108D0*T+.9192D-2*T^2+14.D-6*T^3;
XM=mod(GD/180*pi,2*pi);
GD=350.737486D0+445267.114217D0*T-.1436D-2*T^2+2.D-6*T^3;
D=mod(GD/180*pi,2*pi);
GD=11.D0+15.D0/60+3.2D0/3600 ...
    +(1739527290.54D0*T-11.56D0*T^2-.12D-2*T^3)/3600;
F=mod(GD/180*pi,2*pi);
GD=259+10.D0/60+59.79D0/3600-6962911.23D0*T/3600 ...
    +(7.48*T^2+.0086*T^3)/3600;
US=L-mod(GD*pi/180,2*pi);
XLMO=6.454*sin(D)+.013*sin(3*D)+.177*sin(D+XM) ...
    -.424*sin(D-XM)+.039*sin(3*D-XM)-.064*sin(D+GE) ...
    +.172*sin(D-GE)-.013*sin(D-XM-GE)-.013*sin(2*US);
RMO=1336*cos(D)+3*cos(3*D)+37*cos(D+XM)-133*cos(D-XM) ...
    +8*cos(3*D-XM)-14*cos(D+GE)+36*cos(D-GE) ...
    -3*cos(D-GE-XM)+3*cos(2*US);
BMO=.576*sin(F)+.016*sin(F+XM)-.047*sin(F-XM) ...
    +.021*sin(F-2*US)+.005*sin(F-2*US-XM) ...
    +.005*sin(F+GE)+.005*sin(F-GE);
L=L+(XLME+XLV+XLM+XLJ+XLS)/3600000./180*pi;
L=L+(C+XLMO)/3600*pi/180;
B=(XBV+XBM+XBJ+XBS)/3600000./180*pi;
B=-B+BMO/3600*pi/180;
R=(R+(RME+RV+RM+RJ+RS)/10+RMO)*1.D-8;
R=10.D0^R;
EKL=(23.D0+27.D0/60+8.26D0/3600);
EKL=EKL-(46.845*T+.59D-2*T^2-.181D-2*T^3)/3600;
EKL=EKL/180*pi;
X(1)=R*cos(L)*cos(B);
X(2)=R*sin(L)*cos(B);
X(3)=R*sin(B);
Y(1)=X(1);
Y(2)=X(2)*cos(EKL)-X(3)*sin(EKL);
Y(3)=X(2)*sin(EKL)+X(3)*cos(EKL);
PRAEZ=PRENEW_bern(XMJD);
for I=1:3
    X(I)=0.D0;
    for K=1:3;
        X(I)=X(I)+PRAEZ(K,I)*Y(K);
    end
end
R=R* 149597870691;

JD = XMJD +2400000.5;
T = (JD - 2451545.0) / 36525;
GST  =  280.46061837 + 360.98564736629*(JD-2451545.0) + 0.000387933*T*T - T*T*T/38710000.0;
GST =  GST -floor(GST /360)*360;
degrees2dms(GST);

R3G=[cos(GST*pi/180) sin(GST*pi/180) 0
    -sin(GST*pi/180) cos(GST*pi/180) 0
    0                      0         1 ];

Xct=R3G*X';
end
function PRAEZ=PRENEW_bern(XMJD)
%MATLAB version of bernese PRENEW.for
% 2004 Djamal Asgari
% NAME       :  PRENEW
%
%    CALL PRENEW(XMJD,PRAEZ)
%
% PURPOSE    :  COMPUTATION OF PRECESSION-MATRIX FROM 1950.0
%               TO EPOCH XMJD
%               SEE "USNO CIRCULAR NO 163",1981
%               (IAU RESOLUTIONS)
%
% PARAMETERS :
%         IN :  XMJD   : EPOCH IN MODIFIED JULIAN DATE       R*8
%        OUT :  PRAEZ  : PRECESSION - MATRIX                 R*8(3,3)
%
% SR CALLED  :  ---
%
% REMARKS    :  ---
%
% AUTHOR     :  G.BEUTLER
%
% VERSION    :  3.4  (JAN 93)
%
% CREATED    :  87/10/29 10:00        LAST MODIFIED :  88/11/21 17:14
%
% COPYRIGHT  :  ASTRONOMICAL INSTITUTE
%      1987      UNIVERSITY OF BERNE
%                    SWITZERLAND
%

%      IMPLICIT REAL*8 (A-H,O-Z)
%      REAL*8 PRAEZ(3,3)
%
%  TIME INTERVAL (IN JUL. CENTURIES) BETWEEN EPOCH AND 1950.0
TL=(XMJD-33281.9234D0)/36525 ;
%
%  TIME INTERVAL (IN JUL. CENTURIES) BETWEEN 1950.0 AND J2000.0
TU=(33281.9234D0-51544.5)/36525  ;
%
%  ROTATION ANGLES XI,Z AND THETA
XA=(2306.2181+1.39656*TU-0.000139*TU^2)*TL ...
    +(0.30188-0.000344*TU)*TL^2 ...
    + 0.017998*TL^3;
ZA=(2306.2181+1.39656*TU-0.000139*TU^2)*TL ...
    +(1.09468+0.000066*TU)*TL^2 ...
    + 0.018203*TL^3;
TA=(2004.3109-0.85330*TU-0.000217*TU^2)*TL ...
    -(0.42665+0.000217*TU)*TL^2 ...
    - 0.041833*TL^3;
XA=XA/206264.8;
ZA=ZA/206264.8;
TA=TA/206264.8;
COSXA=cos(XA);
SINXA=sin(XA);
COSZA=cos(ZA);
SINZA=sin(ZA);
COSTA=cos(TA);
SINTA=sin(TA);

%  ROTATION MATRIX
PRAEZ(1,1)= COSXA*COSZA*COSTA-SINXA*SINZA;
PRAEZ(2,1)= COSXA*SINZA*COSTA+SINXA*COSZA;
PRAEZ(3,1)= COSXA*SINTA;
PRAEZ(1,2)=-SINXA*COSZA*COSTA-COSXA*SINZA;
PRAEZ(2,2)=-SINXA*SINZA*COSTA+COSXA*COSZA;
PRAEZ(3,2)=-SINXA*SINTA;
PRAEZ(1,3)=-COSZA*SINTA;
PRAEZ(2,3)=-SINZA*SINTA;
PRAEZ(3,3)= COSTA;
end


