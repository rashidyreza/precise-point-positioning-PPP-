function [SP3_XYZ] = ReadSP3(filename)
field = fopen(filename);
if field == -1
    disp(['The file ''' filename ''' does not exist.']);
    return;
end
current_line=fgetl(field);
a={};
EOH=0;
k=1;
while ischar(current_line)
   a(k,1) = cellstr(current_line);
   EOH = findstr(current_line,'ORB:CMB CLK:CMB');
 if EOH>1
    eoh=k;
 end
   current_line=fgetl(field);
   k=k+1;
end
 nos=eoh+1;
 n=1;
 prn_num = str2num(a{3,1}(5:6));
 while nos < size(a,1)
     Day = str2num(a{nos}(4:13));
     time=str2num(a{nos}(15:31));
     Date = [Day time]
     for i=1:prn_num
        num_s = str2num(a{nos+i,1}(3:4));
        x = str2num(a{nos+i,1}(6:18))*1000;
        y = str2num(a{nos+i,1}(20:32))*1000;
        z = str2num(a{nos+i,1}(34:46))*1000;   
        errort = str2num(a{nos+i,1}(48:60));
        [ gps_week, gps_seconds , JD ] = cal2gpstime(Date);
        Data(n,:) = [num_s , Date , x, y , z , errort ,  gps_week, gps_seconds , JD ];
        n = n+1;
     end
     nos=nos+prn_num+1;
 end
SP3_XYZ.data = Data;
SP3_XYZ.col.PRN = 1;
SP3_XYZ.col.X = 8;
SP3_XYZ.col.Y = 9;
SP3_XYZ.col.Z = 10;
SP3_XYZ.col.dtsat = 11;
SP3_XYZ.col.gps_week = 12;
SP3_XYZ.col.gps_seconds = 13;
SP3_XYZ.col.JD = 14;
function [ gps_week, gps_seconds , JD ] = cal2gpstime(varargin)

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