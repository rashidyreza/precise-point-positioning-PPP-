function [rinex] = read_rinex_obs(fname,nlines)
if (nargin < 2)
    nlines = 1e6;
end
% Initialize variables
rinex_data = [];
line_count = 1;
GNSSType = char([]);
GNSS = char([]);
% Read header
[ fid, rec_xyz, observables ] = read_rinex_header(fname);
num_obs = length(observables);

% Get the first line of the observations.
current_line = fgetl(fid);

% If not at the end of the file, search for the desired information.
while current_line ~= -1 & line_count < nlines
    
    % Get the time for this data epoch.
    current_time = [ str2num(current_line(2:3)) ; str2num(current_line(5:6)) ; ...
        str2num(current_line(8:9)) ; str2num(current_line(11:12)) ; ...
        str2num(current_line(14:15)) ; str2num(current_line(17:27)) ]';
    
    % How many SV's are there?
    current_num_sv = str2num(current_line(30:32));
    
    if current_num_sv>12
        ii=ceil(current_num_sv/12);
        for iii=2:ii
            % Get the next line.
            current_line2 = fgetl(fid);
            line_count = line_count + 1;
            current_line2 = current_line2(33:end);
            current_line = [current_line,current_line2];
        end
    end
    % Figure out which PRN's there are.
    for ii=1:current_num_sv
        GNSS(ii) = current_line(30 + 3*ii);
        current_prn(ii) = str2num(current_line(31 + 3*ii : 32 + 3*ii));
    end
    GNSSType = [GNSSType,GNSS];
    GNSS = char([]);
    % Get the data for all SV's in this epoch.
    for ii=1:current_num_sv
        
        % Get the next line.
        current_line = fgetl(fid);
        line_count = line_count + 1;
        
        % Check the length of the line and pad it with zeros to
        % make sure it is 80 characters long.
        if length(current_line) < 80
            add_spaces = 80 - length(current_line);
            
            for j = 1 : add_spaces
                current_line = [ current_line , '0' ];
            end
        end
        
        % Check if there are any blanks in the data and put a zero there.
        current_line = strrep(current_line,' ', '0');
        
        % Get the observables on this line.
        current_obs = [ str2num(current_line(1:14)) ; str2num(current_line(17:30)) ; ...
            str2num(current_line(33:46)) ; str2num(current_line(49:62)) ; str2num(current_line(65:78)) ];
        
        % If there are > 5 observables, read another line to get the rest of the observables for this SV.
        if num_obs > 5
            % Get the next line.
            current_line = fgetl(fid);
            line_count = line_count + 1;
            
            % Check the length of the line and pad it with zeros to
            % make sure it is 80 characters long.
            if length(current_line) < 80
                add_spaces = 80 - length(current_line);
                
                for j = 1 : add_spaces
                    current_line = [ current_line , '0' ];
                end
            end
            
            % Check if there are any blanks in the data and put a zero there.
            current_line = strrep(current_line,' ', '0');
            % Append the data in this line to the data from previous line.
            current_obs = [ current_obs ; str2num(current_line(1:14)) ; ...
                str2num(current_line(17:30)) ; str2num(current_line(33:46)) ; ...
                str2num(current_line(49:62)) ; str2num(current_line(65:78)) ];
            
        end  % if num_obs > 5
        
        % Format the data for this PRN as Date/Time, PRN, Observations.
        current_data = [ current_time , current_prn(ii) , current_obs' ];
                % Keep only data for the specified PRNs
        if nargin == 3 & PRN_list & isempty(find(PRN_list == current_prn(ii)))
            continue
        end    
        %Append to the master rinex data file.
        rinex_data = [ rinex_data ; current_data ];
        
    end  % for ii=1:current_num_sv
    % Get the next line.
    current_line = fgetl(fid);
    line_count = line_count + 1;
    
end  % while current_line ~= -1

% Convert time format
[ gpswk, gpssec ] = cal2gpstime(rinex_data(:,1:6));
rinex.data = [ gpswk gpssec rinex_data(:, 7:end) ];
% Define columns
rinex = define_cols(rinex, observables);
% Convert CP to meters
rinex = convert_rinex_CP(rinex);
rinex.GNSStype = GNSSType;
rinex.data(:,14) = GNSSType;
rinex.data(:,15) = rinex.data(:,14)*0;
rinex.data(:,16) = rinex.data(:,14)*0;
rinex.data(:,17) = rinex.data(:,14)*0;
rinex.data(:,18:21) = rinex.data(:,14:17)*0;
rinex.data(:,22:26) = rinex.data(:,16:20)*0;
rinex.data(:,27:31) = rinex.data(:,14:18)*0;
rinex.col.GNSStype = 14;
rinex.col.dtsat = 15;
rinex.col.dt_sat_rcv = 16;
rinex.col.Temission = 17;
rinex.col.Xsat = 18;
rinex.col.Ysat = 19;
rinex.col.Zsat = 20;
rinex.col.dtrcv = 21;
rinex.col.Azimuth = 22;
rinex.col.Elevation = 23;
rinex.col.delta_ro_rel = 24;
rinex.col.Tr0 = 25;
rinex.col.I1_klobuchar = 26;
rinex.col.Hatch_DF = 27;
rinex.col.Hatch_IO_free = 28;
rinex.col.TGD = 29;
rinex.col.Mwet = 30;
rinex.col.JD = 31;
rinex.r0 = rec_xyz;

% =========================================================================
function rinex = define_cols(rinex, observables)

% Defaults
rinex.col.WEEK = 1;
rinex.col.TOW = 2;
rinex.col.PRN = 3;

col_offset = 3;

for ii=1:length(observables)
    switch observables{ii}
        case {'L1'}
            rinex.col.L1 = ii + col_offset;
        case {'L2'}
            rinex.col.L2 = ii + col_offset;
        case {'P1'}
            rinex.col.P1 = ii + col_offset;
        case {'P2'}
            rinex.col.P2 = ii + col_offset;
        case {'C1'}
            rinex.col.C1 = ii + col_offset;
        case {'C2'}
            rinex.col.C2 = ii + col_offset;
        case {'D1'}
            rinex.col.D1 = ii + col_offset;
        case {'S1'}
            rinex.col.S1 = ii + col_offset;
        case {'S2'}
            rinex.col.S2 = ii + col_offset;
    end  % switch
end  % for ii=1:length(observables)
% =========================================================================
function [ rinex ] = convert_rinex_CP(rinex)

% GENERAL CONSTANTS
% =========================================================================
c = 299792458;          %----> Speed of light (meters/s).
Re = 6378137 ;          %----> Earth Radius (meters)
% =========================================================================
% CONVERSION FACTORS
% =========================================================================
Hz2MHz = 1E-6;
MHz2Hz = 1E6;
s2ns = 1E9;
ns2s = 1E-9;
s2micros = 1E6;
micros2s = 1E-6;
s2ms = 1E3;
ms2s = 1E-3;
dtr = pi / 180;
rtd = 180 / pi;
m2cm = 100;
cm2m = 1 / 100;
m2mm = 1000;
mm2m = 1 / 1000;
ft2m = 0.3048;  % Source: http://www.nodc.noaa.gov/dsdt/ucg/
m2ft = 1 / 0.3048;
ns2m = c * ns2s;  % Converts time in nano-sec to distance,
% assuming the speed of light.
% =========================================================================
% GNSS SPECIFIC CONSTANTS
% =========================================================================
L1 = 1575.42e6;         %----> Freqs in Hz.
L2 = 1227.60e6;
L5 = 1176.45e6;

L1MHz = 1575.42;        %----> Freqs in MHz.
L2MHz = 1227.60;
L5MHz = 1176.45;

L1GHz = 1.57542;        %----> Freqs in GHz.
L2GHz = 1.22760;
L5GHz = 1.17645;

LAMBDA_L1 = c / L1;     %----> Wavelengths in meters.
LAMBDA_L2 = c / L2;
LAMBDA_L5 = c / L5;

CA_CODE_RATE = 1.023e6; %----> C/A and P code chipping rate in chips/s.
P_CODE_RATE = 10.23e6;

CA_CHIP_PERIOD = 1 / CA_CODE_RATE;   %----> C/A & P code chip periods in s.
P_CHIP_PERIOD = 1 / P_CODE_RATE;

CA_CHIP_LENGTH = c / CA_CODE_RATE;  %----> C/A & P code chip lengths in meters.
P_CHIP_LENGTH = c / P_CODE_RATE;

CA_CODE_LENGTH = 1023;  % chips
% =========================================================================

if rinex.col.L1 ~= 0
    rinex.data(:, rinex.col.L1) = rinex.data(:, rinex.col.L1) * LAMBDA_L1;
end
if rinex.col.L2 ~= 0
    rinex.data(:, rinex.col.L2) = rinex.data(:, rinex.col.L2) * LAMBDA_L2;
end
% =========================================================================
function [ fid, rec_xyz, observables ] = read_rinex_header( file_name )
% Initialize the observables variable.
observables={};

% Assign a file ID and open the given header file.
fid=fopen(file_name);

% If the file does not exist, scream bloody murder!
if fid == -1
    display('Error!  Header file does not exist.');
else
    % Set up a flag for when the header file is done.
    end_of_header=0;
    
    % Get the first line of the file.
    current_line = fgetl(fid);
    
    % If not at the end of the file, search for the desired information.
    while end_of_header ~= 1
        % Search for the approximate receiver location line.
        if strfind(current_line,'APPROX POSITION XYZ')
            
            % Read xyz coordinates into a matrix.
            [rec_xyz] = sscanf(current_line,'%f');
        end
        
        % Search for the number/types of observables line.
        if strfind(current_line,'# / TYPES OF OBSERV')
            % Read the non-white space characters into a temp variable.
            [observables_temp] = sscanf(current_line,'%s');
            
            % Read the number of observables space and then create
            % a matrix containing the actual observables.
            for ii = 1:str2num(observables_temp(1))
                observables{ii} = observables_temp( 2*ii : 2*ii+1 );
            end
        end
        
        % Get the next line of the header file.
        current_line = fgetl(fid);
        
        %Check if this line is at the end of the header file.
        if strfind(current_line,'END OF HEADER')
            end_of_header=1;
        end
    end
end
% =========================================================================
function [ gps_week, gps_seconds ] = cal2gpstime(varargin)

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