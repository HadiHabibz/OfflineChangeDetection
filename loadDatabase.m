% Hadi Habibzadeh 04/11/2018
%% Description and Examples:
% This function loads all files in the database
% It either uses a default path or asks the user
% to select the directory of the database, depending
% on the value of its input. For this function to 
% work, the files in the database must be named 
% systematically. The default is S1.mat, S2.mat, ...,
% S35.mat. Change the parameters in the beginning of
% the function according to your file names
% ------------------------------------------------------
% inputs:
% string: if string is 'default path', the functions looks
% for the default directory to load all files. This mode
% is good for debugging. Otherwise, the functions asks the
% user to select the directory of the database using a GUI.
% ------------------------------------------------------
% Example (1):
% subjectArray = loadDatabase( 'default path' );
% This loads all files in the data path and returns an array 
% of structures, where each element holds one file (subjectArray)
%
% Example (2):
% subjectArray = loadDatabase( 'GUI interface' );
% A gui pops up, asking user to enter the location of the
% database
%% 
function data = loadDatabase( string, subjectsNumber )

% Set these parameters according to your filenames
% for example, if files are named S1.mat, S2.mat, ...
% Set the filenameBeginsWith = 'S' and filenameTypeExtension
% = '.mat'. numberOfFiles shows the number of files in the
% database.
filenameBeginsWith = 'S';
filenameTypeExtension = '.mat';

if( nargin < 2 )
numberOfFiles = 35;
subjectsNumber = 1:35;

else 
    numberOfFiles = length( subjectsNumber );
end

% Check if user wants a gui or not
% Change the address below if database is in another
% directory
if( strcmp( string, 'default path' ) == true )
    path = 'D:\Research\Datasets\35 Subject Public Dataset';
    
else
    path = uigetdir();
    
end % if


% Load files one by one. Construct the filenames
% according to the names identified in the beginning
% of the function
parfor i = 1:numberOfFiles
    currentPath = path;
    filename = filenameBeginsWith;
    filename = strcat( filename, int2str( subjectsNumber( i ) ) );
    filename = strcat( filename, filenameTypeExtension );
    currentPath = strcat( currentPath, '\' );
    currentPath = strcat( currentPath, filename );
    data( i ) = load( currentPath );
    fprintf( 'file %s loaded successfully!\n', currentPath );
end % for ( i )