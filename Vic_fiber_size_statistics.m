clear; close all; clc;

[filename, pathname]=uigetfile({'G:\PhD, PMMH, ESPCI\Experimental Data (EXTRACTED)\20220614-SU8_Fibers_Size\*.tif'}, 'Choose a file to be processed');  % input file

% path for the result files
pathout = uigetdir('G:\PhD, PMMH, ESPCI\Processing\20220614-SU8_Fibers_Size\results\', 'Choose the saving folder');
[status,msg,msgID] = mkdir(pathout);

% path of the experiment
basepath=pathname;
% name of the file to read
tifname=filename;