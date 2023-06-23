% Write the file names into the Excel
close all; clear; clc

% choose the experimental date
exp_date = '20230405';

% get the names
Result_folder = dir(['F:\Processing & Results\FSI - Actin &  Individual Obstacle' ...
    '\',exp_date,'-Actin-Individual_triangularPillar_uppoint\results']);
count = 1;

Files = dir(fullfile(Result_folder(1).folder, 'tra*.mat'));

for jj = 1:length(Files)

    the_name = Files(jj).name(1:end-4);
    if ~contains(the_name,'Info')
        Names_tobeWritten{count, 1} = the_name;
        count = count + 1;
    end
end

% choose the Excel and write the names
excelpathname = ['F:\Processing & Results\FSI - Actin &  Individual Obstacle' ...
    '\Exp data triangular pillar lateral pointing\'];
excelname = [exp_date,'-Actin-Individual_triangularPillar_uppoint.xlsx'];
for kk = 1: length(Names_tobeWritten)

    Loc = ['A', num2str(kk+1)];
    writematrix(Names_tobeWritten{kk, 1},[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);

end