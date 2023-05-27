% Write the file names into the Excel
close all; clear; clc

% choose the experimental date
exp_date = '20230313';

% get the names
Result_Groups = dir(['F:\Processing & Results\Actin Filaments in Porous Media' ...
    '\',exp_date,'-Actin\results\G*']);
count = 1;
for ii = 1:length(Result_Groups)

    Files = dir(fullfile(Result_Groups(1).folder, Result_Groups(ii).name, 'tra*.mat'));

    for jj = 1:length(Files)

        the_name = Files(jj).name(1:end-4);
        if ~contains(the_name,'Info')
            Names_tobeWritten{count, 1} = the_name;
            count = count + 1;
        end
    end    
end

% choose the Excel and write the names
excelpathname = 'F:\Processing & Results\Actin Filaments in Porous Media\Dynamics manually classification\';
excelname = ['Dynamics ',exp_date,'-Actin.xlsx'];
for kk = 1: length(Names_tobeWritten)

    Loc = ['A', num2str(kk)];
    writematrix(Names_tobeWritten{kk, 1},[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);

end