function out = clean_csv_read(file)
    % Load data from .csv
    % Remove quotes
    fid = fopen('output/temp.csv','w');
    s = fileread(file);
    s = strrep(s,'"','');
    fprintf(fid,s);
    fclose(fid);
    out = csvread('output/temp.csv');