function matfilename = cond2mat(rawfilename)

    [p n e] = fileparts(rawfilename);
    matfilename = strcat(n, '.mat');
    
    run(fullfile(p, n));
    save(matfilename, 'names', 'durations', 'onsets');

end