%% Subsubfc to write to logfile

function write_log(logstring, logfname, mode)
    if nargin<3
        mode = 'a';
    end
    logfid = fopen(logfname, mode);
    % using + instead of %s for logstring because then fprintf parses \t characters
    fprintf(logfid,"%s  "+logstring+"\r\n",datestr(now,'dd-mm-yyyy HH:MM:SS.FFF'));
    fclose(logfid);
end