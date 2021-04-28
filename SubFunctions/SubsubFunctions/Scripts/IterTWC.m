%% Script to run TWC (Tuba, Wagner Tuba, Cornophone) convergence study

for i = 1:40
    if i<21
        nm = 'Tuba';
    elseif i<41
        nm = 'WagnerTuba';
    else
        nm = 'Cornophone';
    end
    OSCILOS_brass('geom_filename', sprintf('TWC/20TWC%d.txt',i),...
        'run_name', sprintf('Norm%s%d',nm,i)...
    );
    close all
end