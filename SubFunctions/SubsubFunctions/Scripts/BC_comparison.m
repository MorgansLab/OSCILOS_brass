%% Explore BCs for converged geometries

% From Kemp, lips give roughly 0.3/(2200 * 2pi) s of time delay
% And approx 2% of attenuation across freq. spectrum

% Set p = 101325.0, T = 293.2, choice_M1_u1 = 2, M1_u1 = 0.5, min_GR = -400
% and all plots off

addpath('./SubFunctions/SubsubFunctions')

geoms = ["BilbConvTrump.txt", "BilbConvTromb.txt", "BilbConvHorn.txt", "NormConvTuba.txt"];
maxfreq = [2000, 1000, 600, 1000];

BCs = { {'inlet_type', 2, 'outlet_type', 1}, ...
        {'inlet_type', 2, 'outlet_type', 11}, ...
        {'inlet_type', 4, 'inlet_param1', 0.99, 'inlet_param2', -0.3*1000/(2200*2*pi), 'outlet_type', 11}, ...
        {'inlet_type', 13, 'inlet_param1', 0.2, 'outlet_type', 11} };
    
lws = ["-","--",":","-."];
mks = ["o","s","^","*"];

cls = [1,2,3,5];

EigVals = {};
EFP = {};

plts = [];

% for i = 1:length(geoms)
%     for j = 1:length(BCs)
%         run_name = sprintf("BCs%d ",j)+strrep(geoms{i},'.txt','');
%         EigVals{i,j} = OSCILOS_brass('geom_filename',geoms{i}, 'run_name',run_name, BCs{j}{:}, 'max_freq', maxfreq(i));
%         [EFP{i,j},~] = EFP_calc(abs(imag(EigVals{i,j}))./2./pi);
%         close all
%     end
% end
% 
% return

for i = 1:length(geoms)
    f = figure('Name', strrep(geoms{i},'.txt','')+" EFP comparison for 4 BCs");
    set(f, 'Position', [280 150 600 800])
    colours = get(gca, 'colororder');
    hold on
    for j = 1:length(BCs)
        run_name = sprintf("BCs%d ",j)+strrep(geoms{i},'.txt','');
        fid = fopen(strcat('./Outputs/',run_name,'/Results/Eigenvalues.txt'));
        T = textscan(fid,"%[^\r\n]",1);
        A = textscan(fid,"%d %f %f %f %*[^\r\n ]");
        fclose(fid);
        EFP{i,j} = A{4};
        plts(i,j) = plot(EFP{i,j}(2:end),2:length(EFP{i,j}),lws(j)+mks(j),'Color',colours(cls(j),:),'LineWidth',1.5,'MarkerSize',8);
    end
    set(gca,'FontName','Helvetica','FontSize',font_size,'LineWidth',1)
    XLIM = xlim;
    xlim([-1,1].*max(abs(xlim)));
    %ylim([2,length(M{1})])

    lgd = legend(["Closed - Open","Closed - Levine-Schwinger","Time Lag - Levine Schwinger", "Piston - Levine Schwinger"],'Interpreter','none', 'Location', 'southoutside');

    xlabel("EFP (cents)", 'FontSize', font_size, 'Interpreter', 'latex')
    ylabel("Mode number", 'FontSize', font_size, 'Interpreter', 'latex')
    grid on
    set(gca, 'GridLineStyle', '--')
    saveas(gca, strrep(geoms{i},'.txt','.png'));
end
