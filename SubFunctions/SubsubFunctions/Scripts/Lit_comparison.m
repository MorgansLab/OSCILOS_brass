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

i = 3;

TrumpetData = [83.5 234.2 353.1 469.1 591.9 702.7 812.5];
TrumpetCalc = [83.0 233.0 351.5 469.7 588.3 699.6 809.8];

TromboneData = [38.3 113.2 174.0 233.2 296.0 350.9 410.1];
TromboneCalc = [38.2 111.5 171.7 230.1 293.8 353.4 408.6];

HornData = [23.2 64.8 103.7 144.7 181.9 220.3 256.2];
HornCalc = [24.1 65.0 102.9 144.1 180.7 219.1 253.8];

TubaData = [42.8364688856729, 118.37916063675831, 174.8191027496382, 239.07380607814764, 298.1186685962374, 360.63675832127353, 419.68162083936335, 482.1997105643995, 541.2445730824893, 600.289435600579, 659.3342981186688, 718.3791606367583, 775.6874095513749, 832.1273516642548, 893.7771345875544];
WagnerTubaData = [66.76741575345198, 114.24664978614604, 176.5883948743579, 233.9794337627144, 287.36418088503007, 351.28038304860837, 410.2056082455518, 466.39501941549287, 528.6198493743445, 585.6770306155281, 645.280363562752, 709.2978921717759, 769.4292917866055];
CornophoneData = [56.644070889128585, 115.77790823126227, 173.7432173002475, 230.0312956146474, 293.26093164942563, 348.3245541664023, 406.2196563432125, 459.523941740179, 521.0781319413594, 578.9339658564448, 641.7435504854986, 705.3872881893762, 767.4055578473059, 827.7739655391257, 885.6553833216985, 942.7115726343848];
m=[];
ff=[];
[m(:,1), ff(1)] = EFP_calc(HornData);
[m(:,2), ff(2)] = EFP_calc(HornCalc);

f = figure('Name', strrep(geoms{i},'.txt','')+" comparison to results in literature");
set(f, 'Position', [280 150 600 800])
colours = get(gca, 'colororder');
hold on

plot(m(2:end,1),2:length(m(:,1)),lws(1)+mks(2),'Color','k','LineWidth',1.5,'MarkerSize',8);
plot(m(2:end,2),2:length(m(:,2)),lws(3)+mks(3),'Color','k','LineWidth',1.5,'MarkerSize',8);

run_name = "BCs1 "+strrep(geoms{i},'.txt','');
fid = fopen(strcat('./Outputs/',run_name,'/Results/Eigenvalues.txt'));
T = textscan(fid,"%[^\r\n]",1);
A = textscan(fid,"%d %f %f %f %*[^\r\n ]");
fclose(fid);
EFP{i} = A{4};
plot(EFP{i}(2:end),2:length(EFP{i}),lws(1)+mks(1),'Color',colours(2,:),'LineWidth',1.5,'MarkerSize',8);
set(gca,'FontName','Helvetica','FontSize',font_size,'LineWidth',1)
XLIM = xlim;
xlim([-1,1].*max(abs(xlim)));
ylim([2,8])

lgd = legend(["Measurements [1]", "Simulation [1]", "OSCILOS$$_{brass}$$"],'Interpreter','latex', 'Location', 'southoutside');

xlabel("EFP (cents)", 'FontSize', font_size, 'Interpreter', 'latex')
ylabel("Mode number", 'FontSize', font_size, 'Interpreter', 'latex')
grid on
set(gca, 'GridLineStyle', '--')
saveas(gca, strrep(geoms{i},'.txt','.png'));
