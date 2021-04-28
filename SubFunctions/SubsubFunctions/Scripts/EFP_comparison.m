%% EFP comparison plots

% Run without clearing variables from OSCILOS to generate graph from most
% recent run

%clear all
close all
addpath('./SubFunctions/SubsubFunctions')

%projects = {'Trumpet1','Trumpet2','Trumpet3'};
%projects = {'TrumpetBC2','TrumpetBC13','TrumpetBC14','TrumpetBC2LowFlow','TrumpetBC2HighFlow'};

ITER_NUM = 7;
ITER_MIN = -0.003;
ITER_MAX = 0.003;
ITER_PARAM = 'inlet_param2';
GEOM_NAME = '20PBH.txt';
ITER_NAME = '_';

font_size = 16;

iter_vals = linspace(ITER_MIN, ITER_MAX, ITER_NUM);
projects={};

for iter = 1 : ITER_NUM
    %projects{iter} = sprintf("%s %s %d-%d",ITER_NAME,ITER_PARAM,floor(iter_vals(iter)),floor(100*(iter_vals(iter)-floor(iter_vals(iter)))));
    %projects{iter} = sprintf("%s%d",ITER_NAME,iter);
    projects{iter} = strrep(sprintf("%s %s %s",ITER_NAME,ITER_PARAM,num2str(iter_vals(iter))),'.','-');
    %projects{iter} = strrep(sprintf("%s %s %s",ITER_NAME,ITER_PARAM,GEOM_NAME),'.',sprintf("%d-",iter_vals(iter)));
end

N = length(projects);
m = [];
ff = [];
F = [];
M = {};

for i = 1:N
    fid = fopen(strcat('./Outputs/',projects{i},'/Results/Eigenvalues.txt'));
    T = textscan(fid,"%[^\r\n]",1);
    A = textscan(fid,"%d %f %f %f");
    fclose(fid);
    f = A{2};
    [a, F(i)] = EFP_calc(f);
    M{i} = a;
end

TrumpetData = [83.5 234.2 353.1 469.1 591.9 702.7 812.5];
TrumpetCalc = [83.0 233.0 351.5 469.7 588.3 699.6 809.8];

% [m(:,1), ff(1)] = EFP_calc(TrumpetData);
% [m(:,2), ff(2)] = EFP_calc(TrumpetCalc);

TromboneData = [38.3 113.2 174.0 233.2 296.0 350.9 410.1];
TromboneCalc = [38.2 111.5 171.7 230.1 293.8 353.4 408.6];

[m(:,1), ff(1)] = EFP_calc(TromboneData);
[m(:,2), ff(2)] = EFP_calc(TromboneCalc);

HornData = [23.2 64.8 103.7 144.7 181.9 220.3 256.2];
HornCalc = [24.1 65.0 102.9 144.1 180.7 219.1 253.8];

% [m(:,1), ff(1)] = EFP_calc(HornData);
% [m(:,2), ff(2)] = EFP_calc(HornCalc);

fig = figure('Name', 'Equivalent Fundamental Pitch deviation');
set(fig, 'Position', [280 150 1000 800])

subplot(1,2,1)
set(gca, 'Position', [0.102,0.11,0.434,0.815]);
hold on
plot(m(2:end,1)',2:size(m,1),'-xk','linewidth',1, 'MarkerSize', 6, 'DisplayName', 'Measurements')
plot(m(2:end,2)',2:size(m,1),':^k','linewidth',1, 'MarkerSize', 6, 'DisplayName', 'Bilbao Calc')

plnms = strings(1,N);
ccc = flipud(parula(length(plnms)));
for i = 1:N
    plnms(i) = sprintf("%s %d",ITER_NAME,iter_vals(i));
    plot(M{i}(2:end)',2:length(M{1}),'-o','linewidth',1, 'MarkerSize', 6, 'DisplayName', plnms(i), 'Color', ccc(i,:))
end
    
set(gca,'YColor','k','Box','on');
set(gca,'FontName','Helvetica','FontSize',font_size,'LineWidth',1)
XLIM = xlim;
xlim([-1,1].*max(abs(xlim)));
ylim([2,length(M{1})])

lgd = legend('Interpreter','none', 'Location', 'best');
lgd.BoxFace.ColorType='truecoloralpha';
lgd.BoxFace.ColorData=uint8(255*[1 1 1 0.8]');

xlabel("EFP (cents)", 'FontSize', font_size, 'Interpreter', 'latex')
ylabel("Mode number", 'FontSize', font_size, 'Interpreter', 'latex')
grid on
set(gca, 'GridLineStyle', '--')

subplot(1,2,2)
set(gca, 'Position', [0.678,0.170851772125135,0.227,0.754148227874865]);
FF = 1200/log(2).*log([ff, F]./ff(1));
X = {'Bilbao and Chick Meas', 'Bilbao and Chick Calc'};
X = categorical([X, plnms],[X, plnms]);
cc = [[0 0 0]', 0.5*[1 1 1]']';
hold on
for i = 1:length(X)
    bb = bar(X(i),FF(i));
    if i<=2
        bb.FaceColor = cc(i,:);
    else
        bb.FaceColor = ccc(i-2,:);
    end
end
set(gca,'YColor','k','Box','on');
set(gca,'FontName','Helvetica','FontSize',font_size,'LineWidth',1)
%xlabel("Curve index", 'FontSize', font_size, 'Interpreter', 'latex')
ylabel("Reference pitch deviation (cents)", 'FontSize', font_size, 'Interpreter', 'none')
ax = gca;
ax.XAxis.FontSize = font_size-6;
ax.XAxis.TickLabelRotation = 45;


