%% EFP comparison plots

% Run without clearing variables from OSCILOS to generate graph from most
% recent run

%clear all
close all
addpath('./SubFunctions/SubsubFunctions')

PROJ_NAME = 'NormCornophone';

font_size = 16;

projects={};

for iter = 1 : 20
    projects{iter} = sprintf("%s%d",PROJ_NAME,iter+40);
end

TubaData = [42.8364688856729, 118.37916063675831, 174.8191027496382, 239.07380607814764, 298.1186685962374, 360.63675832127353, 419.68162083936335, 482.1997105643995, 541.2445730824893, 600.289435600579, 659.3342981186688, 718.3791606367583, 775.6874095513749, 832.1273516642548, 893.7771345875544];
WagnerTubaData = [66.76741575345198, 114.24664978614604, 176.5883948743579, 233.9794337627144, 287.36418088503007, 351.28038304860837, 410.2056082455518, 466.39501941549287, 528.6198493743445, 585.6770306155281, 645.280363562752, 709.2978921717759, 769.4292917866055];
CornophoneData = [56.644070889128585, 115.77790823126227, 173.7432173002475, 230.0312956146474, 293.26093164942563, 348.3245541664023, 406.2196563432125, 459.523941740179, 521.0781319413594, 578.9339658564448, 641.7435504854986, 705.3872881893762, 767.4055578473059, 827.7739655391257, 885.6553833216985, 942.7115726343848];

m = [];
ff = [];
[m(:,1), ff(1)] = EFP_calc(CornophoneData);

N = length(projects);
F = [];
M = {};

for i = 1:N
    fid = fopen(strcat('./Outputs/',projects{i},'/Results/Eigenvalues.txt'));
    T = textscan(fid,"%[^\r\n]",1);
    A = textscan(fid,"%d %f %f %f %*[^\r\n ]");
    fclose(fid);
    f = A{2};
    [a, F(i)] = EFP_calc(f);
    M{i} = a;
end

fig = figure('Name', 'Equivalent Fundamental Pitch deviation');
set(fig, 'Position', [280 150 1000 800])

subplot(1,2,1)
set(gca, 'Position', [0.102,0.11,0.434,0.815]);
hold on
plot(m(2:end,1)',2:size(m,1),'-xk','linewidth',1, 'MarkerSize', 6, 'DisplayName', 'Meas.')

plnms = strings(1,N);
ccc = flipud(parula(length(plnms)));
for i = 1:N
    plnms(i) = sprintf("%d",i);
    plot(M{i}(2:end)',2:length(M{i}),'-o','linewidth',1, 'MarkerSize', 6, 'DisplayName', plnms(i), 'Color', ccc(i,:))
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
X = {'Norman Meas'};
X = categorical([X, plnms],[X, plnms]);
cc = [[0 0 0]']';
hold on
for i = 1:length(X)
    bb = bar(X(i),FF(i));
    if i<=1
        bb.FaceColor = cc(i,:);
    else
        bb.FaceColor = ccc(i-1,:);
    end
end
set(gca,'YColor','k','Box','on');
set(gca,'FontName','Helvetica','FontSize',font_size,'LineWidth',1)
%xlabel("Curve index", 'FontSize', font_size, 'Interpreter', 'latex')
ylabel("Reference pitch deviation (cents)", 'FontSize', font_size, 'Interpreter', 'none')
ax = gca;
ax.XAxis.FontSize = font_size-6;
ax.XAxis.TickLabelRotation = 45;


