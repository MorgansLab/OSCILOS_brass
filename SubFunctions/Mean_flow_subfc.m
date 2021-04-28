%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MEAN FLOW
%
% This subroutine solves the mean flow conservation equations and the mean
% flow variables are obtained throughout the geometry. The mean velocity
% and mean temperature across the geometry are then plotted and saved to
% the .pdf and .fig format.
% 
% Last update : 10/09/2020
%
% Gamma is constant, no heat addition (mean or fluctuations)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MF] = Mean_flow_subfc(CFG, GEOM)

    %% Unpack Geometry and Mean Flow parameters
    
    x_sample = GEOM.x_sample;
    r_sample = GEOM.r_sample;
    SectionIndex = GEOM.SectionIndex;
    TubeIndex = GEOM.TubeIndex;
    
    P1          = CFG.MF.P1;
    T1          = CFG.MF.T1;
    index_M1_u1 = CFG.MF.index_M1_u1;
    M1_u1       = CFG.MF.M1_u1;
    index_gamma = CFG.MF.index_gamma;

    if CFG.OUT.CL_OUT
        fprintf("\tCalculating mean flow parameters...  ");
    end
    if CFG.OUT.LOG_OUT
        write_log("\tCalculating mean flow parameters",CFG.OUT.log_filename);
    end

    %% Constant physical parameters

    Ru       = 8.3145;
    W_air    = 28.96512;
    R_air    = Ru/W_air.*1000;

    %% Calculation of the mean flow - Inlet conditions

    % The inlet gas is considered as air, since the gamma and W of the mixture
    % is nearly the same as those of air, with a error of 5%.
    % The burned gases are processed by the program to calculate the mass
    % fraction of the mixture

    gamma_cst = 1.4;
    Cp_cst = gamma_cst ./ (gamma_cst - 1) .* R_air;

    numSection = length(x_sample)-1; % number of sections

    % For T_mean(b) b means the index number of the section

    p_mean(1) = P1;
    T_mean(1) = T1;

    if index_gamma==1
        gamma(1)      = gamma_cst;
        Cp(1)         = Cp_cst;
        c_mean(1)     = (gamma(1)*R_air*T_mean(1)).^0.5;
    end

    rho_mean(1)           = p_mean(1)./(R_air.*T_mean(1)); % mean density

    %% Zero-Mach warning - AMacL 03/08/2020

    M_min = 1e-5; % If too low may incur numerical error due to ill-conditioned systems of equations

    if index_M1_u1==1 % M1 selected
        if abs(M1_u1) < M_min
            warning(newline+"Near-zero inlet Mach number encountered ("+sprintf("%.1e",M1_u1)+") - M1 will be set to " + sprintf("%.1e",M_min));
            M_mean(1) = M_min;                % mean Mach number
        else
            M_mean(1) = M1_u1;                % mean Mach number
        end
        u_mean(1)     = M_mean(1).*c_mean(1); % mean velocity

    elseif index_M1_u1==2 % u1 selected
        if abs(M1_u1) < M_min
            warning(newline+"Near-zero inlet velocity encountered ("+sprintf("%.1e",M1_u1)+") - u1 will be set to " + sprintf("%.1e",M_min * c_mean(1)));
            u_mean(1) = M_min * c_mean(1);    % mean velocity
        else
            u_mean(1) = M1_u1;                % mean velocity
        end
        M_mean(1)     = u_mean(1)./c_mean(1); % mean Mach number

    else
        error("Value of config parameter choice_M1_u1 ("+num2str(index_M1_u1)+") unrecognised - supported values are 1 and 2");
    end

    %% Calculation of the mean flow - Different sections

    for ss = 1:numSection-1 
        Theta(ss)  = (r_sample(ss+1)./r_sample(ss)).^2;       % Surface area ratio S2/S1
    end

    indexHA_num = 0;      % set the initial value to zero and it will be increased by 1 after every HA interface
    for ss = 1:numSection-1 
        % At every interface, the changes are split into two steps:
        % 1. Cross sectional surface area change
        % 2. Heat addition [not currently implemented]
        %
        % --------------step 1-------------------------
        %
        [p_mean(ss+1),rho_mean(ss+1),u_mean(ss+1)] = Fcn_calculation_TP_mean_WO_HeatAddition(p_mean(ss),rho_mean(ss),u_mean(ss),Theta(ss),gamma(ss));
        % ----------
        gamma(ss+1)   = gamma(ss);
        Cp(ss+1)      = Cp(ss);
        T_mean(ss+1)  = gamma(ss+1)/(gamma(ss+1)-1)*p_mean(ss+1)./(Cp(ss+1).*rho_mean(ss+1));                     
        c_mean(ss+1)  = ((gamma(ss+1) - 1).*Cp(ss+1).*T_mean(ss+1)).^0.5;
        M_mean(ss+1)  = u_mean(ss+1)./c_mean(ss+1);

    end

    %% Plot for the mean flow velocity and mean temperature

    if CFG.PLT.MF_PLOT

        fig2=figure('Name','Mean Flow Properties');
        if CFG.PLT.SMALL_PLOTS
            set(fig2, 'Position', [110 60 800 500])
        else
            set(fig2, 'Position', [150 75 1200 800])
        end

        ha = tight_subplot(2,1,[0.15 0],[.10 .05],[.10 .05]);

        N = length(SectionIndex);

        x_plots(1,1:N-1) = x_sample(1:N-1);
        x_plots(2,1:N-1) = x_sample(2:N);

        ColorUDF = 'b';  % color of the line

        % Figure 1 : Mean flow velocity
        axes(ha(1)); 
        hold on
        for ss = 1:N-1
            y_plots(1:2,ss) = u_mean(1,ss); % Horizontal lines
            plot([x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],'color',ColorUDF,'linewidth',2,'linestyle','-');
            if ss < N-1
                if TubeIndex(ss)~=0 || TubeIndex(ss+1)~=0 % If interpolated tube, 'staircase' velocities
                    y_plots(3,ss) = u_mean(ss+1); % Vertical lines
                    plot([x_plots(2,ss),x_plots(2,ss)],[y_plots(2,ss),y_plots(3,ss)],'color',ColorUDF,'linewidth',2,'linestyle','-');
                end
            end
        end

        set(gca,'xlim',[x_sample(1), x_sample(end)],'YAxisLocation','left','Color','w');

        % AMacL 07/2020
        yvalue_max  = max(y_plots(1,:));
        yvalue_min  = min(y_plots(1,:));
        ndps = 4; % Max. number of decimal places
        ymax = yvalue_max + 0.1*(yvalue_max-yvalue_min);
        ymin = yvalue_min - 0.1*(yvalue_max-yvalue_min);
        if ymax-ymin < 0.1^ndps
            % Increments graph at min. resolution
            ymin = round(yvalue_min - 2*0.1^ndps,ndps);
            ymax = round(yvalue_max + 2*0.1^ndps,ndps);
            ytickformat(sprintf('%%0.%df',ndps))
        end

        set(gca,'ylim',[ymin ymax])
        set(gca,'Position', [0.13 0.6 0.82 0.35]) % Leaves more room for ylabels

        xlabel('$x$ [m]','Color','k','Interpreter','LaTex');
        xticklabels('auto')

        ylabel('$\overline{u}$ [m/s]','Color','k','Interpreter','LaTex');
        yticklabels('auto')

        title('Mean flow velocity','Interpreter','latex');

        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica');%,'LineWidth',2)
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end

        % Figure 2 : Mean temperature
        axes(ha(2));
        hold on
        for ss = 1:N-1
            y_plots(1:2,ss) = T_mean(1,ss); % Horizontal lines
            plot([x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],'color',ColorUDF,'linewidth',2,'linestyle','-');
            if ss < N-1
                if TubeIndex(ss)~=0 || TubeIndex(ss+1)~=0 % If interpolated tube, 'staircase' velocities
                    y_plots(3,ss) = T_mean(ss+1); % Vertical lines
                    plot([x_plots(2,ss),x_plots(2,ss)],[y_plots(2,ss),y_plots(3,ss)],'color',ColorUDF,'linewidth',2,'linestyle','-');
                end
            end
        end 

        set(gca,'xlim',[x_sample(1), x_sample(end)],'YAxisLocation','left','Color','w'); 

        % AMacL 07/2020
        yvalue_max  = max(y_plots(1,:));
        yvalue_min  = min(y_plots(1,:));
        ndps = 2; % Max. number of decimal places
        ymax = yvalue_max + 0.1*(yvalue_max-yvalue_min);
        ymin = yvalue_min - 0.1*(yvalue_max-yvalue_min);
        if ymax-ymin < 0.1^ndps
            % Increments graph at min. resolution
            ymin = round(yvalue_min - 2*0.1^ndps,ndps);
            ymax = round(yvalue_max + 2*0.1^ndps,ndps);
            ytickformat(sprintf('%%0.%df',ndps))
        end

        set(gca,'ylim',[ymin ymax])
        set(gca,'Position', [0.13 0.1 0.82 0.35]) % Leaves more room for yticks

        xlabel('$x$ [m]','Color','k','Interpreter','LaTex');
        xticklabels('auto')

        ylabel('$\overline{T}$ [K]','Color','k','Interpreter','Latex');
        yticklabels('auto')

        title('Mean Temperature','Interpreter','latex');

        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica');%,'LineWidth',2);
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end

        % Saving the figure
        flpth = fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Initialization','Mean_flow');
        if CFG.PLT.SAVE_FIGS
            saveas(fig2,flpth,'fig')
        end
        if CFG.PLT.SAVE_PDFS
            save2pdf(flpth,fig2,300)
        end
        drawnow
    end

    if CFG.OUT.CL_OUT; fprintf("done\n"); end
    
    MF.c_mean   = c_mean;
    MF.M_mean   = M_mean;
    MF.Cp       = Cp;
    MF.p_mean   = p_mean;
    MF.T_mean   = T_mean;
    MF.u_mean   = u_mean;
    MF.rho_mean = rho_mean;
    MF.gamma    = gamma;
    MF.Theta    = Theta;
end


%% Functions used to compute the mean flow parameters throughout in the domain

function [p_mean2,rho_mean2,u_mean2] =Fcn_calculation_TP_mean_WO_HeatAddition(p1,rho1,u1,Theta,gamma)
% This function is used to calculate the mean thermoproperties changing
% across the interface between two sections with different section surface
% area.

    options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
    x0      = [1,1,1];  
    F       = @(x)myfun_p_rho_u(x, p1,rho1,u1,Theta,gamma);
    x2      = fsolve(F, x0, options);   % solve equation

    iternum     = 0;
    exitflag    = 0;
    while exitflag > 4 || exitflag < 1 && iternum <20
        [x2, ~, exitflag]   = fsolve(F, x0, options);   % solve equation
        iternum = iternum + 1;
        x0      = x2;
    end
    if iternum>19
        h = msgbox('The mean flow calculation is not converged!'); 
    end

    p_mean2     = x2(1)*p1;
    rho_mean2   = x2(2)*rho1;
    u_mean2     = x2(3)*u1;

end

function F = myfun_p_rho_u(x, p1,rho1,u1,Theta,gamma)
    p2      = x(1)*p1;
    rho2    = x(2)*rho1;
    u2      = x(3)*u1;

    F(1) = Theta*rho2*u2 - rho1*u1;
    if Theta>=1              % Area increasing
        F(2) = p2+rho2*u2^2 - (p1 + rho1*u1^2/Theta);
    elseif Theta<1          % Area decreasing
        F(2) = p2*rho1^gamma - p1*rho2^gamma;
    end
    F(3) = gamma/(gamma-1)*(p2/rho2 - p1/rho1) + 0.5*(u2^2 - u1^2);
end