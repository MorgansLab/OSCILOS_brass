%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BOUNDARY CONDITIONS
%
% This subroutine loads the acoustic boundary conditions at the inlet and
% outlet respectively. The acoustic reflection coefficient at the inlet and
% outlet are then plotted and saved to the .pdf and .fig format.
% 
% Last update : 10/09/2020
% 
% Boundary conditions of type
% 1 (Open end),
% 2 (Closed end),
% 4 (Given amplitude and time delay) and
% 5 (Given amplitude and phase lag)
% are retained from OSCILOS_lite cold flow code
% 
% BC type 6 (Polynomial from coefficients) reinstated from OSCILOS_long
%
% BC type 9 (Interpolated from user gain/phase data),
% 10 (User-defined gain-phase functions),
% 11&12 (Levine-Schwinger),
% 13 (Vibrating plate at pipe end) and
% 14 (Vibrating plate behind orifice) added for OSCILOS_brass
%
% Definition of reflection coefficent changed to implement anonymous
% functions and solver arguments updated accordingly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BC] = BC_subfc(CFG, GEOM, MF)

    %% Unpack BC parameters

    Inlet_type = CFG.BC.Inlet_type;
    Inlet_param1 = CFG.BC.Inlet_param1;
    Inlet_param2 = CFG.BC.Inlet_param2;
    Inlet_param3 = CFG.BC.Inlet_param3;
    Outlet_type = CFG.BC.Outlet_type;
    Outlet_param1 = CFG.BC.Outlet_param1;
    Outlet_param2 = CFG.BC.Outlet_param2;
    Outlet_param3 = CFG.BC.Outlet_param3;

    if CFG.OUT.CL_OUT
        fprintf("\tCalculating boundary conditions...  ");
    end
    if CFG.OUT.LOG_OUT
        write_log("\tCalculating boundary conditions",CFG.OUT.log_filename);
    end

    %% Frequency sampling

    f_sample    = linspace(CFG.SR.FreqMin, CFG.SR.FreqMax, CFG.SR.FreqNum);  % sampling frequency
    s           = 2*pi*f_sample*1i;                     % s = i omega

    %% Inlet setup (BCs numbered as in OSCILOS_long)

    switch Inlet_type
        case 1                                      % open end, R=-1;
            R1h = @(sf) -1 .* ones(1,length(sf));

        case 2                                      % closed end, R=1; 
            R1h = @(sf) 1 .* ones(1,length(sf));

    % 	 case 3                                     % choked (not implemented)
    %        gamma           = MF.gamma(1,1);
    %        M1              = MF.M_mean(1,1);
    %        TEMP            = gamma*M1/(1+(gamma-1)*M1^2);
    %        num1      = (1-TEMP)/(1+TEMP);
    %        R1h = @(sf) num1 .* ones(1,length(sf));

        case 4                                      % time lag
            tau_d1 = Inlet_param2./1000; % IN MILLISECONDS
            R1h = @(sf) Inlet_param1.*exp(-sf.*tau_d1);

        case 5                                      % phase lag
            Phi = Inlet_param2;
            if Inlet_param3 == 1 % Use param3 == 1 to indicate phase in degrees
                Phi = Phi*pi/180;
            end
            R1h = @(sf) Inlet_param1.*exp(-Phi.*1i).*ones(1,length(sf));

        case 6                                      % Polynomial with time delay from coefficients
            C_inlet_BCpoly = getBCpoly(CFG.config_filename, Inlet_param2, Inlet_param3, "INLET_POLYNOMIAL_COEFFICIENTS");

            num1 = cell2mat(C_inlet_BCpoly{1});
            den1 = cell2mat(C_inlet_BCpoly{2});
            tau_d1 = Inlet_param1./1000;

            R1h = @(sf) polyval(num1,sf)./polyval(den1,sf).*exp(-sf*tau_d1);

    %     case 7                                    % Polynomial with time delay from user data (not yet implemented)

    %     case 8                                    % Something to do with heat exchangers (not applicable to inlet)

        case 9                                      % Gain and phase interpolated from data
            C_inlet_BCdata = getBCdata(CFG.config_filename, "INLET_GAIN_PHASE_DATA");

            In_freq      = C_inlet_BCdata(:,1); % Frequency
            In_gain      = C_inlet_BCdata(:,2); % Gain        
            In_phase     = C_inlet_BCdata(:,3); % Phase

            if Inlet_param1 == 1 % Use param1 == 1 to indicate phase in degrees
                In_phase = In_phase.*pi./180;
            end

            interp_type = 'makima';
            R1h = @(sf) interp1(In_freq, In_gain, imag(sf)./(2*pi), interp_type).*exp(-1i.*interp1(In_freq, In_phase, imag(sf)./(2*pi), interp_type));

        case 10                                     % Arbitrary function for gain and phase
            C_inlet_funcs = getBCfuncs(CFG.config_filename, "INLET_GAIN_PHASE_FUNCTIONS");

            InletBCfunc_Gain = eval(C_inlet_funcs{1});
            InletBCfunc_Phase = eval(C_inlet_funcs{2});

            if Inlet_param1 == 1 % Use param1 == 1 to indicate phase in degrees
                Inlet_phase_mod = pi/180;
            else
                Inlet_phase_mod = 1;
            end

            R1h = @(sf) InletBCfunc_Gain(imag(sf)./(2*pi)).*exp(-1i.*InletBCfunc_Phase(imag(sf)./(2*pi)).*Inlet_phase_mod);

        case 11                                     % Levine-Schwinger by Norris-Sheng polynomial
            a_LS = GEOM.r_sample(1); % inlet radius
            R1h = @(sf) polyLevineSchwingerBC(sf, a_LS, MF.c_mean(1));
            He_max = a_LS * CFG.SR.FreqMax * 2 * pi / MF.c_mean(1);
            if He_max > 3.832
                warning(newline+"Inlet Levine-Schwinger condition, and polynomial approximation thereof, breaks down inside the frequency scan range"...
                    + newline + "Maximum viable frequency for inlet radius is "+sprintf("%.2f Hz", 3.832.*MF.c_mean(1)./a_LS./2./pi))
            end

        case 12                                     % Flanged Levine-Schwinger by Norris-Sheng polynomial
            a_LS = GEOM.r_sample(1); % inlet radius
            R1h = @(sf) polyFlangedLevineSchwingerBC(sf, a_LS, MF.c_mean(1));
            He_max = a_LS * CFG.SR.FreqMax * 2 * pi / MF.c_mean(1);
            if He_max > 3.832
                warning(newline+"Inlet flanged Levine-Schwinger condition, and polynomial approximation thereof, breaks down inside the frequency scan range"...
                    + newline + "Maximum viable frequency for inlet radius is "+sprintf("%.2f Hz", 3.832.*MF.c_mean(1)./a_LS./2./pi))
            end

        case 13                                     % Oscillating Piston at pipe end (Mawardi 1951)
            r_PB = GEOM.r_sample(1); % inlet radius
            a_PB = Inlet_param1.*r_PB; % piston radius ratio
            R1h = @(sf) pistonBC(sf, r_PB, a_PB, MF.c_mean(1));

        case 14                                     % Oscillating Masses of orifice in front of piston (Iwanov-Schitz/Rscherkin 1963)
            r_OM = GEOM.r_sample(1); % inlet radius
            a_OM = Inlet_param1.*r_OM; % param1 is piston radius ratio
            b_OM = Inlet_param2.*r_OM; % param2 is orifice radius ratio
            d_OM = Inlet_param3.*r_OM; % param3 is nondimensional orifice offset
            R1h = @(sf) orificeBC(sf, r_OM, a_OM, b_OM, d_OM, MF.c_mean(1));

    end


    %% Outlet setup (BCs numbered as in OSCILOS_long)

    switch Outlet_type
        case 1                            % open end, R=-1;
            R2h = @(sf) -1 .* ones(1,length(sf));

        case 2                            % closed end, R=1; 
            R2h = @(sf) 1 .* ones(1,length(sf));

    %     case 3                          % choked (not implemented)
    %        gamma   = MF.gamma(1,end);
    %        M2      = MF.M_mean(1,end);
    %        TEMP    = (gamma-1)*M2/2;
    %        num2  = (1-TEMP)/(1+TEMP);
    %        R2h = @(sf) num2 .* ones(1,length(sf));

        case 4                            % time lag
            tau_d2 = Outlet_param2./1000;
            R2h = @(sf) Outlet_param1.*exp(-sf.*tau_d2);

        case 5                            % phase lag
            Phi = Outlet_param2;
            if Outlet_param3 == 1 % Use param3 == 1 to indicate phase in degrees
                Phi = Phi*pi/180;
            end
            R2h = @(sf) Outlet_param1.*exp(-Phi.*1i).*ones(1,length(sf));

        case 6                            % Polynomial with time delay from coefficients
            C_outlet_BCpoly = getBCpoly(CFG.config_filename, Outlet_param2, Outlet_param3, "OUTLET_POLYNOMIAL_COEFFICIENTS");

            num2 = cell2mat(C_outlet_BCpoly{1});
            den2 = cell2mat(C_outlet_BCpoly{2});
            tau_d2 = Outlet_param1./1000;

            R2h = @(sf) polyval(num2,sf)./polyval(den2,sf).*exp(-sf*tau_d2);

    %     case 7                          % Polynomial with time delay from user data (not yet implemented)

    %     case 8                          % Something to do with heat exchangers (not implemented)

        case 9                            % Gain and phase interpolated from data
            C_outlet_BCdata = getBCdata(CFG.config_filename, "OUTLET_GAIN_PHASE_DATA");

            Out_freq      = C_outlet_BCdata(:,1); % Frequency
            Out_gain      = C_outlet_BCdata(:,2); % Gain
            Out_phase     = C_outlet_BCdata(:,3); % Phase

            if Outlet_param1 == 1 % Use param1 == 1 to indicate phase in degrees
                Out_phase = Out_phase.*pi./180;
            end

            interp_type = 'makima';
            R2h = @(sf) interp1(Out_freq, Out_gain, imag(sf)./(2*pi), interp_type).*exp(-1i.*interp1(Out_freq, Out_phase, imag(sf)./(2*pi), interp_type));

        case 10                           % Arbitrary function for gain and phase
            C_outlet_funcs = getBCfuncs(CFG.config_filename, "OUTLET_GAIN_PHASE_FUNCTIONS");

            OutletBCfunc_Gain = eval(C_outlet_funcs{1});
            OutletBCfunc_Phase = eval(C_outlet_funcs{1});

            fclose(fid);

            if Outlet_param1 == 1 % Use param1 == 1 to indicate phase in degrees
                Outlet_phase_mod = pi/180;
            else
                Outlet_phase_mod = 1;
            end

            R2h = @(sf) OutletBCfunc_Gain(imag(sf)./(2*pi)).*exp(-1i.*OutletBCfunc_Phase(imag(sf)./(2*pi)).*Outlet_phase_mod);

        case 11                           % Levine-Schwinger by Norris-Sheng polynomial
            a_LS = GEOM.r_sample(end); % outlet radius
            R2h = @(sf) polyLevineSchwingerBC(sf, a_LS, MF.c_mean(end));
            He_max = a_LS * CFG.SR.FreqMax * 2 * pi / MF.c_mean(end);
            if He_max > 3.832
                warning(newline+"Outlet Levine-Schwinger condition, and polynomial approximation thereof, breaks down inside the frequency scan range"...
                    + newline + "Maximum frequency for outlet radius is "+sprintf("%.2f Hz", 3.832.*MF.c_mean(1)./a_LS./2./pi))
            end

        case 12                           % Flanged Levine-Schwinger by Norris-Sheng polynomial
            a_LS = GEOM.r_sample(end); % outlet radius
            R2h = @(sf) polyFlangedLevineSchwingerBC(sf, a_LS, MF.c_mean(end));
            He_max = a_LS * CFG.SR.FreqMax * 2 * pi / MF.c_mean(end);
            if He_max > 3.832
                warning(newline+"Outlet flanged Levine-Schwinger condition, and polynomial approximation thereof, breaks down inside the frequency scan range"...
                    + newline + "Maximum frequency for outlet radius is "+sprintf("%.2f Hz", 3.832.*MF.c_mean(1)./a_LS./2./pi))
            end

        case 13
            r_PB = GEOM.r_sample(end); % outlet radius
            a_PB = Outlet_param1*r_PB; % piston radius
            R2h = @(sf) pistonBC(sf, r_PB, a_PB, MF.c_mean(end));

        case 14                                     % Oscillating Masses of orifice in front of piston (Iwanov-Schitz/Rscherkin 1963)
            r_OM = GEOM.r_sample(end); % inlet radius
            a_OM = Outlet_param1.*r_OM; % param 1 is piston radius ratio
            b_OM = Outlet_param2.*r_OM; % param 2 is orifice radius ratio
            d_OM = Outlet_param3.*r_OM; % param 3 is nondimensional orifice offset
            R2h = @(sf) orificeBC(sf, r_OM, a_OM, b_OM, d_OM, MF.c_mean(end));
    end

    %% Plotting the reflection coefficient

    % Setup

    R1 = R1h(s);
    R2 = R2h(s);

    if CFG.PLT.BC_PLOT
        fig3 = figure('Name','Boundary Conditions - Reflection Coefficients');
        if CFG.PLT.SMALL_PLOTS
            set(fig3, 'Position', [120 70 800 600])
        else
            set(fig3, 'Position', [200 100 1200 800])
        end

        ha = tight_subplot(2,2,[0.06 0.12],[.10 .10],[.10 .05]);

        % Modulus of the upstream reflection coefficient
        axes(ha(1)); 

        hold on
        hLine1 = plot(f_sample,abs(R1),'-','color','b','Linewidth',2);
        set(gca,'Position',[0.1,0.56,0.365,0.37]);
        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica')%,'LineWidth',1)
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
        xlabel(gca,'','Color','k','Interpreter','LaTex');
        ylabel(gca,'Gain [-]','Color','k','Interpreter','LaTex')
        %set(gca,'xlim',[0 1000],'xTick',200:200:1000,'xticklabel',{},'YAxisLocation','left','Color','w');
        %set(gca,'ylim',[0 1.2],'yTick',0:0.2:1.2,'yticklabel',{'',0.2:0.2:1})
        if (Inlet_type==9) % plot interpolation points
            pts1=scatter(In_freq, In_gain, 80, 'g', 'LineWidth', 2);
            legend("Interpolated", "Data Points", 'Location', 'best');
        end

        % AMacL 07/2020
        Ylim=get(gca,'ylim');
        ymin = Ylim(1);
        ymax = Ylim(2);
        ndps = 2; % Max. number of decimal places
        if ymax-ymin < 0.1^ndps
            % Increments graph at min. resolution
            ymin = round(ymin - 2*0.1^ndps,ndps);
            ymax = round(ymax + 2*0.1^ndps,ndps);
            ytickformat(sprintf('%%0.%df',ndps))
        end
        ylim([ymin,ymax]);
        xticklabels('auto')
        yticklabels('auto')
        %set(gca,'box','on','linewidth',2,'FontSize',16);
        title1 = title({'$\mathcal{R}_{inlet}$'},'interpreter','latex');
        if ~CFG.PLT.SMALL_PLOTS; set(title1, 'fontsize', 20); end
        hold off

        % Phase of the upstream reflection coefficient
        axes(ha(3)); 

        hold on
        hLine2=plot(f_sample,unwrap(angle(R1),1.9*pi)./pi,'-','color', 'b','Linewidth',2);
        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica');%,'LineWidth',1)
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
        xlabel(gca,'$f$ [Hz]','Color','k','Interpreter','LaTex');%,'Position',[1000,-1.35,-1]);
        ylabel(gca,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex')
        %set(gca,'xlim',get(ha(1),'xlim'),'xTick',get(ha(1),'xTick'),'YAxisLocation','left','Color','w');
        if (Inlet_type==9) % plot interpolation points
            pts2=scatter(In_freq, unwrap(In_phase./-pi,1.9*pi), 80, 'g', 'LineWidth', 2);
            legend("Interpolated", "Data Points", 'Location', 'best');
        end
        xticklabels('auto')
        yticklabels('auto')
        %set(gca,'box','on','linewidth',2,'FontSize',16);
        hold off

        % Modulus of the dowstream reflection coefficient
        axes(ha(2)); 

        hold on
        hLine1=plot(f_sample,abs(R2),'-','color','r','Linewidth',2);
        set(gca,'Position',[0.585,0.56,0.365,0.37]);
        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica');%,'LineWidth',1)
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
        xlabel(gca,'','Color','k','Interpreter','LaTex');
        ylabel(gca,'Gain [-]','Color','k','Interpreter','LaTex')
        %set(gca,'xlim',[0 1000],'xTick',200:200:1000,'xticklabel',{},'YAxisLocation','left','Color','w');
        %set(gca,'ylim',[0 1.2],'yTick',0:0.2:1.2,'yticklabel',{'',0.2:0.2:1})
        if (Outlet_type==9) % plot interpolation points
            pts3=scatter(Out_freq, Out_gain, 80, 'g', 'LineWidth', 2);
            legend("Interpolated", "Data Points", 'Location', 'best');
        end

        % AMacL 07/2020
        Ylim=get(gca,'ylim');
        ymin = Ylim(1);
        ymax = Ylim(2);
        ndps = 2; % Max. number of decimal places
        if ymax-ymin < 0.1^ndps
            % Increments graph at min. resolution
            ymin = round(ymin - 2*0.1^ndps,ndps);
            ymax = round(ymax + 2*0.1^ndps,ndps);
            ytickformat(sprintf('%%0.%df',ndps))
        end
        ylim([ymin,ymax]);
        xticklabels('auto')
        yticklabels('auto')
        %set(gca,'box','on','linewidth',2,'FontSize',16);
        title2 = title({'$\mathcal{R}_{outlet}$'},'interpreter','latex');
        if ~CFG.PLT.SMALL_PLOTS; set(title2, 'fontsize', 20); end
        hold off

        % Phase of the downstream reflection coefficient
        axes(ha(4)); 

        hold on
        hLine2=plot(f_sample,unwrap(angle(R2),1.9*pi)./pi,'-','color', 'r','Linewidth',2);
        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica');%,'LineWidth',1)
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
        xlabel(gca,'$f$ [Hz]','Color','k','Interpreter','LaTex');%, 'Position',[1000,0.531,-1]);
        ylabel(gca,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex')
        %set(gca,'xlim',get(ha(1),'xlim'),'xTick',get(ha(1),'xTick'),'YAxisLocation','left','Color','w');
        if (Outlet_type==9) % plot interpolation points
            pts2=scatter(Out_freq, unwrap(Out_phase./-pi,1.9*pi), 80, 'g', 'LineWidth', 2);
            legend("Interpolated", "Data Points", 'Location', 'best');
        end
        xticklabels('auto')
        yticklabels('auto')
        %set(gca,'box','on','linewidth',2,'FontSize',16);
        hold off

        % Saving the figure
        flpth = fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Initialization','Boundary_Conditions');
        if CFG.PLT.SAVE_FIGS
            saveas(fig3,flpth,'fig')
        end
        if CFG.PLT.SAVE_PDFS
            save2pdf(flpth,fig3,300)
        end
        drawnow
    end

    if CFG.OUT.CL_OUT; fprintf("done\n"); end
    
    BC.R1h = R1h;
    BC.R2h = R2h;
end

%% BC File IO Functions
function [D_mat] = getBCdata(fname, header_string)
    C     = {[]};
    D_mat = [];
    fid = fopen(fname);
    while ~strcmp(fgetl(fid),header_string) && ~feof(fid)
    end
    if feof(fid)
        warning(newline+"BC data could not be read from file: "+header_string+" not found");
    else
        while isempty(C{1}) && ~feof(fid)
            C = textscan(fid,'%f %f %f'); % Get data
            if isempty(C{1})
                fgets(fid); % skip line, try again
            end
        end
        if ~isempty(C{1})
            D_mat=cell2mat(C);
        else
            warning(newline+"BC data could not be read from file: data not found");
        end
    end
    fclose(fid);
end

function [C] = getBCfuncs(fname, header_string)
    fid = fopen(fname);
    while ~strcmp(fgetl(fid),header_string) && ~feof(fid)
    end
    if feof(fid)
        warning(newline+"BC functions could not be read from file: " + header_string + " not found");
    else
        C = {fgetl(fid),fgetl(fid)};
        if ~(C{1}(1)=='@') || ~(C{2}(1)=='@')
            warning(newline+"BC functions are not of expected format");
        end
    end
    fclose(fid);
end

function [C] = getBCpoly(fname, Nnum, Nden, header_string)
    C = {{[]},{[]}};
    L = "";
    pos = 0;
    fid = fopen(fname);
    while ~strcmp(fgetl(fid),header_string) && ~feof(fid)
    end
    if feof(fid)
        warning(newline+"BC polynomial coefficients could not be read from file: "+header_string+" not found");
    else
        % Numerator
        while isempty(C{1}{1}) && ~feof(fid)
            L = fgetl(fid);
            [C{1},pos] = textscan(L,'%f', Nnum+1); % Get data
        end
        if length(C{1}{1})<Nnum+1
            warning(newline+"Expected BC polynomial numerator coefficient(s) not found: order "+num2str(Nnum)+" specified so "+num2str(Nnum+1)+" coefficients expected");
        end
        if pos+1<length(L)
            D = textscan(L(pos+1:end),'%f');
            if ~isempty(D{1})
                warning(newline+num2str(length(D{1}))+" unexpected BC polynomial numerator coefficient(s) found: order "+num2str(Nnum)+" specified so "+num2str(Nnum+1)+" coefficients expected");
            end
        end
        
        %Denominator
        while isempty(C{2}{1}) && ~feof(fid)
            L = fgetl(fid);
            [C{2},pos] = textscan(L,'%f', Nden+1); % Get data
        end
        if length(C{2}{1})<Nden+1
            warning(newline+"Expected BC polynomial denominator coefficient(s) not found: order "+num2str(Nden)+" specified so "+num2str(Nden+1)+" coefficients expected");
        end
        if pos+1<length(L)
            D = textscan(L(pos+1:end),'%f');
            if ~isempty(D{1})
                warning(newline+num2str(length(D{1}))+" unexpected BC polynomial denominator coefficient(s) found: order "+num2str(Nden)+" specified so "+num2str(Nden+1)+" coefficients expected");
            end
        end
    end
    fclose(fid);
end

%% Levine-Schwinger functions

function [R] = polyLevineSchwingerBC(sf_LS, a_LS, c) % (sf = i omega values, a_LS = radius, c = soundspeed)
    % See Norris and Sheng 1989
    k_LS = imag(sf_LS) / c; % wavenumber
    He_LS = a_LS.*k_LS;
    alpha_LS = 0.2;
    beta_LS = -0.084;
    Rmod_LS = (1 + alpha_LS.*He_LS + beta_LS.*(He_LS.^2)) ./ (1 + alpha_LS.*He_LS + (0.5 + beta_LS).*(He_LS.^2));
    l_a_LS = (0.6133 + 0.027.*(He_LS.^2)) ./ (1 + 0.19.*(He_LS.^2));
    R = -1 .* Rmod_LS .* exp(-2i .* He_LS .* l_a_LS); % -ve exponent harmonic convention
end

function [R] = polyFlangedLevineSchwingerBC(sf_LS, a_LS, c) % (sf = i omega values, a_LS = radius, c = soundspeed)
    % See Norris and Sheng 1989
    k_LS = imag(sf_LS) / c; % wavenumber
    He_LS = a_LS.*k_LS;
    alpha_LS = 0.323;
    beta_LS = -0.077;
    Rmod_LS = (1 + alpha_LS.*He_LS + beta_LS.*(He_LS.^2)) ./ (1 + alpha_LS.*He_LS + (1 + beta_LS).*(He_LS.^2));
    l_a_LS = (0.82159 - 0.49.*(He_LS.^2)) ./ (1 - 0.46.*(He_LS.^3));
    R = -1 .* Rmod_LS .* exp(-2i .* He_LS .* l_a_LS); % -ve exponent harmonic convention
end

%% Piston BC functions

function [R] = pistonBC(sf, r, a, c)
    k = imag(sf) ./ c;
    Zs_Z0 = 1 + 1i.*imag(sf)./c.*sumSeries(@(n) eval_piston_sum_term(n, r, a, k)); % Plus because -ve harmonic convention and phase lag?
    ZA = Zs_Z0.*(r.^2/a.^2);
    R = (ZA - 1)./(ZA + 1);
end

function [R] = orificeBC(sf, r, a, b, d, c)
    S = pi .* b.^2;
    A = pi .* r.^2;
    k = imag(sf) ./ c;
    persistent M_sum
    if isempty(M_sum)
        M_I_rho  = 4 .* S .* r .* sumSeries(@(n) eval_orificeI_sum_term(n, r, b, d), 1e-12, 1e6, 100);
        M_II_rho = 4 .* S .* r .* sumSeries(@(n) eval_orificeII_sum_term(n, r, b), 1e-12, 1e6, 100);
        M_sum = M_I_rho + M_II_rho;
    end
    Z_m = 1i .* k .* M_sum ./ A;
    % Z_i = (1i .* tan(k.*t) + Z_m) ./ (1 + 1i.*Z_m.*tan(k.*t));
    % If piston at end is full pipe diameter then impedance there is 1
    % so add to contribution from the referred orifice impedance
    % Z_t = Z_i + 1;
    
    Z1 = 1i .* (1 - cos(k.*d)) ./ sin(k.*d);
    Z2 = -1i ./ sin(k.*d);
    if k(1)==0
        Z1(1)=0;
        Z2(1)=Inf;
    end
    
    Zs_Z0 = 1 + 1i.*imag(sf)./c.*sumSeries(@(n) eval_piston_sum_term(n, r, a, k)); % Plus because -ve harmonic convention and phase lag?
    Zpiston = Zs_Z0.*(r.^2/a.^2);

    Z_t = Z_m + Z1 + parallelComb(Z2, Z1+Zpiston);
    
    R = (Z_t - 1)./(Z_t + 1);
end

function [f] = eval_piston_sum_term(n, r, a, k)
    k_0nR = besselderivzero(0, n);
    k_0n = k_0nR ./ r;
    f = besselj(1,k_0n.*a).^2 ./ (k_0n.^2 .* sqrt(k_0n.^2 - k.^2) .* besselj(0,k_0n.*a).^2);
end

function [f_I] = eval_orificeI_sum_term(n, r, a, t)
    xmR = besselderivzero(0, n);
    xm = xmR ./ r;
    f_I = ( besselj(1, xm.*a).^2 .* coth(xm.*t)) ./ ( (xm.*r).^3 .* besselj(0,xm.*r).^2 );
end

function [f_II] = eval_orificeII_sum_term(n, r, a)
    xmR = besselderivzero(0, n);
    xm = xmR ./ r;
    f_II = ( besselj(1, xm.*a).^2 ) ./ ( (xm.*r).^3 .* besselj(0,xm.*r).^2 );
end

function [S] = sumSeries(f, t, nmax, nmin)
    if nargin<2, t = 10^-12; end % default tolerance
    if nargin<3, nmax = 100; end % default max. iterations
    if nargin<4, nmin = 2; end % default min. iterations
    S=0;
    A=t;
    for i=1:nmax
        S = S + f(i);
        if all(abs(S-A) < t) && (i >= nmin)
            return
        end
        A = S;
    end
    warning(newline+"Sum expression for piston BC impedance did not converge to specified tolerance");
end

function [Z] = parallelComb(Z1, Z2, varargin)
    Z = 1./Z1 + 1./Z2;
    for i=1:length(varargin)
        Z = Z + 1./varargin{i};
    end
    Z = 1./Z;
end
