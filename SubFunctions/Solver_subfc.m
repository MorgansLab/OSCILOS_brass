%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FREQUENCY-DOMAIN SOLVER
%
% This subroutine finds the eigenvalues of the system, and the
% corresponding eigenmodes. The eigenvalues are saved as an output file
% called 'Eigenvalues.txt'. A map of the eigenvalues in the range of
% interest is then plotted and saved to the .pdf and .fig format, and
% Equivalent Fundamental Pitch is calculated and plotted if there are more
% than 2 modes found. The mode shapes for the first plot_modes eigenmodes
% (or all of the modes found, whichever is smaller), represented as the
% modulus of the acoustic pressure and velocity, are also plotted and saved
% to the .pdf and .fig format.
% 
% Last update : 10/09/2020 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eigval_unique] = Solver_subfc(CFG, GEOM, MF, BCA)

    %% Unpack Geometry, Mean Flow, BC and Scan Range parameters
    
    x_sample     = GEOM.x_sample;
    r_sample     = GEOM.r_sample;
    SectionIndex = GEOM.SectionIndex;

    c_mean   = MF.c_mean;
    M_mean   = MF.M_mean;
    Cp       = MF.Cp;
    p_mean   = MF.p_mean;
    T_mean   = MF.T_mean;
    u_mean   = MF.u_mean;
    rho_mean = MF.rho_mean;
    gamma    = MF.gamma;
    Theta    = MF.Theta;
    
    R1h = BCA.R1h;
    R2h = BCA.R2h;
    
    FreqMax = CFG.SR.FreqMax;
    FreqMin = CFG.SR.FreqMin;
    FreqNum = CFG.SR.FreqNum;
    GRMax   = CFG.SR.GRMax;
    GRMin   = CFG.SR.GRMin;
    GRNum   = CFG.SR.GRNum;
    
    %% Pre-processing

    C1 =  [1 1 0;1 -1 0;1 1 0];
    C2 =  [1 1 0;1 -1 0;1 1 -1];
    
    numSection = length(x_sample)-1; % number of sections

    for ss = 1:numSection-1 
        cRatio  = c_mean(ss+1)./c_mean(ss);
        M1      = M_mean(ss);
        M2      = M_mean(ss+1);
        gamma1  = gamma(ss);
        gamma2  = gamma(ss+1);
        Thetaf   = Theta(ss);
        rho1    = rho_mean(ss);
        rho2    = rho_mean(ss+1);

        B1_temp(1,1) =   0;
        B1_temp(1,2) =   cRatio;
        B1_temp(1,3) =   M1.*cRatio;
        B1_temp(3,1) =   M1.*gamma1./(gamma1-1);
        B1_temp(3,2) =   M1.^2;
        B1_temp(3,3) =  -M1./(gamma1-1);

        B2_temp(1,1) =   0;
        B2_temp(1,2) =   Thetaf;
        B2_temp(1,3) =   Thetaf.*M2;
        B2_temp(3,1) =   cRatio.* Thetaf.*M2.*gamma2./(gamma2-1);
        B2_temp(3,2) =   cRatio.* Thetaf.*M2.^2;
        B2_temp(3,3) =  -cRatio.* Thetaf.*M2./(gamma2-1); 

        if Thetaf>=1 % Area increase
            B1_temp(2,1) =   Thetaf;
            B1_temp(2,2) = 2*M1;
            B1_temp(2,3) =   M1.^2;
            B2_temp(2,1) =   Thetaf;
            B2_temp(2,2) = 2*Thetaf.*M2;
            B2_temp(2,3) =   Thetaf.*M2.^2;
        elseif Thetaf<1 % Area decrease
            B1_temp(2,1) = 1/rho1.^gamma1;
            B1_temp(2,2) = 0;
            B1_temp(2,3) =-1/rho1.^gamma1;
            B2_temp(2,1) = 1/rho2.^gamma2;
            B2_temp(2,2) = 0;
            B2_temp(2,3) =-1/rho2.^gamma2;
        end  

    %     B1{1,ss}         = B1_temp;          % first step
    %     B2{1,ss}         = B2_temp;
    %     B1{2,ss}         = eye(3);           % second step
    %     B2{2,ss}         = eye(3);
        BC{ss}           = (B2_temp*C2)\(B1_temp*C1);
    end

    x_diff = diff(x_sample);
    tau_plus      = x_diff./(c_mean(1,:)+u_mean(1,:));
    tau_minus     = x_diff./(c_mean(1,:)-u_mean(1,:));
    tau_c         = x_diff./ u_mean(1,:);       % convection time delay

    %% Helmholtz number warning

    [a_max,a_max_i] = max(r_sample);
    if a_max_i==length(r_sample); a_max_i = a_max_i-1; end
    FreqCutoff = c_mean(a_max_i)/(2*pi*a_max)*1.8412;

    if FreqMax >= FreqCutoff
        warning(newline+"Radius encountered (%f m) at which higher-order transverse modes propagate within the frequency scan range"...
            + newline + "Cutoff frequency for this radius is %f Hz"...
            + newline + "'Low frequency' and 'narrow bore' assumptions break down under these conditions",a_max,FreqCutoff);
    end

    %% Eigenvalues calculation

    if CFG.OUT.CL_OUT
        fprintf("\tCalculating eigenvalues... ");
    end
    if CFG.OUT.LOG_OUT
        write_log("\tCalculating eigenvalues", CFG.OUT.log_filename);
    end

    % Parameters (AMacL 12/08/2020 - previously hard-coded)
    eig_prec = 10;     % Precision to which eigenvalues are evaluated [Re 1/s, Imag rad/s]
    eig_min_diff = 10; % Minimum spacing (absval) between consecutive eigenvalues before they are considered identical
    eig_min_freq = 1;  % Minimum eigenvalue frequency considered [rad/s]

    % Split the scan domain

    Scan_GRSp        = linspace( GRMin, GRMax, GRNum);
    Scan_FreqSp      = linspace( FreqMin, FreqMax, FreqNum);

    % Split the contour domain

    Cont_GRSp        = linspace( GRMin, GRMax, 10*GRNum);
    Cont_FreqSp      = linspace( FreqMin, FreqMax, 10*FreqNum);

    % Linear solver

    % AMacL 01/09/2020 Added waitbar
    if CFG.OUT.WAIT_BARS
        wx = 0;
        wb = waitbar(wx,sprintf("%d%%",wx),'Name','Calculating modes');
        nss = length(Scan_GRSp);
        nkk = length(Scan_FreqSp);
        tnum = nss*nkk;
        nxi = 100; % Number of 'increments' on waitbar
        nxx = 0;
    end

    eigen_num=1;
    eigenvalue=[];
    eigenvalue_prec=[];
    options = optimset('Display','off'); % the calculation results by fsolve will not be shown in the workspace
    for ss = 1:length(Scan_GRSp)
        GR = Scan_GRSp(ss); % Current growth rate
        for kk = 1:length(Scan_FreqSp)
            omega = 2*pi*Scan_FreqSp(kk); % Current frequency
            if GR == 0
                GR = 1;
            end
            if omega == 0
                omega =10;
            end
            s0 = GR+1i*omega;                                              % initial value
            [x,fval,exitflag] = fsolve(@(sf)Fcn_DetEqn_Linear(sf,R1h,R2h,M_mean,gamma,numSection,tau_plus,tau_minus,tau_c,BC),s0,options);   % linear                          
            if exitflag==1
                eigenvalue(eigen_num) = x;
                eigenvalue_prec(eigen_num) = floor(eigenvalue(eigen_num)./eig_prec).*eig_prec;      % this is used to set the precision
                eigen_num = eigen_num + 1;                
            end
            if CFG.OUT.WAIT_BARS
                wx = (nkk*(ss-1) + kk - 1)/(tnum-1);
                if (wx >= nxx/nxi)
                    waitbar(wx,wb,sprintf("%d%%",floor(wx*100)));
                    nxx = nxx + 1;
                end
            end
        end
    end
    if CFG.OUT.WAIT_BARS; close(wb); end


    % AMacL 12/08/2020 re-expressed these 3 operations to use logical indexing (previously loops)
    % Substituted EigValCol{1} for eigval_unique

    % Getting rid of duplicate values
    [eigval_unique_prec,m,n] = unique(eigenvalue_prec); 
    eigval_unique = eigenvalue(m); % these are the eigenvalues
    eigval_unique = sort(eigval_unique);

    % Getting rid of near-zero values
    eigval_unique = eigval_unique( abs(imag(eigval_unique)) >= eig_min_freq );

    % Getting rid of redundant values
    eigenvalue_unique_diff = diff(eigval_unique);
    eigval_unique = [eigval_unique( abs(eigenvalue_unique_diff) >= eig_min_diff ), eigval_unique(end)]; % this leaves the eigenvalues we want

    % Getting rid of values out of the range of interest
    eigval_unique = eigval_unique( ~( real(eigval_unique)<GRMin | real(eigval_unique)>GRMax...
                                      | abs(imag(eigval_unique))<(FreqMin*2*pi) | abs(imag(eigval_unique))>(FreqMax*2*pi) ) );

    if isempty(eigval_unique)
        error("Solver error: No meaningful eigenvalues were found."...
            + newline + "This can happen if:"...
            + newline + "No modes lie between min_freq and max_freq (try increasing max_freq)"...
            + newline + "No modes lie between min_GR and max_GR (try widening GR range)"...
            + newline + "The phase of a boundary reflectance is too strong a function of frequency (check BCs)");
    end

    if length(eigval_unique) < 4
        EFPref = abs(imag(eigval_unique(end))./2./pi)./length(eigval_unique);
        warning(newline+"Only %d modes found\nEFP deviation will be calculated with reference frequency %f Hz corresponding to highest mode", length(eigval_unique), EFPref);
    else
        EFPref = abs(imag(eigval_unique(4))./2./pi)./4;
    end

    EFP = EFP_calc(imag(eigval_unique)/(2*pi),EFPref); % Equivalent fundamental pitch values

    %% Store the data in the file 'Eigenvalues.txt'

    data = [1:length(eigval_unique); abs(imag(eigval_unique)./2./pi); real(eigval_unique); EFP; Which_note(abs(imag(eigval_unique)./2./pi))];

    fid=fopen(CFG.OUT.eig_filename,'w');
    fprintf(fid, 'Mode number \t Frequency [Hz] \t Growth rate [1/s] \t EFP [cents] \t Note\n');
    fprintf(fid, '%1.0f \t \t %.2f  \t \t % .2f  \t \t %+.2f  \t %s\n', data);
    fclose(fid);

    %% Map of eigenvalues

    if CFG.PLT.EIGMAP_PLOT
        if CFG.OUT.CL_OUT
            fprintf("done\n\tCalculating contour plot values... ");
        end
        if CFG.OUT.LOG_OUT
            write_log("\tCalculating contour plot values", CFG.OUT.log_filename);
        end

        % AMacL 01/09/2020 Added waitbar
        if CFG.OUT.WAIT_BARS
            wx = 0;
            wb = waitbar(wx,sprintf("%d%%",wx),'Name','Calculating contour plot');
            nuu = length(Cont_GRSp);
            nkk = length(Cont_FreqSp);
            tnum = nuu*nkk;
            nxi = 100; % Number of 'increments' on waitbar
            nxx = 0;
        end

        % Contour values calculation
        F = zeros(length(Cont_FreqSp),length(Cont_GRSp)); % AMacL 11/08/2020 defined F as transpose of previous
        for uu = 1:length(Cont_GRSp)
            GR = Cont_GRSp(uu);
            for kk  = 1:length(Cont_FreqSp)
                omega       = 2*pi*Cont_FreqSp(kk);
                s           = GR+1i*omega;
                F(kk,uu) = Fcn_DetEqn_Linear(s,R1h,R2h,M_mean,gamma,numSection,tau_plus,tau_minus,tau_c,BC); % AMacL 11/08/2020 swapped kk and uu
            end
            if CFG.OUT.WAIT_BARS
                wx = (nkk*(uu-1) + kk - 1)/(tnum-1);
                if (wx >= nxx/nxi)
                    waitbar(wx,wb,sprintf("%d%%",floor(wx*100)));
                    nxx = nxx + 1;
                end
            end
        end
        if CFG.OUT.WAIT_BARS; close(wb); end

        fig4=figure('Name','Eigenvalue Map');
        if CFG.PLT.SMALL_PLOTS
            set(fig4, 'Position', [130 80 600 600])
        else
            set(fig4, 'Position', [250 125 800 800])
        end

        hold on
        Cont_F = 20*log10(abs(F));
        Cont_F(Cont_F==-Inf) = realmin; % In the unlikely (but possible, has happened) event that F has a zero element, remove -Inf dB so saveas function doesn't break
        contourf(Cont_GRSp./100,Cont_FreqSp,Cont_F); % B.B. 05/07/2019 Suppressed output % AMacL 11/08/2020 substituted ValCol{1}' for F
        plot(real(eigval_unique)./100,abs(imag(eigval_unique))./2./pi,'p','markersize',8,'color','k','markerfacecolor',[1,1,1],'LineWidth',0.5);
        hold off
        % IMPORTANT NOTE : THERE IS NO ABS(IMAG(...)) IN THE PLOT ABOVE. THIS MAY
        % LEAD CERTAIN EIGENVALUES TO BE ABSENT OF THE CONTOUR PLOT
        % AMacL 15/07/2020 have added this in but eigenvalues are also missed if
        % scan numbers are too low

        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica');%,'LineWidth',2)
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
        xlabel(gca,  '$ \textrm{Growth rate}/100~~$ [$1/s$] ','Color','k','Interpreter','LaTex');
        ylabel(gca,'$ \textrm{Frequency}~~$ [Hz]','Color','k','Interpreter','LaTex');
        set(gca,'ylim',[FreqMin FreqMax],'YAxisLocation','left','Color','w');
        set(gca,'xlim',[GRMin   GRMax]./100);
        set(gca,'layer','top');
        colorbar 
        colormap(hot);
        hcb=colorbar;
        set(hcb,'box','on','Unit','points','location','eastoutside')

        hTitle = title('Eigenvalues are located at minima');
        set(hTitle, 'interpreter','latex', 'fontunits','points')
        if ~CFG.PLT.SMALL_PLOTS; set(hTitle, 'fontsize', 20); end
        set(hcb,'box','on','Unit','points','location','eastoutside') % Workaround to reinitialize the colorbar

        % Saving the figure
        flpth = fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Results','Eigenvalues_map');
        if CFG.PLT.SAVE_FIGS
            saveas(fig4,flpth,'fig')
        end
        if CFG.PLT.SAVE_PDFS
            save2pdf(flpth,fig4,300)
        end
        drawnow
    end

    if CFG.OUT.CL_OUT; fprintf("done\n"); end

    %% Plot Equivalent Fundamental Pitch deviation graph

    if length(eigval_unique) < 3
        warning(newline+"Only %d modes found, EFP will not be plotted",length(eigval_unique));
        CFG.PLT.EFP_PLOT = false;
    end

    if CFG.PLT.EFP_PLOT
        fig6 = figure('Name', 'Equivalent Fundamental Pitch deviation');
        if CFG.PLT.SMALL_PLOTS
            set(fig6, 'Position', [140 90 400 480])
        else
            set(fig6, 'Position', [280 150 560 700])
        end
        hold on
        plot(EFP(2:end),2:length(EFP),'-xk','linewidth',2, 'MarkerSize', 14)
        set(gca,'YColor','k','Box','on');
        set(gca,'FontName','Helvetica');%,'LineWidth',2)
        if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
        ylim([2,length(EFP)])

        xlabel("EFP (cents)", 'Interpreter', 'latex')
        ylabel("Mode number", 'Interpreter', 'latex')
        grid on
        set(gca, 'GridLineStyle', '--')

        xval_max  = max(EFP(2:end));
        xval_min  = min(EFP(2:end));
        ndps = 0; % Max. number of decimal places
        xmax = xval_max + 0.1*(xval_max-xval_min);
        xmin = xval_min - 0.1*(xval_max-xval_min);
        if xmax-xmin < 0.1^ndps
            % Increments graph at min. resolution
            xmin = round(xval_min - 2*0.1^ndps,ndps);
            xmax = round(xval_max + 2*0.1^ndps,ndps);
            xtickformat(sprintf('%%0.%df',ndps))
        end
        xlim([xmin xmax]);
        xlim([-1,1].*max(abs(xlim)));

        if max(diff(yticks))<1
            yticks(2:length(EFP));
        end

        % Saving the figure
        drawnow
        flpth=fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Results','EFP_plot');
        if CFG.PLT.SAVE_FIGS
            saveas(fig6,flpth,'fig')
        end
        if CFG.PLT.SAVE_PDFS
            save2pdf(flpth,fig6,300)
        end
        drawnow
    end

    %% Figures corresponding to the different eigenmodes

    if CFG.PLT.PLOT_MODES>0
        if CFG.OUT.LOG_OUT
            write_log("\tPlotting mode shapes:",CFG.OUT.log_filename);
        end
        n_plot_mod=min(CFG.PLT.PLOT_MODES,length(eigval_unique)); % A maximum of 10 modes are represented
        for numMode=1:n_plot_mod
            if CFG.OUT.LOG_OUT
                write_log("\t\tmode "+num2str(numMode),CFG.OUT.log_filename);
            end
            s_star = eigval_unique(numMode);             % eigenvalue

            [x_resample,p,u] = Fcn_calculation_eigenmode_Linear(s_star,R1h,tau_plus,tau_minus,tau_c,BC,c_mean,u_mean,rho_mean,numSection,x_sample);

            fig5(numMode)=figure('Name',"Mode Shape "+num2str(numMode)+" ("+num2str(sprintf('%.2f',abs(imag(eigval_unique(numMode))./2./pi)))+" Hz)",'NumberTitle','off');
            if CFG.PLT.SMALL_PLOTS
                set(fig5(numMode), 'Position', [140+10*numMode 100 800 400])
            else
                set(fig5(numMode), 'Position', [280+20*numMode 150 1200 800])
            end
            ha = tight_subplot(2,1,[0.05 0],[.10 .10],[.10 .05]);

            axes(ha(2)); 

            hold on

            for lll=1:length(x_sample)-1
                plot(x_resample(lll,:),abs(p(lll,:)),'-','color','k','Linewidth',2)
            end

            valueLevel = round(log10(max(max(abs(p)))));
            ymax1=ceil(max(max(abs(p)))./10^valueLevel).*10^valueLevel;
            ymin1=floor(min(min(abs(p)))./10^valueLevel).*10^valueLevel;
            ylimitUD=[ymin1 ymax1+0.25*(ymax1-ymin1)];
            ytickUD=linspace(ylimitUD(1),ylimitUD(2),6);

            for zz=1:length(ytickUD)
                yticklabelUD{zz}=num2str(ytickUD(zz));
            end

            yticklabelUD{end}='';
            xmax1=max(max(x_resample));
            xmin1=min(min(x_resample));
            xlimitUD=[xmin1 xmax1];
            xtickUD=linspace(xlimitUD(1),xlimitUD(2),6);

            set(gca,'YColor','k','Box','on');
            set(gca,'FontName','Helvetica');%,'LineWidth',2)
            if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
            xlabel('$x $ [m]','Color','k','Interpreter','LaTex');
            ylabel('$|~\hat{p}~|$ [Pa] ','Color','k','Interpreter','LaTex');
            set(gca,'ylim', ylimitUD,'yTick',ytickUD,'yticklabel',yticklabelUD,'YAxisLocation','left');
            set(gca,'xlim',xlimitUD,'xtick',xtickUD);
            ylimit=get(gca,'ylim');
            yticklabels('auto');
            xticklabels('auto');
            set(gca,'layer','top');

            NSp = length(x_sample);
            if NSp < 10 %&& Outlet_type ~= 8 % B.B. 07/07/2019 - added test to avoid section lines when using HX (reserved to indicate HX extent)
                for k=2:length(x_sample)-1
                   plot([x_sample(k),x_sample(k)],ylimit,'--','linewidth',0.5,'color','k') 
                end
            % NOTE : k=1:length(x_sample) was changed to k=2:length(x_sample)-1 to remove the dotted lines placed on the y-axis     
            end % B.B. 05/07/2019 STOP
            grid on
            hold off

            axes(ha(1)); 

            hold on
            for k=1:length(x_sample)-1
                plot(gca,x_resample(k,:),abs(u(k,:)),'-','color','k','Linewidth',2)
            end
            set(gca,'YColor','k','Box','on');
            set(gca,'FontName','Helvetica');%,'LineWidth',2)
            if ~CFG.PLT.SMALL_PLOTS; set(gca,'FontSize',16); end
            set(gca,'xlim',xlimitUD,'xtick',xtickUD,'xticklabel',[]);
            xlabel(gca,'','Color','k','Interpreter','LaTex');
            ylabel(gca,'$|~\hat{u}~|$ [m/s] ','Color','k','Interpreter','LaTex');
            ylimit=get(gca,'ylim');
            yticklabels('auto');
            set(gca,'layer','top');

            if NSp < 10 && CFG.BC.Outlet_type ~= 8 % B.B. 07/07/2019 - added test to avoid section lines when using HX (reserved to indicate HX extent)
                for k=2:length(x_sample)-1
                    plot([x_sample(k),x_sample(k)],ylimit,'--','linewidth',1,'color','k') 
                end
            % NOTE : k=1:length(x_sample) was changed to k=2:length(x_sample)-1 to remove the dotted lines placed on the y-axis        
            end   % B.B. 05/07/2019 STOP
            grid on
            hold off

            str = ['Mode ',num2str(numMode),' - Frequency = ',num2str(sprintf('%.2f',abs(imag(eigval_unique(numMode))./2./pi))), ' Hz'];
            hTitle_bis = title(str);
            set(hTitle_bis, 'interpreter','latex', 'fontunits','points')
            if ~CFG.PLT.SMALL_PLOTS; set(hTitle,'fontsize',20); end

            % Saving the figure
            path_solver=fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Results',strcat('Mode_',num2str(numMode)));
            if CFG.PLT.SAVE_FIGS
                saveas(fig5(numMode),path_solver,'fig')
            end
            if CFG.PLT.SAVE_PDFS
                save2pdf(path_solver,fig5(numMode),300)
            end
            drawnow
        end
    end
    
    %% FINAL OUTPUTS

    time = toc;
    text = sprintf("Calculation completed in %.2f seconds\n\tEigenvalues saved to %s",time,CFG.OUT.eig_filename);
    if CFG.OUT.FIN_MSG
        msgbox(text,'OSCILOS_brass','help');
    end
    if CFG.OUT.CL_OUT
        fprintf("%s\n",text)
    end
    if CFG.OUT.LOG_OUT
        write_log(sprintf("Calculation completed in %.2f seconds",time),CFG.OUT.log_filename);
        write_log(sprintf("\tEigenvalues saved to %s",strrep(CFG.OUT.eig_filename,'\','\\')),CFG.OUT.log_filename);
    end
    
end

%% Solver Functions

function F = Fcn_DetEqn_Linear(s,R1h,R2h,M_mean,gamma,numSection,tau_plus,tau_minus,tau_c,BC)

    % Boundary conditions - Pressure reflection coefficients
    R1_temp = R1h(s);
    R2_temp = R2h(s);
    % Bertie's entropy wave reflection coefficient
    Rs = -0.5*M_mean(end)./(1 + 0.5*(gamma(end) - 1 ).*M_mean(end)); % B.B. 10/07/2019
    
    Te = 0; % No indirect noise

    G = eye(3);
    for ll = 1:numSection-1
        if Te==0 % AMacL 07/2020 - remove -> Inf contributions from exp(-s*tau_c) for low mean flow
            D1 = diag([ exp(-s*tau_plus(ll)), exp( s*tau_minus(ll)), 0]);
        else
            D1 = diag([ exp(-s*tau_plus(ll)), exp( s*tau_minus(ll)), exp(-s*tau_c(ll))]);
        end
        Z = BC{ll}*D1; % BC{k} is (B_k2*C_k2)^-1 * (B_k1*C_k1)
        G = Z*G;
    end

    A1_minus        = 1000;
    A1_plus         = R1_temp.*A1_minus;
    E1 = 0;
    Array_LeftBD    = [A1_plus, A1_minus, E1].';
    
    if Te==0 % AMacL 07/2020
        % Prevent Te.*exp(-s*tau_c) going -> 0*Inf (=NaN) limit for very large tau_c caused by low mean flow
        D1End = diag([exp(-s*tau_plus(end)), exp( s*tau_minus(end)), 0]);
    else
        D1End = diag([exp(-s*tau_plus(end)), exp( s*tau_minus(end)), Te.*exp(-s*tau_c(end))]);
    end

    Array_RightBD   = D1End*G*Array_LeftBD;
    AN_plus         = Array_RightBD(1);
    AN_minus        = Array_RightBD(2);
    EN_plus         = Array_RightBD(3);
    
    F = (R2_temp.*AN_plus + Rs.*EN_plus) - AN_minus;  

end

function [x_resample,p,u]=Fcn_calculation_eigenmode_Linear(s_star,R1h,tau_plus,tau_minus,tau_c,BC,c_mean,u_mean,rho_mean,numSection,x_sample)
    % This function is used to plot the modeshape of selected mode, 
    % Linear

    R1 = R1h(s_star);

    A_minus(1)  = 1;
    A_plus(1)   = R1.*A_minus(1);
    E(1)        = 0;
    Array(:,1)  = [A_plus(1),A_minus(1), E(1)].'; % B.B. 05/07/2019 - Suppressed output

    for pp = 1:numSection-1 
        D1 = diag([ exp(-s_star*tau_plus(pp)), exp( s_star*tau_minus(pp)), exp(-s_star*tau_c(pp))]);
        if E(1)==0 % AMacL 07/2020 - remove -> Inf contributions from exp(-s*tau_c) for low mean flow
            D1(3,3) = 0;
        end
        Z{pp} = BC{pp}*D1;
        Array(:,pp+1) = Z{pp}*Array(:,pp); % B.B. 05/07/2019 - output suppressed
    end

    for k = 1:numSection-1
        A_plus(k+1)     = Array(1,k+1);
        A_minus(k+1)    = Array(2,k+1);
        E(k+1)          = Array(3,k+1);
    end

    for k=1:length(x_sample)-1
        x_resample(k,:)=linspace(x_sample(k),x_sample(k+1),201);    % resample of x-coordinate
    end

    for k = 1:length(x_sample)-1
        kw1_plus(k)  = s_star./(c_mean(k)+u_mean(k));
        kw1_minus(k) = s_star./(c_mean(k)-u_mean(k));

        p(k,:) =    A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))+...
                    A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1)));
        u(k,:) =    (A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))-...
                    A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1))))./rho_mean(k)./c_mean(k);
    end
end
