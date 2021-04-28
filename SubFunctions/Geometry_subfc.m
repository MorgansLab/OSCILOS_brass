%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GEOMETRY
%
% This subroutine loads the geometry from the file 'Geometry.txt', plots it
% and saves it to the .pdf and .fig format.
%
% Support for the interpolation of several other functions has been added
% to the OSCILOS_lite parser to ease the specification of sections with
% analytically varying geometries, the sections for which are interpolated
% to achieve a representative 'staircase' function.
% 
% Last update : 10/09/2020
% 
% No Liners, no Helmholtz resonators, no heat source (mean or fluctuations)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GEOM] = Geometry_subfc(CFG)

    %% Retrieving the data from the input file

    if CFG.OUT.CL_OUT
        fprintf("\tParsing geometry file %s...  ",CFG.GEOM.geom_filename);
    end
    if CFG.OUT.LOG_OUT
        write_log("\tParsing geometry file "+strrep(CFG.GEOM.geom_filename,'\','\\'),CFG.OUT.log_filename);
    end

    fid=fopen(CFG.GEOM.geom_filename);

    C_title         = textscan(fid, '%[^0-9\t\r\n, ]'); % read title
    C_cell          = textscan(fid, repmat('%f ',[1,size(C_title{1},1)])); % read numeric data
    fclose(fid);

    x_sample      = C_cell{1}; % Axial references                   
    r_sample      = C_cell{2}; % Radius              
    SectionIndex  = C_cell{3}; % Section Index     
    TubeIndex     = C_cell{4}; % Tube Index

    if size(C_title{1},1) > 4
        TubeSplit   = C_cell{5}; % No. of tube subdivisions (0 will use default interpolation if TubeIndex > 0)
    else
        TubeSplit   = zeros(size(C_cell{1},1),1); % If unspecified, use default interpolation
    end

    %% Sorting geometry in order of ascending axial co-ordinate

    [x_sample, g_ord] = sort(x_sample);

    r_sample     = r_sample(g_ord);
    SectionIndex = SectionIndex(g_ord);
    TubeIndex    = TubeIndex(g_ord);
    TubeSplit    = TubeSplit(g_ord);

    %% Updating the geometry for a gradually changing cross section area

    % Warn if final TubeIndex is nonzero
    if TubeIndex(end)~=0
        warning(newline+"Geometry file: Final section should have TubeIndex 0 (currently " ...
        + num2str(TubeIndex(end)) ...
        + ") - this will be taken as zero.")
        TubeIndex(end)=0;
    end

    % Warn of unsupported section indices
    indexUS = find(SectionIndex~=0); % & SectionIndex~=10 & SectionIndex~=11); % add these in to support heat addition sections
    if ~isempty(indexUS)
        warning(newline+"Geometry file: SectionIndex value(s) " + num2str(SectionIndex(indexUS)) ...
            + " unsupported in this version and will be taken as zero.")
        SectionIndex(indexUS) = 0;
    end

    % Warn if only 2 sections provided (3 currently required)
    if length(x_sample)<2
        error(newline+"Geometry file: insufficient number of sections to define geometry"...
            + newline + "Only " + num2str(length(x_sample)) + " provided.")
    elseif (length(x_sample)==2 && TubeIndex(1)==0)
        warning(newline+"Geometry file: only 2 sections provided - solver requires 3."...
            + newline + "A 3rd section will be interpolated at the midpoint with TubeIndex matching the inlet.")
        x_sample(3) = x_sample(2);
        r_sample(3) = r_sample(2);
        SectionIndex(3) = SectionIndex(2);
        TubeIndex(3) = TubeIndex(2);

        x_sample(2) = 0.5*(x_sample(1)+x_sample(3));
        r_sample(2) = 0.5*(r_sample(1)+r_sample(3));
        SectionIndex(2) = 0;
        TubeIndex(2) = TubeIndex(1);
    end

    indexGC   = find(TubeIndex ~= 0);
    isNoGC    = isempty(indexGC);
    TubeProf  = cell(length(indexGC),5);

    if isNoGC==0 % only runs if TubeIndex ~= 0 somewhere
        % TubeSplit defines no. of tubes and NumSplit defines no. of sections bounding these tubes so is 1 greater
        NumSplit = 1 + TubeSplit(indexGC);
        NumSplit(NumSplit==1) = 50; % Default interpolation number
        for ss = 1:length(indexGC)
        % Shift the section arrays forward by NumSplit(ss)-2 (preserving initial and final sections)
            k = indexGC(ss);
            N = length(x_sample); % get the current length of x_sample
            x_sample((k+1:N) + NumSplit(ss)-2)        = x_sample(k+1:N);
            r_sample((k+1:N) + NumSplit(ss)-2)        = r_sample(k+1:N);
            SectionIndex((k+1:N) + NumSplit(ss)-2)    = SectionIndex(k+1:N);
            TubeIndex((k+1:N) + NumSplit(ss)-2)       = TubeIndex(k+1:N);

            if ss < length(indexGC)
                % Shift subsequent indices
                indexGC(ss+1:end) = indexGC(ss+1:end) + NumSplit(ss)-2;
            end

        % Populate x_sample and r_sample in the gap left by shift
            r1 = r_sample(k);
            rN = r_sample(k+NumSplit(ss)-1);
            x1 = x_sample(k);
            xN = x_sample(k+NumSplit(ss)-1);
        
        % Specify shape function y(z) for abstract axial co-ordinate z
            switch TubeIndex(indexGC(ss))
                case 1 % conical
                    z1 = x1;
                    zN = xN;
                    P    = @(z) z; % Shape function
                    Pinv = @(y) y; % Inverse shape function

                case 2 % exponential horn (concave outwards diverging, zero gradient LHS)
                    z1 = -5;
                    zN = 0;
                    P    = @(z) exp(z); % Shape function
                    Pinv = @(y) log(y); % Inverse shape function

                case 3 % exponential nozzle (concave outwards converging, zero gradient RHS)
                    z1 = 0;
                    zN = -5;
                    P    = @(z) exp(z); % Shape function
                    Pinv = @(y) log(y); % Inverse shape function

                case 4 % parabolic horn (concave outwards diverging, zero gradient LHS)
                    z1 = 0;
                    zN = 1;
                    P    = @(z) z.^2; % Shape function
                    Pinv = @(y) sqrt(y); % Inverse shape function

                case 5 % parabolic nozzle (concave outwards converging, zero gradient RHS)
                    z1 = 1;
                    zN = 0;
                    P    = @(z) z.^2; % Shape function
                    Pinv = @(y) sqrt(y); % Inverse shape function

                case 6 % sigmoid
                    z1 = -5;
                    zN = 5;
                    P    = @(z) 1./(1+exp(-1.*z)); % Shape function
                    Pinv = @(y) -1.*log(1./y - 1); % Inverse shape function

                case 7 % sinusoid
                    z1 = 0;
                    zN = pi;
                    P    = @(z) 1-cos(z); % Shape function
                    Pinv = @(y) acos(1-y); % Inverse shape function

                otherwise % default to conical
                    z1 = x1;
                    zN = xN;
                    P    = @(z) z; % Shape function
                    Pinv = @(y) y; % Inverse shape function
            end

            ZmapX = @(z) x1 + (xN - x1)./(zN - z1).*(z - z1); % Maps z to x
            R    = @(z) r1 + (rN - r1)./(P(zN) - P(z1)).*(P(z) - P(z1)); % Radius as a function of z
            
            z_sample = Pinv( linspace(P(z1), P(zN), 2*NumSplit(ss)-3) ); % Create z-space discretisation
            x_sample(k+1:k+NumSplit(ss)-2) = ZmapX(z_sample(2:2:length(z_sample)-1)); % Distribute x-sample using z2, z4, ... zN-3, zN-1
            r_sample(k+1:k+NumSplit(ss)-2) = R(z_sample(3:2:length(z_sample))); % Calculate radii at intermediate z-points z3, z5, ... zN-4, zN-2
            % Staggering the sampling of x and r allows the analytical line
            % to pass through the centre of the staircase, matching the
            % volumes of the staircased and analytical geometries

            SectionIndex(k+1:k+NumSplit(ss)-2) = 0; % Sections within interpolated tube have zero SectionIndex
            TubeIndex(k+1:k+NumSplit(ss)-2) = TubeIndex(indexGC(ss)); % Copy TubeIndex to new sections as a record of how they were created
            TubeIndex(k) = 0; % First section considered cylindrical - black line will be drawn
            
            % Store info required to reconstruct analytical curve later for plot
            TubeProf{ss,1} = indexGC(ss); % Initial index
            TubeProf{ss,2} = indexGC(ss) + NumSplit(ss) - 1; % Final index
            TubeProf{ss,3} = z1; % Lower z limit
            TubeProf{ss,4} = zN; % Upper z limit
            TubeProf{ss,5} = R; % Radius function
        end
    end

    %% Plotting the geometry

    if CFG.PLT.GEOM_PLOT

        fig1=figure('Name','Geometry');
        if CFG.PLT.SMALL_PLOTS
            set(fig1, 'Position', [100 50 800 300]);
        else
            set(fig1, 'Position', [100 50 1200 400]);
        end

        hold on

        WallLineWidth = 2;
        W           = abs(x_sample(end) - x_sample(1));             % Length of the combustor
        H           = 2.1*max(r_sample);                            % Diameter of the combustor
        plot_ratio  = 1.1;                                          % Ratio of the axes limit to the combustor dimension
        axes_W      = plot_ratio*W;                                 % axes x width
        axes_H      = 2*plot_ratio*H;                               % axes y width        
        x_min       = x_sample(1) - (axes_W-W)./2;                  % axes x min
        y_min       = -axes_H./2;                                   % axes y min

        % These are used for the legend
        xMax = max(get(gca,'xlim'));
        yMax = max(get(gca,'ylim'));
        plot(gca,xMax,yMax,'-b','linewidth',3)
        plot(gca,xMax,yMax,'-g','linewidth',3)
        plot(gca,xMax,yMax,'-r','linewidth',3)
        plot(gca,xMax,yMax,'-m','linewidth',3)

        % Combustor outline
        for s=1:length(x_sample)-1

            x_plot1=1000*[x_sample(s),x_sample(s+1)];
            x_plot5=1000*[x_sample(s+1),x_sample(s+1)];
            y_plot1=1000*[r_sample(s),r_sample(s)];
            y_plot2=1000*[-r_sample(s),-r_sample(s)];
            y_plot3=1000*[r_sample(s),r_sample(s+1)];
            y_plot4=1000*[-r_sample(s),-r_sample(s+1)];

            plot(x_plot1, y_plot1,'-k','linewidth',WallLineWidth);
            plot(x_plot1, y_plot2,'-k','linewidth',WallLineWidth);

            if r_sample(s) ~=r_sample(s+1)
                plot(x_plot5, y_plot3,'-k','linewidth',WallLineWidth);
                plot(x_plot5, y_plot4,'-k','linewidth',WallLineWidth);
            end

        end

        % Combustor's sections  
        for s = 1:length(SectionIndex)

            if SectionIndex(s)==0          % just interface
                indexColor = 'k';
                indexLinestyle  = '-';
                indexLineWidth  = 1;
            elseif SectionIndex(s)==10        % with heat addition but perturbation
                indexColor = 'm';
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            elseif SectionIndex(s)==11         % with heat addition and heat perturbations
                indexColor = 'r';
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            end 

            if s == 1
                indexColor = 'b';   % inlet
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            elseif s == length(SectionIndex)
                indexColor = 'g';   % outlet
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            end

            if s > 1 && s < length(SectionIndex)
                if TubeIndex(s) >= 1
                    indexLinestyle  = 'none';
                end
            end

            plot(1e3*[x_sample(s),x_sample(s)],-1e3*[-r_sample(s),r_sample(s)],...
                'linestyle',indexLinestyle,...
                'color',indexColor,'linewidth',indexLineWidth);
        end
        
        % Draw analytical curves for interpolated geometries with few steps
        for s = 1 : size(TubeProf,1)
            if (TubeProf{s,2} - TubeProf{s,1}) <=20
                x_prof = 1000.*linspace(x_sample(TubeProf{s,1}), x_sample(TubeProf{s,2}), 100);
                z_prof = linspace(TubeProf{s,3}, TubeProf{s,4}, 100);
                r_prof = 1000.*TubeProf{s,5}(z_prof);
                plot(x_prof, r_prof, '--k');
                plot(x_prof, -1*r_prof, '--k');
            end
        end

        % Handle properties
        set(gca,'xlim',1000*[x_min, x_min+axes_W]);
        set(gca,'box','on','linewidth',0.5,'gridlinestyle','-.');
        set(gca,'xgrid','on','ygrid','on');
        set(gca,'ylim',1000*[y_min, y_min+axes_H]);
        xlabel('x~ [mm]','Color','k','Interpreter','LaTex');
        ylabel('r~ [mm]','Color','k','Interpreter','LaTex');

        % Legend properties
        %legend1 = ['Mean heat addition', newline, ' and heat perturbations'];
        %legend2 = ['Mean heat addition', newline, 'but no heat perturbation'];
        hlegend = legend('Inlet','Outlet'); %,legend1,legend2);
        set(hlegend,'location','eastoutside');

        if ~CFG.PLT.SMALL_PLOTS
            set(gca,'fontsize',16);
            set(hlegend,'fontsize',14);
        end

        % Saving the figure
        flpth = fullfile(CFG.pthstr,'Outputs',CFG.OUT.RUN_NAME,'Initialization','Geometry');
        if CFG.PLT.SAVE_FIGS
            saveas(fig1,flpth,'fig')
        end
        if CFG.PLT.SAVE_PDFS
            save2pdf(flpth,fig1,300)
        end
        drawnow
    end

    if CFG.OUT.CL_OUT 
        fprintf("done\n");
        fprintf("\t\t%d sections listed, %d generated\n",length(TubeSplit),length(TubeIndex));
    end
    if CFG.OUT.LOG_OUT 
        write_log(sprintf("\t\t%d sections listed, %d generated",length(TubeSplit),length(TubeIndex)),CFG.OUT.log_filename);
    end

    GEOM.x_sample = x_sample;
    GEOM.r_sample = r_sample;
    GEOM.SectionIndex = SectionIndex;
    GEOM.TubeIndex = TubeIndex;
    GEOM.TubeSplit = TubeSplit;

end
