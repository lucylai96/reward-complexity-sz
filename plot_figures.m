function plot_figures(fig,results,data)
    
    % Plot figures.
    %
    % USAGE: plot_figures(fig,[results],[data])
    %
    % INPUTS:
    %   fig - {'fig2', 'fig3', 'fig4'}
    
    if nargin < 2; load results_collins14.mat; end
    if nargin < 3; data = load_data; end
    
    C = linspecer(2);   % colors
    
    switch fig
        
        case 'fig2'
            
            figure;
            T = {'A' 'B' 'C' 'D' 'E'};
            R = squeeze(nanmean(results.R));
            V = squeeze(nanmean(results.V));
            ylim = [0.25 1.1];
            xlim = [0 0.9];
            cond = [data.cond];
            for j = 1:size(R,2)
                subplot(2,3,j);
                h(1) = plot(R(:,j),V(:,j),'-k','LineWidth',4);
                hold on
                xlabel('Policy complexity','FontSize',25);
                ylabel('Average reward','FontSize',25);
                set(gca,'FontSize',25,'YLim',ylim,'XLim',xlim);
                for i = 1:2
                    ix = cond==i-1;
                    h(i+1) = plot(results.R_data(ix,j),results.V_data(ix,j),'o','Color',C(i,:),'MarkerSize',10,'LineWidth',3,'MarkerFaceColor',C(i,:));
                end
                if j==1
                    legend(h,{'Theory' 'HC' 'SZ'},'FontSize',25,'Location','SouthEast');
                end
                mytitle([T{j},')   Set size: ',num2str(j+1)],'Left','FontSize',30,'FontWeight','Bold');
            end
            
            subplot(2,3,6);
            x = 2:6;
            for i=1:2
                [mu,~,ci] = normfit(results.R_data(cond==i-1,:));
                err = diff(ci)/2;
                errorbar(x',mu,err,'-o','Color',C(i,:),'MarkerSize',10,'LineWidth',4,'MarkerFaceColor',C(i,:));
                hold on;
            end
            set(gca,'FontSize',25,'XLim',[1.5 6.5],'XTick',2:6);
            ylabel('Policy complexity','FontSize',25);
            xlabel('Set size','FontSize',25);
            mytitle('F)','Left','FontSize',30,'FontWeight','Bold');
            
            set(gcf,'Position',[200 200 1200 800])
            
        case 'fig3'
            
            figure;
            
            cond = [data.cond];
            for i = 1:size(results.bias,2)
                for j = 1:2
                    [m(i,j),~,ci] = normfit(results.bias(cond==j-1,i));
                    err(i,j) = diff(ci)/2;
                end
            end
            
            subplot(1,2,1);
            x = 2:6;
            for i = 1:2
                errorbar(x',m(:,i),err(:,i),'-o','Color',C(i,:),'MarkerSize',10,'LineWidth',4,'MarkerFaceColor',C(i,:));
                hold on;
            end
            legend({'HC' 'SZ'},'FontSize',25,'Location','NorthWest');
            set(gca,'FontSize',25,'XLim',[1.5 6.5],'XTick',2:6);
            xlabel('Set size','FontSize',25);
            ylabel('Bias','FontSize',25);
            mytitle('A)','Left','FontSize',30,'FontWeight','Bold');
            
            T = {'HC' 'SZ'};
            subplot(1,2,2);
            for j = 1:2
                y = results.bias(cond==j-1,:);
                x = results.R_data(cond==j-1,:);
                plot(x(:),y(:),'o','Color',C(j,:),'MarkerSize',10,'LineWidth',3,'MarkerFaceColor',C(j,:));
                H = lsline; set(H,'LineWidth',4);
                hold on;
                [r,p,rl,ru] = corrcoef(x(:),y(:));
                disp([T{j},': r = ',num2str(r(2,1)),', p = ',num2str(p(2,1)),', CI = [',num2str(rl(2,1)),',',num2str(ru(2,1)),']']);
                [r,p] = corr(x(:),y(:),'type','spearman')
            end
            mytitle('B)','Left','FontSize',30,'FontWeight','Bold');
            set(gca,'FontSize',25);
            xlabel('Policy complexity','FontSize',25);
            ylabel('Bias','FontSize',25);
            
            set(gcf,'Position',[200 200 1200 400])
            
        case 'fig4'
            
            T = {'A' 'B' 'C' 'D' 'E'};
            for c = 1:5
                subplot(2,3,c)
                m = squeeze(results.b_sep(c,:,:));
                err = squeeze(results.bci_sep(c,:,:));
                h = barerrorbar(m',err');
                h(1).CData = C(1,:); h(2).CData = C(2,:);
                ylabel('Parameter value','FontSize',25);
                set(gca,'XTickLabel',{'\beta_0' '\beta_1' '\beta_2'},'FontSize',25,'YLim',[-3 3]);
                mytitle([T{c},')   Set size: ',num2str(c+1)],'Left','FontSize',30,'FontWeight','Bold');
                
                if c==1
                    legend(h,{'HC' 'SZ'},'FontSize',25,'Location','South');
                end
            end
            
            subplot(2,3,6);
            C = linspecer(3);   % colors
            h = bar(results.bic(:,1)-results.bic(:,2)); set(h,'FaceColor',C(3,:))
            set(gca,'FontSize',25);
            xlabel('Set size','FontSize',25);
            ylabel('\Delta BIC','FontSize',25);
            mytitle('F)','Left','FontSize',30,'FontWeight','Bold');
            
            set(gcf,'Position',[200 200 1200 800])
            
        case 'fig7'
            
            load results_collins14.mat
            R = results;
            load model_fits
            
            figure;
            T = {'A' 'B' 'C'};
            L = {'\beta' '\alpha_\theta' '\alpha_P'};
            bias = mean(R.bias,2);
            ylim = [0 0.3];
            cond = [data.cond];
            q = [1 2 4];
            for j = 1:length(q)
                subplot(1,3,j);
                for i = 1:2
                    ix = cond==i-1;
                    plot(results.x(ix,q(j)),bias(ix),'o','Color',C(i,:),'MarkerSize',10,'LineWidth',3,'MarkerFaceColor',C(i,:));
                    hold on;
                end
                H = lsline; set(H,'LineWidth',4);
                if j==1
                    legend({'HC' 'SZ'},'FontSize',25,'Location','NorthEast');
                end
                mytitle(T{j},'Left','FontSize',30,'FontWeight','Bold');
                xlabel(L{j},'FontSize',25);
                ylabel('Bias','FontSize',25);
                set(gca,'FontSize',25,'YLim',ylim);
            end
            
            set(gcf,'Position',[200 200 1200 400]);
            
            bias = mean(R.bias,2);
            [b,~,stats] = glmfit(results.x,bias,'normal');
            
            [b,~,stats] = glmfit(results.x(c==1,:),bias(c==1),'normal');
            [b,~,stats] = glmfit(results.x(c==0,:),bias(c==0),'normal');
            
            [~,p,~,stats] = ttest2(results.x(c==0,:),results.x(c==1,:))
            
        case 'fig8'
            
            
            data = load_data;
            load simresults_m3
            
            figure;
            for s = 1:length(data)
                for k = 2:6
                    for i = 1:7
                        ix = data(s).iter==i & data(s).ns==k;
                        r_data(s,i,k-1) = mean(data(s).reward(ix));
                        r_model(s,i,k-1) = mean(simdata(s).expreward(ix));
                    end
                end
            end
            
            cond = [data.cond];
            for i = 1:2
                m{i} = squeeze(mean(r_data(cond==i-1,:,:)));
                mm{i} = squeeze(mean(r_model(cond==i-1,:,:)));
            end
            
            T = {'A) HC: data' 'B) SZ: data' 'C) HC: model' 'D) SZ: model'};
            C = linspecer(5);
            
            subplot(2,2,1);
            for i = 1:5
                plot(m{1}(:,i),'-','LineWidth',4,'Color',C(i,:));
                set(gca,'FontSize',25,'XLim',[0 8],'XTick',1:7,'YLim',[0 1]);
                xlabel('Iteration','FontSize',25);
                ylabel('Average reward','FontSize',25);
                mytitle(T{1},'Left','FontSize',30,'FontWeight','Bold');
                hold on;
            end
            L = legend({'2' '3' '4' '5' '6'},'FontSize',25,'Location','SouthEast');
            title(L,'Set size','FontSize',25);
            
            subplot(2,2,2);
            for i = 1:5
                plot(m{2}(:,i),'-','LineWidth',4,'Color',C(i,:));
                set(gca,'FontSize',25,'XLim',[0 8],'XTick',1:7,'YLim',[0 1]);
                xlabel('Iteration','FontSize',25);
                ylabel('Average reward','FontSize',25);
                mytitle(T{2},'Left','FontSize',30,'FontWeight','Bold');
                hold on;
            end
            
            subplot(2,2,3);
            for i = 1:5
                plot(mm{1}(:,i),'-','LineWidth',4,'Color',C(i,:));
                set(gca,'FontSize',25,'XLim',[0 8],'XTick',1:7,'YLim',[0 1]);
                xlabel('Iteration','FontSize',25);
                ylabel('Average reward','FontSize',25);
                mytitle(T{3},'Left','FontSize',30,'FontWeight','Bold');
                hold on;
            end
            
            subplot(2,2,4);
            for i = 1:5
                plot(mm{2}(:,i),'-','LineWidth',4,'Color',C(i,:));
                set(gca,'FontSize',25,'XLim',[0 8],'XTick',1:7,'YLim',[0 1]);
                xlabel('Iteration','FontSize',25);
                ylabel('Average reward','FontSize',25);
                mytitle(T{4},'Left','FontSize',30,'FontWeight','Bold');
                hold on;
            end
            
            set(gcf,'Position',[200 200 800 800])
            
        case 'fig9'
            
            load results_collins14
            
            beta = linspace(0.1,15,50);
            R = squeeze(mean(results.R));
            m = mean(results.R_data);
            
            C = linspecer(5);
            for i = 1:5
                h(i) = plot(R(:,i),beta','-','LineWidth',4,'Color',C(i,:));
                hold on;
                ix = find(R(:,i)>=m(i),1);
                plot(R(ix,i),beta(ix),'o','LineWidth',4,'Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',12);
            end
            lgd = legend(h,{'2' '3' '4' '5' '6'},'FontSize',25,'Location','NorthWest');
            title(lgd,'Set size','FontSize',25);
            set(gca,'FontSize',25);
            xlabel('Policy complexity','FontSize',25);
            ylabel('\beta','FontSize',25);
            
            for s=1:length(data); for k=2:6; beta(s,k-1) = mean(simdata(s).beta(simdata(s).ns==k)); end; end
            
        case 'RT'
            
            cond = [data.cond];
            ix = find(cond==1);
            for s=1:length(ix)
                for j=2:6
                    rt(s,j-1) = mean(data(ix(s)).rt(data(ix(s)).ns==j));
                end
            end
            
            [se,m] = wse(rt);
            errorbar(log(2:6)',m,se,'-k','LineWidth',4)
            set(gca,'FontSize',25);
            ylabel('Response time (sec)','FontSize',25);
            xlabel('Set size (log)','FontSize',25);
            
        case 'entropy'
            
            cond = [data.cond];
            ix = find(cond==1);
            for s=1:length(ix)
                B = unique(data(s).learningblock);
                H = zeros(length(B),1);
                setsize = zeros(length(B),1);
                for b = 1:length(B)
                    ix = data(s).learningblock==B(b);
                    state = data(s).state(ix);
                    action = data(s).action(ix);
                    ns = data(s).ns(ix); ns = ns(1);
                    p = zeros(ns,3);
                    for i = 1:ns
                        for a = 1:3
                            p(i,a) = mean(action(state==i)==a);
                        end
                    end
                    H(b) = mean(sum(-safelog(p).*p,2));
                    setsize(b) = ns;
                end
                
                for i = 2:6
                    h(s,i-1) = mean(H(setsize==i));
                end
            end
            
            [se,m] = wse(h);
            errorbar(2:6,m,se,'-k','LineWidth',4);
            set(gca,'FontSize',25,'XLim',[1 7],'XTick',2:6);
            ylabel('H(A|S)','FontSize',25);
            xlabel('Set size','FontSize',25);
            
    end
    
end

function lineStyles=linspecer(N,varargin)
    
    % function lineStyles = linspecer(N)
    % This function creates an Nx3 array of N [R B G] colors
    % These can be used to plot lots of lines with distinguishable and nice
    % looking colors.
    %
    % lineStyles = linspecer(N);  makes N colors for you to use: lineStyles(ii,:)
    %
    % colormap(linspecer); set your colormap to have easily distinguishable
    %                      colors and a pleasing aesthetic
    %
    % lineStyles = linspecer(N,'qualitative'); forces the colors to all be distinguishable (up to 12)
    % lineStyles = linspecer(N,'sequential'); forces the colors to vary along a spectrum
    %
    % % Examples demonstrating the colors.
    %
    % LINE COLORS
    % N=6;
    % X = linspace(0,pi*3,1000);
    % Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);
    % C = linspecer(N);
    % axes('NextPlot','replacechildren', 'ColorOrder',C);
    % plot(X,Y,'linewidth',5)
    % ylim([-1.1 1.1]);
    %
    % SIMPLER LINE COLOR EXAMPLE
    % N = 6; X = linspace(0,pi*3,1000);
    % C = linspecer(N)
    % hold off;
    % for ii=1:N
    %     Y = sin(X+2*ii*pi/N);
    %     plot(X,Y,'color',C(ii,:),'linewidth',3);
    %     hold on;
    % end
    %
    % COLORMAP EXAMPLE
    % A = rand(15);
    % figure; imagesc(A); % default colormap
    % figure; imagesc(A); colormap(linspecer); % linspecer colormap
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % by Jonathan Lansey, March 2009-2013 – Lansey at gmail.com               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % credits and where the function came from
    % The colors are largely taken from:
    % http://colorbrewer2.org and Cynthia Brewer, Mark Harrower and The Pennsylvania State University
    %
    %
    % She studied this from a phsychometric perspective and crafted the colors
    % beautifully.
    %
    % I made choices from the many there to decide the nicest once for plotting
    % lines in Matlab. I also made a small change to one of the colors I
    % thought was a bit too bright. In addition some interpolation is going on
    % for the sequential line styles.
    %
    %
    %
    
    if nargin==0 % return a colormap
        lineStyles = linspecer(64);
        %     temp = [temp{:}];
        %     lineStyles = reshape(temp,3,255)';
        return;
    end
    
    if N<=0 % its empty, nothing else to do here
        lineStyles=[];
        return;
    end
    
    % interperet varagin
    qualFlag = 0;
    
    if ~isempty(varargin)>0 % you set a parameter?
        switch lower(varargin{1})
            case {'qualitative','qua'}
                if N>12 % go home, you just can't get this.
                    warning('qualitiative is not possible for greater than 12 items, please reconsider');
                else
                    if N>9
                        warning(['Default may be nicer for ' num2str(N) ' for clearer colors use: whitebg(''black''); ']);
                    end
                end
                qualFlag = 1;
            case {'sequential','seq'}
                lineStyles = colorm(N);
                return;
            otherwise
                warning(['parameter ''' varargin{1} ''' not recognized']);
        end
    end
    
    % predefine some colormaps
    set3 = colorBrew2mat({[141, 211, 199];[ 255, 237, 111];[ 190, 186, 218];[ 251, 128, 114];[ 128, 177, 211];[ 253, 180, 98];...
        [ 179, 222, 105];[ 188, 128, 189];[ 217, 217, 217];[ 204, 235, 197];[ 252, 205, 229];[ 255, 255, 179]}');
    set1JL = brighten(colorBrew2mat({[228, 26, 28];[ 55, 126, 184];[ 77, 175, 74];[ 255, 127, 0];[ 255, 237, 111]*.95;[ 166, 86, 40];...
        [ 247, 129, 191];[ 153, 153, 153];[ 152, 78, 163]}'));
    set1 = brighten(colorBrew2mat({[ 55, 126, 184]*.95;[228, 26, 28];[ 77, 175, 74];[ 255, 127, 0];[ 152, 78, 163]}),.8);
    
    set3 = dim(set3,.93);
    
    switch N
        case 1
            lineStyles = { [  55, 126, 184]/255};
        case {2, 3, 4, 5 }
            lineStyles = set1(1:N);
        case {6 , 7, 8, 9}
            lineStyles = set1JL(1:N)';
        case {10, 11, 12}
            if qualFlag % force qualitative graphs
                lineStyles = set3(1:N)';
            else % 10 is a good number to start with the sequential ones.
                lineStyles = cmap2linspecer(colorm(N));
            end
        otherwise % any old case where I need a quick job done.
            lineStyles = cmap2linspecer(colorm(N));
    end
    lineStyles = cell2mat(lineStyles);
end

% extra functions
function varIn = colorBrew2mat(varIn)
    for ii=1:length(varIn) % just divide by 255
        varIn{ii}=varIn{ii}/255;
    end
end

function varIn = brighten(varIn,varargin) % increase the brightness
    
    if isempty(varargin),
        frac = .9;
    else
        frac = varargin{1};
    end
    
    for ii=1:length(varIn)
        varIn{ii}=varIn{ii}*frac+(1-frac);
    end
end

function varIn = dim(varIn,f)
    for ii=1:length(varIn)
        varIn{ii} = f*varIn{ii};
    end
end

function vOut = cmap2linspecer(vIn) % changes the format from a double array to a cell array with the right format
    vOut = cell(size(vIn,1),1);
    for ii=1:size(vIn,1)
        vOut{ii} = vIn(ii,:);
    end
end

% colorm returns a colormap which is really good for creating informative
% heatmap style figures.
% No particular color stands out and it doesn't do too badly for colorblind people either.
% It works by interpolating the data from the
% 'spectral' setting on http://colorbrewer2.org/ set to 11 colors
% It is modified a little to make the brightest yellow a little less bright.
function cmap = colorm(varargin)
    n = 100;
    if ~isempty(varargin)
        n = varargin{1};
    end
    
    if n==1
        cmap =  [0.2005    0.5593    0.7380];
        return;
    end
    if n==2
        cmap =  [0.2005    0.5593    0.7380;
            0.9684    0.4799    0.2723];
        return;
    end
    
    frac=.95; % Slight modification from colorbrewer here to make the yellows in the center just a bit darker
    cmapp = [158, 1, 66; 213, 62, 79; 244, 109, 67; 253, 174, 97; 254, 224, 139; 255*frac, 255*frac, 191*frac; 230, 245, 152;...
        171, 221, 164; 102, 194, 165; 50, 136, 189; 94, 79, 162];
    x = linspace(1,n,size(cmapp,1));
    xi = 1:n;
    cmap = zeros(n,3);
    for ii=1:3
        cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
    end
    cmap = flipud(cmap/255);
end

function h = mytitle(txt,just,varargin)
    
    % Plot titles with left or right justification.
    %
    % USAGE: h = mytitle(txt,[just],[varargin])
    %
    % INPUTS:
    %   txt - title string
    %   just (optional) - justification: 'Left' (default) or 'Right'
    %   varargin (optional) - additional arguments for title
    %
    % OUTPUTS:
    %   h - handle to title object
    %
    % Sam Gershman, Dec 2011
    
    if nargin < 2 || isempty(just)
        just = 'Left';
    end
    
    switch just
        case 'Left'
            h = title(txt,'HorizontalAlignment','left','Units','normalized','Position',[0 1],varargin{:});
        case 'Right'
            h = title(txt,'HorizontalAlignment','right','Units','normalized','Position',[1 1],varargin{:});
    end
    
end

function varargout = barerrorbar(varargin)
    % BARERRORBAR   Create a bar plot with error bars. BARERRORBAR() uses the
    %   MATLAB functions BAR() and ERRORBAR() and changes the 'XData' property
    %   of the errorbar plot so that the error bars are plotted at the center
    %   of each bar. This does not support "stack" bar plots.
    %
    % Syntax: varargout = barerrorbar(varargin)
    %
    % Inputs:
    %   Method 1: Cell array method
    %       varargin{1} - A cell array containing all of the inputs to be fed
    %           to bar().
    %       varargin{2} - A cell array containing all of the inputs to be fed
    %           to errorbar().
    %   Method 2: Simple Method
    %       varargin{1} - The data to be plotted by bar().
    %       varargin{2} - The E input to errorbar().
    %
    % Outputs:
    %   varargout{1} - The handle to the bar plot.
    %   varargout{2} - The handle to the errorbar plot.
    %
    % Examples:
    %   x = 0.2*randn(3,4,100) + 1;
    %   xMeans = mean(x,3);
    %   xMeansConf = repmat(2*0.2/10,3,4);
    %   xMeansL = repmat(3*0.2/10,3,4);
    %   xMeansU = repmat(4*0.2/10,3,4);
    %
    %   figure
    %   barerrorbar(xMeans,xMeansConf);
    %
    %   figure
    %   barerrorbar({3:5,xMeans,'m'}, {repmat((3:5)',1,4),xMeans, xMeansL,xMeansU,'bx'});
    %
    %   figure
    %   barerrorbar({3:5,xMeans,'k'}, {repmat((3:5)',1,4),xMeans, xMeansL,xMeansU,'bx'});
    %   hold on
    %   barerrorbar({7:9,xMeans}, {repmat((7:9)',1,4),xMeans, 2*xMeansL,4*xMeansU,'d','color',[0.7 0.7 0.7]});
    %
    % Other m-files required: none
    % Subfunctions: interpretinputs
    % MAT-files required: none
    %
    % See also: bar.m errorbar.m
    
    % Author: Kenneth D. Morton Jr. and J. Simeon Stohl
    % Revised by: Kenneth D. Morton Jr.
    % Duke University, Department of Electrical and Computer Engineering
    % Email Address: kennethmorton@ieee.org
    % Created: 07-July-2005
    % Last revision: 06-January-2006
    
    % Check if hold is on
    startedWithHoldOn = ishold;
    
    [data, barInputs, errorbarInputs] = interpinputs(varargin);
    
    % Create the bar plot and keep the handles for later use.
    barHandles = bar(barInputs{:},'FaceColor','flat');
    
    for k = 1:length(barHandles)
        barHandles(1).CData = [0.7 0.7 0.7];
        barHandles(2).CData = [0.2 0.2 0.2];
    end
    barHandlesStruct = get(barHandles);
    
    if ~startedWithHoldOn
        hold on
    else
        % Hold is already on so we need to check the XTick and save them.
        oldXTicks = get(gca,'XTick');
        oldXTickLabels = get(gca,'XTickLabel');
    end
    
    % Find out the bar width which is dependent on the number of bar sets and
    % the number of bars in each set.
    barWidth = barHandlesStruct(1).BarWidth;
    
    errorbarHandles = errorbar(errorbarInputs{:});
    
    % The crazy stuff to find the bar centers. Some of it is how bar.m does it.
    [nGroups, nBarsPerGroup] = size(data);
    for iBpg = 1:nBarsPerGroup
        groupMembership(:,iBpg) = barHandlesStruct(iBpg).XData;
    end
    groupWidth = min(0.8, nBarsPerGroup/(nBarsPerGroup+1.5));
    groupLocs = groupMembership(:,1);
    distanceToNearestGroup = zeros(1,nGroups);
    if nGroups > 1
        for iGroup = 1:nGroups
            distanceToNearestGroup(iGroup) = ...
                min(abs(groupLocs(iGroup)-...
                groupLocs(groupLocs~=groupLocs(iGroup))));
        end
    end
    groupWidth = groupWidth*min(distanceToNearestGroup);
    barGap = (nGroups - 1) * groupWidth / (nGroups-1) / nBarsPerGroup;
    almostCenters  = (0:nBarsPerGroup-1)'.*barGap - 0.5*barGap*nBarsPerGroup;
    relativeCenters = almostCenters + mean([(1-barWidth)/2.*barGap; (1+barWidth)/2.*barGap]);
    
    centers = repmat(relativeCenters',nGroups,1) + groupMembership;
    % Change the XData of the errorbars to be at our bar centers.
    for iBpg = 1:nBarsPerGroup
        set(errorbarHandles(iBpg),'XData',centers(:,iBpg));
    end
    
    % Turn hold off if it wasn't on to begin with
    if ~startedWithHoldOn
        hold off
    else
        % Set the XTick and XTickLabels to inlcude the old and new information.
        newXTicks = groupMembership(:,1);
        newXTickLabels = num2str(newXTicks);
        set(gca,'XTick',sort(unique([oldXTicks newXTicks'])));
        if ischar(oldXTickLabels)
            % If this is a string then they are probably default so update with
            % the new additions.
            set(gca,'XTickLabel',[oldXTickLabels; newXTickLabels]);
        end
    end
    
    % Prepare both sets of handles as outputs, if requested.
    if nargout > 0
        varargout{1} = barHandles;
    end
    if nargout > 1
        varargout{2} = errorbarHandles;
    end
    
end

function [data, barInputs, errorbarInputs] = interpinputs(inCell)
    % Interperate the inputs so that they can be used in barerrorbar(). Make
    % two different possibilities based on the type of inputs.
    %   Method 1: Cell array method.
    %       varargin{1} - A cell array containing all of the inputs to be fed
    %           to bar().
    %       varargin{2} - A cell array containing all of the inputs to be fed
    %           to errorbar().
    %   Method 2: Simple Method
    %       varargin{1} - The data to be plotted by bar().
    %       varargin{2} - The data to be plotted by errorbar().
    
    if iscell(inCell{1})
        % We have entered Method 1
        barInputs = inCell{1};
        if length(barInputs) > 2
            data = barInputs{2};
        elseif length(barInputs) < 2
            data = barInputs{1};
        elseif length(barInputs) == 2
            if ischar(barInputs{2})
                data = barInputs{1};
            else
                data = barInputs{2};
            end
        end
    else
        barInputs = {inCell{1}};
        data = inCell{1};
        nRows = size(barInputs{1},1);
        if nRows == 1
            barInputs{1} = barInputs{1}';
        end
    end
    
    %Plot black dot errorbars for the default.
    if iscell(inCell{2})
        errorbarInputs = inCell{2};
    else
        errorbarInputs = {inCell{1}, inCell{2}, 'k.','LineWidth',2};
    end
    
    if length(inCell) > 2
        error(['Too many input arguments.' ...
            ' Try using the cell array input method.']);
    end
    
    % Search for the 'stack' option in the bar inputs
    for iBI = 1:length(barInputs)
        if ischar(barInputs{iBI})
            if strcmpi(barInputs{iBI},'stack')
                error('barerrorbar() does not support "stack" bar plots.')
            end
        end
    end
end