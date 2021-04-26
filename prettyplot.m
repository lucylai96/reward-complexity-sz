function prettyplot(varargin)
% purpose: makes things aestheticly pleasing
% written by: lucy lai
% todo: add defaults later (like default font size and stuff, if not
% specified

%varargin: 1-fontsize
for i = 1:length(varargin)  
    if i ==1
        fz = varargin{1};
    elseif i ==2
    end
end

set(gcf,'color','w'); % make matlab plot background white
%set(gca,'LineWidth',2); 

if exist('fz')
set(gca,'FontSize',fz); % make font bigger (input can be font size)
else
    set(gca,'FontSize',15); 
end

box off


end