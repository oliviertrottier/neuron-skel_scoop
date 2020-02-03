function hfig = tightax(hfig)
% Function to resize the axis handle of a figure to fill the entire figure
% without excess space around them. The figure size is kept constant and 
% the axes are enlarged to fit the figure.

% The colorbars are moved, after finding their tightinset.

% Inspired from tightfig.

% Input
% hfig = handle to figure. If not supplied, the current figure is used.
if nargin == 0
    hfig = gcf;
end

% Process all graphics callbacks.
drawnow;
pause(1);

origwindowstyle = get(hfig, 'WindowStyle');
set(hfig, 'WindowStyle', 'normal');
%% Find all the colorbars and overlay a phantom axis on top of them to measure their tight inset.
cbarhandles = findall(hfig,'type', 'ColorBar');
N_cbars = numel(cbarhandles);
cbars_pos = zeros(N_cbars,4);
cbars_ti = zeros(N_cbars,4);
cbars_tightpos = zeros(N_cbars,4);
for i=1:N_cbars
    cbars_pos(i,:) = cbarhandles(i).Position;
    cbars_ti(i,:) = get_colorbar_tightinset(cbarhandles(i));
    
    % If the colorbar tight width is out of bounds, reduce the width of the
    % associated axis. This is necessary because the width of colorbars is
    % not changed at the end. Therefore, colorbars need to fit completely
    % in the figure, before tight measurements are taken.
    width_overflow = cbars_pos(i,1) + cbars_ti(i,3) - 1;
    if width_overflow > 0
        peer_axis_h = get_colorbar_peer_axis(cbarhandles(i));
        axis_width = peer_axis_h.Position(3);
        resize_factor = 0.95*(axis_width - width_overflow)/axis_width;
        peer_axis_h.Position(3:4) = peer_axis_h.Position(3:4)*resize_factor;
        
        % Record the new colorbar position.
        cbars_pos(i,:) = cbarhandles(i).Position;
        cbars_ti(i,:) = get_colorbar_tightinset(cbarhandles(i));
    end
    
    % Find the colorbar tight position.
    cbars_tightpos(i,1:2) = cbars_pos(i,1:2) - cbars_ti(i,1:2);
    cbars_tightpos(i,3:4) = cbars_pos(i,3:4) + [cbars_ti(i,1) + cbars_ti(i,3) cbars_ti(i,2) + cbars_ti(i,4)];
end
%%
% Get all axes.
hax = findall(hfig, 'type', 'axes');

% Get positions of the axes.
N_axes = numel(hax);
axes_ti = zeros(N_axes,4);
axes_pos = zeros(N_axes,4);
axes_tightpos = zeros(N_axes,4);

% Margin size in normalized units.
%Margin_size = 0.025;
Margin_size = 0.01;
for i=1:N_axes
    axes_pos(i,:) = get(hax(i), 'Position');
    if strcmp(hax(i).Visible,'on')
        axes_ti(i,:) = get(hax(i),'TightInset');
    else
        axes_ti(i,:) = zeros(1,4);
    end
    % Increase the borders slightly to make sure that no labels are
    % cropped.
    axes_ti(i,:) = axes_ti(i,:) + Margin_size;

    % Find the axis tight position.
    axes_tightpos(i,1:2) = axes_pos(i,1:2) - axes_ti(i,1:2);
    axes_tightpos(i,3:4) = axes_pos(i,3:4) + [axes_ti(i,1) + axes_ti(i,3) axes_ti(i,2) + axes_ti(i,4)];
end

% Ensure very tiny border so outer box always appears
% axes_ti_min = 0.5;
% axes_ti(axes_ti < axes_ti_min) = axes_ti_min;

% we will check if any 3d axes are zoomed, to do this we will check if
% they are not being viewed in any of the 2d directions
views2d = [0,90; 0,0; 90,0];

for i = 1:numel(hax)
    % Get the current viewing angle of the axes.
    [az,el] = view(hax(i));
    
    % Determine if the axes are zoomed.
    iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
    
    % Test if we are viewing in 2d mode or a 3d view
    is2d = all(bsxfun(@eq, [az,el], views2d), 2);
    
    if iszoomed && ~any(is2d)
        error('Cannot make figures containing zoomed 3D axes tight.')
    end
end
%% Add colorbars to the axis handles array.
hax = [hax(:); cbarhandles(:)];
N_axes = numel(hax);
axes_ti = [axes_ti; cbars_ti];
axes_pos = [axes_pos; cbars_pos];
axes_tightpos = [axes_tightpos; cbars_tightpos];
%% Move all axes towards the bottom left corner.
axis_displacement = min(axes_tightpos(:,1:2),[],1);

for i = 1:N_axes
    % Change axes' positions.
    axes_pos(i,1:2) = axes_pos(i,1:2) - axis_displacement;
    axes_tightpos(i,1:2) = axes_pos(i,1:2) - axes_ti(i,1:2);
    axes_tightpos(i,3:4) = axes_pos(i,3:4) + [axes_ti(i,1) + axes_ti(i,3) axes_ti(i,2) + axes_ti(i,4)];
end

% Determine the amount of non-strechable space in each dimension. This is
% the sum of all the tightinsets (or margins).
NonStretchable_size = [sum([axes_ti(:,1);axes_ti(:,3)]), sum([axes_ti(:,2);axes_ti(:,4)])];

% Determine the width and height expansion ratio.
fig_tight_width = max(axes_tightpos(:,1) + axes_tightpos(:,3)) - min(axes_tightpos(:,1));
fig_tight_height = max(axes_tightpos(:,2) + axes_tightpos(:,4)) - min(axes_tightpos(:,2));
%exp_ratios = 1./[fig_tight_width fig_tight_height];
exp_ratios = (1 - NonStretchable_size)./([fig_tight_width fig_tight_height] - NonStretchable_size);

% Determine the new position of the axes.
for i = 1:numel(hax)
    %axes_outerpos(i,1:2) = exp_ratios.*(axes_outerpos(i,1:2) - axes_outerpos_min(1:2)) + axes_outerpos_min(1:2);
    axes_tightpos(i,1:2) = exp_ratios.*axes_tightpos(i,1:2);
    axes_tightpos(i,3:4) = exp_ratios.*axes_tightpos(i,3:4);
    axes_pos(i,1:2) = axes_tightpos(i,1:2) + axes_ti(i,1:2);
    axes_pos(i,3:4) = axes_tightpos(i,3:4) - [axes_ti(i,1) + axes_ti(i,3) axes_ti(i,2) + axes_ti(i,4)];
end

% Change the size of all axes and colorbars.
ax_handles = findall(gcf,'type','axes');
ax_colorbars = [ax_handles.Colorbar];
for i = 1:numel(hax)
    if ~isa(hax(i),'matlab.graphics.illustration.ColorBar')
        hax(i).Position = axes_pos(i,:);
    else
        peer_ax = ax_handles(ax_colorbars==hax(i));
        % Do not change the width of colorbars.
        colorbar_loc = hax(i).Location;
        %hax(i).Position([4]) = axes_pos(i,[4]);
        hax(i).Position(1) = axes_pos(i,1);
        hax(i).Position(3) = axes_pos(i,3);
        hax(i).Position([2 4]) = peer_ax.Position([2 4]);
        hax(i).Location = colorbar_loc;
    end
end

set(hfig, 'WindowStyle', origwindowstyle);
end
%% Function to find the tight inset of a colorbar.
function tightinset = get_colorbar_tightinset(cbar_handles)
N_cbar = numel(cbar_handles);
tightinset = zeros(N_cbar,4);

% Add a phantom axis for measuring text labels.
Parent_fig = cbar_handles(1).Parent;
phantom_axis_h = axes(Parent_fig,'Visible','off');
for i=1:N_cbar
    % Initialize the bounding box positions ([xmin ymin xmax ymax]).
    cbar_pos = cbar_handles(i).Position;
    bounding_box_pos = [cbar_pos(1:2) cbar_pos(1:2) + cbar_pos(3:4)];
    
    % Find the position of the colorbar label in units of the parent figure.
    label =  cbar_handles(i).Label;
    label_units = label.Units;
    label.Units = 'normalized';
    label_pos = label.Extent.*repmat(cbar_pos(3:4),1,2) + [cbar_pos(1:2) 0 0];
    cbar_handles(i).Label.Units = label_units;
    bounding_box_pos(1:2) = min([bounding_box_pos(1:2);label_pos(1:2)],[],1);
    bounding_box_pos(3:4) = max([bounding_box_pos(3:4);label_pos(1:2) + label_pos(3:4)],[],1);
    
    % Find the position of the tick labels.
    % Add a textbox for each tick label to measure their extent.
    N_ticks = numel(cbar_handles.Ticks);
    Tick_Pos = zeros(N_ticks,4);
    cbar_limits = cbar_handles(i).Limits;
    
    % Determine if the ticks are on the left/right of the colorbar.
    are_ticks_outside = contains(cbar_handles(i).Location,'outside');
    for j=1:N_ticks
        % Insert textboxes on tick marks to measure their extent.
        ticklabel_textbox = text(phantom_axis_h,0,0,cbar_handles(i).TickLabels(j),'Units','normalized');
        ticklabel_textbox.FontSize = cbar_handles(i).FontSize;
        
        Tick_norm_pos = (cbar_handles(i).Ticks(j)-cbar_limits(1))/(cbar_limits(2)-cbar_limits(1));
        if are_ticks_outside
            Tick_Pos(j,1) = cbar_pos(1) + cbar_pos(3);
        else
            Tick_Pos(j,1) = cbar_pos(1) - ticklabel_textbox.Extent(3)*cbar_pos(3);
        end
        
        Tick_Pos(j,2) = cbar_pos(2) + Tick_norm_pos*cbar_pos(4);
        Tick_Pos(j,3:4) = ticklabel_textbox.Extent(3:4).*cbar_pos(3:4);
    end
    bounding_box_pos(1:2) = min([bounding_box_pos(1:2);Tick_Pos(:,1:2)],[],1);
    bounding_box_pos(3:4) = max([bounding_box_pos(3:4);Tick_Pos(:,1:2) + Tick_Pos(:,3:4)],[],1);
    
    % Define the tightinset from the bounding box.
    tightinset(i,1:2) = cbar_pos(1:2) - bounding_box_pos(1:2);
    tightinset(i,3) = bounding_box_pos(3) - (cbar_pos(1) + cbar_pos(3));
    tightinset(i,4) = bounding_box_pos(4) - (cbar_pos(2) + cbar_pos(4));
end

% Delete the phantom axis.
delete(phantom_axis_h);
end
%% Function to find the axis handle associated to a given colorbar handle.
function peer_axis_h = get_colorbar_peer_axis(colorbar_handle)
% Find all axes handle in the parent figure
axes_handles = findall(colorbar_handle.Parent,'type','axes');
for i=1:numel(axes_handles)
    if axes_handles(i).Colorbar==colorbar_handle
        peer_axis_h = axes_handles(i);
        break
    end
end
end