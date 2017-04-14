function wilson_cowan_model_GUI()
% wilson_cowan_model_GUI.m
% MS 2016.10.06

% Wilson-Cowan model graphic user interface for visually testing different
% parameter values (synaptic weights, inputs, firing rate functions)

close all

% Parameters to be fed into Wilson-Cowan model
T = [0 600];
dt = 0.01;
c1_ee = 16; % Default for Gaussian firing rate (Meijer et al., 2015)
c2_ei = 18;
c3_ie = 12;
c4_ii = 3;
P = 4.00;
Q = 0.00;
inact_E = 0.25; % Inactiviation rate (tau_inact) of excitatory population
inact_I = 0.05; % Inactiviation rate (tau_inact) of inhibitory population
alpha_E = 2*1.5828; % Meijer et al. (2015) for sigmoid firing rate
alpha_I = 2*2.2201;
theta_E = 5.2516 * 1.1;
theta_I = 3.7512 * 0.9;

% Initialize parameters for loadng a series of synaptic weights
synT = [];
synW = [];
cHandles = [];
cMarkers = [];

% Set up figure, panels, and axes
h = 400;
w = 1200;
f = figure('Visible','on','Position',[20 700-h w h]);
f_panel = uipanel('Parent',f,'Position',[0.01 0.01 0.18 0.98],...
    'Title','Control Panel');
f_trace = uipanel('Parent',f,'Position',[0.20 0.01 0.26 0.98],...
    'Title','Time Series');
f_phase = uipanel('Parent',f,'Position',[0.465 0.01 0.26 0.98],...
    'Title','Phase Plane');
f_other = uipanel('Parent',f,'Position',[0.73 0.01 0.26 0.98],...
    'Title','Other');
f_trace_ax = axes('Parent',f_trace,'FontSize',8);
f_phase_ax = axes('Parent',f_phase,'FontSize',8);
f_input_ax = axes('Parent',f_other,'FontSize',8,...
    'Position',[.15 .55 .80 .35]);
f_firingrate_ax = axes('Parent',f_other,'FontSize',8,...
    'Position',[.15 .12 .80 .35]);

f_weights_ax = axes('Parent',f_panel,'FontSize',4,...
    'Position',[.04 .14 .92 .30]);

% Set up GUI buttons
f_panel_T = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.15 .93 .3 .05],...
    'Value',T(2),'String',num2str(T(2)),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.05 .93 .1 .05],...
    'String','T=')
f_panel_dt = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.65 .93 .3 .05],...
    'Value',dt,'String',num2str(dt),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.55 .93 .1 .05],...
    'String','dt=')
f_panel_c1 = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.65 .83 .3 .05],...
    'Value',c1_ee,'String',num2str(c1_ee),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.05 .83 .45 .05],...
    'String','c1 ( w_{E <- E} ) =')
f_panel_c2 = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.65 .77 .3 .05],...
    'Value',c2_ei,'String',num2str(c2_ei),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.05 .77 .45 .05],...
    'String','c2 ( w_{E <- I} ) =')
f_panel_c3 = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.65 .71 .3 .05],...
    'Value',c3_ie,'String',num2str(c3_ie),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.05 .71 .45 .05],...
    'String','c3 ( w_{I <- E} ) =')
f_panel_c4 = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.65 .65 .3 .05],...
    'Value',c4_ii,'String',num2str(c4_ii),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.05 .65 .45 .05],...
    'String','c4 ( w_{I <- I} ) =')
f_panel_P = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.15 .59 .30 .05],...
    'Value',P,'String',num2str(P),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.05 .59 .10 .05],...
    'String','P =')
f_panel_Q = uicontrol('Parent',f_panel,'Style','edit',...
    'Units','Normalized','Position',[.65 .59 .30 .05],...
    'Value',Q,'String',num2str(Q),...
    'Callback',@edit_callback);
uicontrol('Parent',f_panel,'Style','text',...
    'Units','Normalized','Position',[.55 .59 .10 .05],...
    'String','Q =')
f_panel_run = uicontrol('Parent',f_panel,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .02 .96 .09],...
    'String','Run Wilson-Cowan Model','Callback',@pushbutton_callback);

f_panel_load = uicontrol('Parent',f_panel,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .45 .96 .06],...
    'String','Load synaptic weight set','Callback',@pushbutton_callback);

% GUI callback functions

    function pushbutton_callback(source,~)
        switch source
            case f_panel_run
                weights.ee = c1_ee;
                weights.ei = c2_ei;
                weights.ie = c3_ie;
                weights.ii = c4_ii;
                params.T = T;
                params.dt = dt;
                params.P = P;
                params.Q = Q;
                params.alpha_E = alpha_E;
                params.alpha_I = alpha_I;
                params.theta_E = theta_E;
                params.theta_I = theta_I;
                params.inact_E = inact_E;
                params.inact_I = inact_I;
                
                set(f, 'HandleVisibility', 'off');
                close all;
                set(f, 'HandleVisibility', 'on');
                [t, E, I, J_E, J_I, F_E, F_I] = get_WC_deterministic_for_GUI(...
                                                weights, params);
                update_plots(t,E,I,J_E,J_I,F_E,F_I)
                get_WC_phase_plane(weights, params)
                
            case f_panel_load
                [FileName,PathName,FilterIndex] = uigetfile('.mat');
                if FilterIndex == 1
                    S = load([PathName,FileName]);
                    synT = S.sampleTimes;
                    synW = [S.W_EE, S.W_EI, S.W_IE, S.W_II];
                    
                    cHandles = gobjects(length(synT),4);
                    hold(f_weights_ax,'off')
                    plot(f_weights_ax,synT,synW(:,1),'r')
                    hold(f_weights_ax,'on')
                    plot(f_weights_ax,synT,synW(:,2),'m')
                    plot(f_weights_ax,synT,synW(:,3),'g')
                    plot(f_weights_ax,synT,synW(:,4),'b')
                    for i = 1:length(synT)
                        cHandles(i,1) = plot(f_weights_ax,synT(i),synW(i,1),'r.','MarkerSize',12);
                        cHandles(i,2) = plot(f_weights_ax,synT(i),synW(i,2),'m.','MarkerSize',12);
                        cHandles(i,3) = plot(f_weights_ax,synT(i),synW(i,3),'g.','MarkerSize',12);
                        cHandles(i,4) = plot(f_weights_ax,synT(i),synW(i,4),'b.','MarkerSize',12);
                    end
                    set(cHandles,'ButtonDownFcn',@weight_select_callback)
                end           
        end
    end

    function edit_callback(source,~)
        val = str2double(get(source,'String'));
        if isempty(val)
            warning('numerical input required')
            return
        end
        switch source
            case f_panel_T
                T(2) = val;
            case f_panel_dt
                dt = val;
            case f_panel_c1
                c1_ee = val;
            case f_panel_c2
                c2_ei = val;
            case f_panel_c3
                c3_ie = val;
            case f_panel_c4
                c4_ii = val;
            case f_panel_P
                P = val;
            case f_panel_Q
                Q = val;
        end
    end

    function update_plots(t,E,I,J_E,J_I,F_E,F_I)
        plot(f_trace_ax,t,E,'r');
        hold(f_trace_ax,'on');
        plot(f_trace_ax,t,I,'b');
        hold(f_trace_ax,'off');
        legend(f_trace_ax,'E','I')
        
        plot(f_phase_ax,E,I)
        set(f_phase_ax,'xlim',[0 0.5],'ylim',[0 0.5])
        
        plot(f_input_ax,t,J_E,'r');
        hold(f_input_ax,'on');
        plot(f_input_ax,t,J_I,'b');
        hold(f_input_ax,'off');
        legend(f_input_ax,'J_E','J_I')
        
        plot(f_firingrate_ax,t,F_E,'r');
        hold(f_firingrate_ax,'on');
        plot(f_firingrate_ax,t,F_I,'b');
        hold(f_firingrate_ax,'off');
        legend(f_firingrate_ax,'F_E','F_I')
        
        f_trace_ax.XLabel.String = 'time';
        f_trace_ax.YLabel.String = 'activity';
        f_trace_ax.YLim = [0 1];
        f_phase_ax.XLabel.String = 'E activity';
        f_phase_ax.YLabel.String = 'I activity';
        f_phase_ax.XLim = [0 1];
        f_phase_ax.YLim = [0 1];
        
        f_firingrate_ax.XLabel.String = 'time';
        f_firingrate_ax.YLabel.String = 'firing rate';
        f_input_ax.YLabel.String = 'input current';
    end
end