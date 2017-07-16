    function MS_waveform(ax,t,y,varargin)
        % Function plots and formats waveform on the provided axes
        
        % Convert time vector into seconds
        t = (t-t(1))*24*60*60;
        % Plot waveform and format axis
        plot(ax,t,y,'k')
        set(ax,'xlim',[t(1) t(end)])
        
        % Variable input arguments:
        % YLIM
        if ~isempty(varargin)
            if varargin{1}(1) < varargin{1}(2)
                set(ax,'ylim',varargin{1})
            end
        end
    end