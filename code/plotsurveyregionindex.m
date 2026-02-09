function h = plotsurveyregionindex(Tbl, var, opt)
%PLOTSURVEYREGIONINDEX Plot survey region index timeseries
%
% h = plotsurveyregionindex(Tbl, var, ...)
%
% This function plots 1D timeseries data.
%
% Input variables:
%
%   Tbl:    index data structure of tables (see readindexdata.m)
%
%   var:    name of variable to plot
%
% Optional input variables (passed as parameter/value pairs, default in []):
%
%   style:  plot style, corresponding to one of the following:
%
%           doy_vs_val_shaded: ESR-style plots of data over time versus day
%           of year, showing... 
%           - shaded percentiles of hindcast across years (in progressively
%             darker 10%-bin intervals)
%           - climatological hindcast median value in black solid line
%           - climatological hindcast mean in black dashed line
%           - most recent year of hindcast in red
%           - most recent persistence forecast in red dotted line
%           - black text labels of hindcast climatology value at start of
%             each month
%           - red text labels of hindcast most recent year value at start
%             of each month
%
%           other options coming (PEEC slides, etc.)
%           ['doy_vs_val_shaded']
%
%   col:    data sub-variable column to use for plots [1]
%
%   year:   year to highlight (not fully developed yet, currently mostly
%           relies on what data is in the Tbl input and highlights last
%           year)
%           [max(year(Tbl.Hc.t))]

% Copyright 2025 Kelly Kearney

arguments
    Tbl
    var {mustBeTextScalar}
    opt.style {mustBeTextScalar} ="doy_vs_val_shaded"
    opt.col (1,1) {mustBeInteger} =1
    opt.year (1,1) {mustBeInteger} =max(year(Tbl.Hc.t))
end

switch opt.style

    case "doy_vs_val_shaded" % ESR-style

        y = Tbl.Hc.(var)(:,opt.col);

        [yg, yy, tmid] = reshapetimeseries(Tbl.Hc.t, y, 'bin', 'date');

        h = plotgrid('size', [1 1], 'mar', 0.08, 'mr', 0.01);
        h.fig.Position(3:4) = [6 4]*72;
        h.fig.Color = 'w';

        [A, h.e] = ensemble2bnd(doy(tmid), yg, ...
            'dims', 'xey', ...
            'center', 'median', ...
            'prc', [0:10:40 60:10:100], ...
            'plot', 'boundedline', ...
            ...
            'alpha', true, ...
            'cmap', [0 0 0]);
        set(h.e.patch, 'facealpha', 0.1)
        hold on;
        h.yr = plot(doy(tmid), yg(:,end), 'r');
        h.clim = plot(doy(Tbl.Clim.t), Tbl.Clim.(var)(:,opt.col), '--k');
        h.fc = plot(doy(Tbl.Fc.t), Tbl.Fc.(var)(:,opt.col), ':r');
        
        tk = datetime(opt.year,1:12,1);
        tk.Format = 'MMM';
        set(h.ax, 'xtick', doy(tk), 'xticklabel', string(tk), 'xgrid', 'on', 'xlim', [0 365]);

        dypx = (10./(h.ax.Position(4)*h.fig.Position(4))).*diff(h.ax.YLim);

        tprop = {'fontname', 'andale mono', 'vert', 'top', 'fontsize', 8};
        for im = 1:12
            if im == 1
                tannotate = datetime(opt.year,im,2);
            else
                tannotate = datetime(opt.year,im,1);
            end

            hcval = interp1(doy(tmid), yg, doy(tannotate));
            [~,isrt] = sort(hcval);
            
            climval = interp1(doy(Tbl.Clim.t), Tbl.Clim.(var)(:,opt.col), doy(tannotate));

            ytext1 = h.ax.YLim(2) - dypx;
            ytext2 = h.ax.YLim(2) - dypx*2;
            if ~isnan(hcval(end))
                text(doy(tannotate), ytext1, sprintf('%0.2f', hcval(end)), 'color', 'r', tprop{:});
            end
            text(doy(tannotate), ytext2, sprintf('%0.2f', climval), 'color', 'k', tprop{:});

        end

end
