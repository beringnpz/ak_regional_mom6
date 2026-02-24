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
%
%   ax:     handle to target axis [gca]
%
% Output variables:
%
%   h:      structure of graphics object handles
%           doy_vs_val_shaded:
%               e.ln:    line objects, median (boundedline line)
%               e.patch: patch objects, percentiles (boundedline patch)
%               yr:      line object, target year
%               clim:    line object, climatological mean
%               fc:      line object, forecast

% Copyright 2025 Kelly Kearney

arguments
    Tbl
    var {mustBeTextScalar}
    opt.style {mustBeTextScalar} ="doy_vs_val_shaded"
    opt.col (1,1) {mustBeInteger} =1
    opt.year (1,1) {mustBeInteger} =max(year(Tbl.Hc.t))
    opt.ax =gca
end

if ~(year(Tbl.Fc.t(1)) == opt.year)
    warning('Forecast data year does not match target year... plots may not make sense!');
end

switch opt.style

    case "doy_vs_val_shaded" % ESR-style

        % Reshape to day-of year

        y = Tbl.Hc.(var)(:,opt.col);
        [yg, yy, tmid] = reshapetimeseries(Tbl.Hc.t, y, 'bin', 'date');
        iscurrent = yy == opt.year;

        % Plot shaded percentiles around median

        [A, h.e] = ensemble2bnd(doy(tmid), yg, ...
            'dims', 'xey', ...
            'center', 'median', ...
            'prc', [0:10:40 60:10:100], ...
            'plot', 'boundedline', ...
            'alpha', true, ...
            'cmap', [0 0 0]);
        set(h.e.patch, 'facealpha', 0.1)
        hold on;
        h.yr = plot(doy(tmid), yg(:,iscurrent), 'r');
        h.clim = plot(doy(Tbl.Clim.t), Tbl.Clim.(var)(:,opt.col), '--k');
        h.fc = plot(doy(Tbl.Fc.t), Tbl.Fc.(var)(:,opt.col), ':r');
        
        tk = datetime(opt.year,1:12,1);
        tk.Format = 'MMM';
        set(opt.ax, 'xtick', doy(tk), 'xticklabel', string(tk), 'xgrid', 'on', 'xlim', [0 365]);

        hfig = ancestor(opt.ax, 'figure'); % this does assume default axis (normalized) and figure (pixels) units
        dypx = (10./(opt.ax.Position(4)*hfig.Position(4))).*diff(opt.ax.YLim);

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

            ytext1 = opt.ax.YLim(2) - dypx;
            ytext2 = opt.ax.YLim(2) - dypx*2;
            if ~isnan(hcval(iscurrent))
                text(doy(tannotate), ytext1, sprintf('%0.2f', hcval(iscurrent)), 'color', 'r', tprop{:});
            end
            text(doy(tannotate), ytext2, sprintf('%0.2f', climval), 'color', 'k', tprop{:});

        end

    case "doy_vs_val_lines" % PEEC slide style

        y = Tbl.Hc.(var)(:,opt.col);
        [yg, yy, tmid] = reshapetimeseries(Tbl.Hc.t, y, 'bin', 'date');

        tk = datetime(opt.year,1:12,1);

        isavg = yy >= 1993 & yy < 2023; % baseline period for avg 1993-2022
        iscurrent = yy == opt.year;
        isprev3 = ismember(yy, opt.year - (3:-1:1));

        cmap = rgb(["teal";"purple";"medium blue"]);
        cidx = mod(yy(isprev3), 3);
        cidx(cidx==0) = 3;

        axes(opt.ax);
        h.lndoy = plot(doy(tmid), yg, 'color', rgb('light gray'));
        set(h.lndoy(iscurrent), 'color', 'r');
        set(h.lndoy(isprev3), {'color'}, num2cell(cmap(cidx,:),2));

        hold(opt.ax, 'on');
        h.lndoymean = plot(doy(tmid), mean(yg(:,isavg),2), '--k');

        h.lndoyfc = plot(doy(Tbl.Fc.t), Tbl.Fc.(var)(:,opt.col), '--r');

        % set(opt.ax, 'box', 'off', 'tickdir', 'out');
        set(opt.ax, 'xlim', [0 365], 'xtick', doy(tk), 'xticklabel', datestr(tk,'m'));

        % linkaxes(h.ax, 'y');
        % set([h.lndaily h.lnsummer h.lndoy' h.lndoyfc h.lnfc h.lnfcsummer], 'linewidth', 1.5)

    case "date_vs_val_annual_and_daily" % PEEC slide style

        y = Tbl.Hc.(var)(:,opt.col);

        [tlim(1), tlim(2)] = bounds([Tbl.Hc.t; Tbl.Fc.t]);
        yrlim = year(tlim);
        
        tsummer = datetime(yrlim(1):yrlim(2),7,1);
        ysummer = interp1(Tbl.Hc.t, y, tsummer);

        yfcsummer = interp1(Tbl.Fc.t, Tbl.Fc.(var)(:,opt.col), tsummer);

        isavg = tsummer >= datetime(1993,1,1) & tsummer < datetime(2023,1,1); % baseline period for avg 1993-2022

        yprc = prctile(ysummer(isavg), [25 50 75]);
        yavg = mean(ysummer(isavg));
        ystd = std(ysummer(isavg));

        h.lndaily = plot(opt.ax, Tbl.Hc.t, y, 'color', rgb('light gray'));
        hold(opt.ax, 'on');
        h.lnsummer = plot(opt.ax, tsummer, ysummer, ...
            'marker', 'o', 'color', rgb('green'), 'markerfacecolor', 'w');

        h.yavg = plot(opt.ax, tlim, ones(2,1).*yavg, 'color', rgb('gray'));
        % h.yprc = plot(h.ax(1), minmax(t), ones(2,1).*yprc, 'color', rgb('gray'), 'linestyle', ':');
        h.ystd = plot(opt.ax, tlim, ones(2,1).*(yavg + [-1 1].*ystd), 'color', rgb('gray'), 'linestyle', '--');


        h.lnfc = plot(opt.ax, Tbl.Fc.t, Tbl.Fc.(var)(:,opt.col), 'r');
        h.lnfcsummer = plot(opt.ax, tsummer, yfcsummer, ...
            'marker', 'o', 'color', 'r', 'markerfacecolor', 'w');


end
