import os
import numpy as np
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
from matplotlib.ticker import MultipleLocator

def bar_plot (data,
              error,
              xlabels = None,
              ylabels = None,
              xlabel = None,
              ylabel = None,
              barwidth = 0.75,
              colors = 'turquoise',
              capsize = 0,
              fmt = 'none',
              msize = 12,
              ewidth = 1.5,
              edgecolor = 'none',
              ecolors = 'k',
              fontsize = 12,
              xlim = None,
              ylim = None,
              yMinorTicks = False,
              adjustBottom = False,
              shiftBottomAxis = None,
              xbounds = None,
              leg = None,
              show = True,
              figdir = None,
              figname = None):
    """Bar plot with one group.

    Args:
        data (list): data of one group to plot.
        error (list): error on each bar in one group, empty list for no error bars.
        xlabels (list): label for each bar on the x-axis.
        ylabels (list(numeric)): label for each tick on the y-axis.
        xlabel (str): label on the x-axis.
        ylabel (str): label on the y-axis.
        barwidth (numeric): width of each bar.
        colors (list): list of colors for each bar, or one color for all bars.
        capsize (numeric): width of upper and lower caps on error bars.
        fmt (str): error bar center point format.
        msize (numeric): size for error marker.
        ewidth (numeric): thickness for error bars, including caps. 
        edgecolor (str): color for bar edges.
        ecolors (list):  list of colors for each error bar, or one color for all bars.
        fontsize (numeric): font size for x and y tick labels and axis labels.
        xlim (list): limits for x-axis.
        ylim (list): limits for y-axis.
        yMinorTicks (boolean): True to show minor ticks on y-axis.
        adjustBottom (boolean) = True if bottom margin of figure needs adjustment.
        shiftBottomAxis (numeric) = shift x-axis by this much.
        xbounds (list): chop off x-axis beyond these bounds.
        leg (list): legend labels for two groups.
        show (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()
    ph = plt_fig.add_subplot(111)
    
    if adjustBottom:
        plt_fig.subplots_adjust(bottom = adjustBottom)
    if shiftBottomAxis is not None:
        ph.spines["bottom"].set_position(("axes", shiftBottomAxis))
    if xbounds is not None:
        ph.spines['bottom'].set_bounds(xbounds[0], xbounds[1])
    if xlim is not None:
        ph.set_xlim(xlim)
    if ylim is not None:
        ph.set_ylim(ylim)
    
    numBars = len(data)
    if isinstance(colors, list):
        if len(colors) < numBars:
            colors.extend(['b'] * (numBars - len(colors)))
    ind = np.array(range(numBars)) + 1
    ph.bar(ind,
           data,
           barwidth,
           color = colors,
           edgecolor = edgecolor)
    
    if len(error) > 0:
        if len(error) < numBars:
            error.extend([0] * (numBars - len(error)))
        if isinstance(ecolors, str):
            ecolors = [ ecolors ] * numBars
        elif len(ecolors) < numBars:
            ecolors.extend(['k'] * (numBars - len(ecolors)))
        for pos, d, err, ecolor in zip(ind, data, error, ecolors):
            if type(err) in (list, tuple):
                err = [ [e] for e in err ]
            ph.errorbar(pos,
                        d,
                        yerr = err,
                        elinewidth = ewidth,
                        fmt = fmt,
                        markersize = msize,
                        capsize = capsize,
                        capthick = ewidth,
                        ecolor = ecolor)
    
    ph.tick_params('both', length=10, which='major')
    ph.tick_params('both', length=5, which='minor')
    ph.get_xaxis().set_tick_params(which='both', bottom=False, top=False, labelbottom=False)
    if xlabels is not None:
        ph.set_xticks(ind)
        ph.xaxis.set_ticks_position('bottom')
        ph.get_xaxis().set_tick_params(which='both', direction='out')
        ph.set_xticklabels(xlabels)
    if xlabel is not None:
        ph.set_xlabel(xlabel)
    ph.get_yaxis().set_tick_params(which='both', right=False, direction='out')
    ph.yaxis.set_ticks_position('left')
    if ylabels is not None:
        ph.set_yticks(ylabels)
        if all(n % 1 == 0 for n in ylabels):
            ylabels = [int(n) for n in ylabels]
        ph.set_yticklabels(ylabels)
    if yMinorTicks:
        yticks = ph.get_yticks(minor=False)
        minorLocator = MultipleLocator( (yticks[1] - yticks[0]) / 5)
        ph.yaxis.set_minor_locator(minorLocator)
    if ylabel is not None:
        ph.set_ylabel(ylabel)
    ph.spines["top"].set_visible(False)
    ph.spines['right'].set_visible(False)
    ph.spines['bottom'].set_position('zero')
    for item in ([ph.xaxis.label, ph.yaxis.label] + ph.get_xticklabels() + ph.get_yticklabels()):
        item.set_fontsize(fontsize)
    if leg is not None:
        ph.legend(leg, bbox_to_anchor = (0.5,1.02), loc = 'lower center', ncol=2, frameon=False)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def multi_bar_plot (data,
                    errors,
                    xlabels = None,
                    ylabels = None,
                    xlabel = None,
                    ylabel = None,
                    colors = 'b',
                    barwidth = 0.3,
                    capsize = 0,
                    fmt = 'none',
                    msize = 12,
                    ewidth = 1.5,
                    edgecolor = 'none',
                    ecolors = 'k',
                    fontsize = 12,
                    xlim = None,
                    ylim = None,
                    adjustBottom = False,
                    shiftBottomAxis = None,
                    xbounds = None,
                    leg = None,
                    show = True,
                    figdir = None,
                    figname = None):
    """Bar plot with multiple groups.

    Args:
        data (list): list of data for each group to plot.
        errors (list): list of errors for each group, empty lists for no error bars.
        xlabels (list): label for each bar on the x-axis.
        ylabels (list(numeric)): label for each tick on the y-axis.
        xlabel (str): label on the x-axis.
        ylabel (str): label on the y-axis.
        colors (list): list of colors for each group, optionally one color for all bars per group.
        barwidth (numeric): width of each bar.
        capsize (numeric): width of upper and lower caps on error bars.
        fmt (str): error bar center point format.
        msize (numeric): size for error marker.
        ewidth (numeric): thickness for error bars, including caps. 
        edgecolor (str): color for bar edges.
        ecolors (list):  list of error bar colors for rach group, optionally one color for all bars per group.
        fontsize (numeric): font size for x and y tick labels and axis labels.
        xlim (list): limits for x-axis.
        ylim (list): limits for y-axis.
        adjustBottom (boolean) = True if bottom margin of figure needs adjustment.
        shiftBottomAxis (numeric) = shift x-axis by this much.
        xbounds (list): chop off x-axis beyond these bounds.
        leg (list): legend labels for two groups.
        show (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()
    ph = plt_fig.add_subplot(111)
    
    if adjustBottom:
        plt_fig.subplots_adjust(bottom = adjustBottom)
    if shiftBottomAxis is not None:
        ph.spines["bottom"].set_position(("axes", shiftBottomAxis))
    if xbounds is not None:
        ph.spines['bottom'].set_bounds(xbounds[0], xbounds[1])
    
    numGroups = len(data)
    if len(errors) < numGroups:
        errors.extend([[]] * (numGroups - len(errors)))
    if leg is None:
        leg = ['Group ' + str(i+1) for i in range(numGroups)]
    if isinstance(colors, str):
        colors = [colors] * numGroups
    elif len(colors) < numGroups:
        colors.extend(['b'] * (numGroups - len(colors)))
    if isinstance(ecolors, str):
        ecolors = [ecolors] * numGroups
    elif len(ecolors) < numGroups:
        ecolors.extend(['k'] * (numGroups - len(ecolors)))
    if xlim is not None:
        ph.set_xlim(xlim)
    if ylim is not None:
        ph.set_ylim(ylim)
    
    for i, d, error, color, ecolor, label in zip(range(numGroups), data, errors, colors, ecolors, leg):
        numBars = len(d)
        ind = np.array(range(numBars)) + (i * barwidth) + 1
        ph.bar(ind,
               d,
               barwidth,
               color = color,
               edgecolor = edgecolor,
               label = label)
        
        if len(error) > 0:
            if len(error) < numBars:
                error.extend([0] * (numBars - len(error)))
            if isinstance(ecolor, str):
                ecolor = [ecolor] * numBars
            elif len(ecolor) < numBars:
                ecolor.extend(['k'] * (numBars - len(ecolor)))
            for pos, d, err, ecol in zip(ind, d, error, ecolor):
                ph.errorbar(pos,
                            d,
                            err,
                            elinewidth = ewidth,
                            fmt = fmt,
                            markersize = msize,
                            capsize = capsize,
                            capthick = ewidth,
                            ecolor = ecol)
    
    ind = np.array(range(numBars)) + ((numGroups - 1) * barwidth / 2.) + 1
    ph.tick_params('both', length=10, which='major')
    ph.get_xaxis().set_tick_params(which='both', bottom=False, top=False, labelbottom=False)
    if xlabels is not None:
        ph.set_xticks(ind)
        ph.xaxis.set_ticks_position('bottom')
        ph.get_xaxis().set_tick_params(which='both', direction='out')
        ph.set_xticklabels(xlabels)
    ph.get_yaxis().set_tick_params(which='both', right=False, direction='out')
    ph.yaxis.set_ticks_position('left')
    if ylabels is not None:
        ph.set_yticks(ylabels)
        ph.set_yticklabels(ylabels)
    if xlabel is not None:
        ph.set_xlabel(xlabel)
    if ylabel is not None:
        ph.set_ylabel(ylabel)
    ph.spines["top"].set_visible(False)
    ph.spines['right'].set_visible(False)
    ph.spines['bottom'].set_position('zero')
    for item in ([ph.xaxis.label, ph.yaxis.label] + ph.get_xticklabels() + ph.get_yticklabels()):
        item.set_fontsize(fontsize)
    if leg is not None:
        ph.legend(bbox_to_anchor = (0.5,1.02), loc = 'lower center', ncol=2, frameon=False)
        #ph.legend(bbox_to_anchor = (2,0), loc = 'lower center', ncol=2, frameon=False)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def curve_plot (xdata,
                ydata,
                xrange = None,
                yrange = None,
                styles = None,
                msize = 12,
                mwidth = 1,
                xlabel = None,
                ylabel = None,
                xticklabels = None,
                fontsize = 12,
                leg = None,
                show = True,
                figdir = None,
                figname = None):
    """Plot multiple curves.

    Args:
        xdata (list): x-axis data points.
        ydata (list): y-axis data points, list for each curve.
        xrange (list): limits on x-axis.
        yrange (list): limits on y-axis.
        styles (list): plot style for each curve.
        msize (numeric): marker size.
        mwidth (numeric): marker edge width.
        xlabel (str): label for x-axis.
        ylabel (str): label for y-axis.
        xticklabels (list): label for each tick on the x-axis.
        fontsize (numeric): font size for x and y tick labels and axis labels.
        leg (list): label for each group.
        show (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()
    ph = plt_fig.add_subplot(111)
    if styles is None:
        styles = ['.'] * len(ydata)
    for i, sample in enumerate(ydata):
        if len(xdata) > 0:
            ph.plot(xdata[i], sample, styles[i], markersize=msize, markeredgewidth=mwidth)
        else:
            ph.plot(sample, styles[i])
    ph.spines["top"].set_visible(False)
    ph.spines['right'].set_visible(False)
    ph.xaxis.set_ticks_position('bottom')
    ph.yaxis.set_ticks_position('left')
    ph.get_xaxis().set_tick_params(which='both', direction='out')
    ph.get_yaxis().set_tick_params(which='both', direction='out')
    if xticklabels is not None:
        ph.set_xticklabels(xticklabels)
    if xrange is not None:
        plt.xlim(xrange)
    if yrange is not None:
        plt.ylim(yrange)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if leg is not None:
        plt.legend(leg)
    for item in ([ph.xaxis.label, ph.yaxis.label] + ph.get_xticklabels() + ph.get_yticklabels()):
        item.set_fontsize(fontsize)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def box_plot (data,
              xlabels,
              xlabel = None,
              ylabel = None,
              fontsize = 12,
              xlim = None,
              ylim = None,
              adjustBottom = False,
              shiftBottomAxis = None,
              xbounds = None,
              ybounds = None,
              colors = 'steelblue',
              show = True,
              figdir = None,
              figname = None):
    """Box plot.

    Args:
        data (list): data to plot, each group is a list.
        xlabels (list): label for each box on the x-axis.
        xlabel (str): label for x-axis.
        ylabel (str): label for y-axis.
        fontsize (numeric): font size for x and y tick labels and axis labels.
        xlim (list): limits for x-axis.
        ylim (list): limits for y-axis.
        adjustBottom (boolean) = True if bottom margin of figure needs adjustment.
        shiftBottomAxis (numeric) = shift x-axis by this much.
        xbounds (list): chop off x-axis beyond these bounds.
        ybounds (list): chop off y-axis beyond these bounds.
        colors (list): color for each box, optionally one color for all boxes.
        show (boolean): 'yes' to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()
    ph = plt_fig.add_subplot(111)
    bp = ph.boxplot(data, patch_artist=True)
    
    numBoxes = len(data)
    if isinstance(colors, str):
        colors = [colors] * numBoxes
    elif len(colors) < numBoxes:
        colors.extend(['b'] * (numBoxs - len(colors)))
    
    ## change outline color, fill color and linewidth of the boxes
    for box, color in zip(bp['boxes'], colors):
        # change outline color
        box.set(color='black', linewidth=1)
        # change fill color
        box.set(facecolor=color)

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='black', linestyle='-', linewidth=1)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='black', linewidth=1)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='black', linewidth=1)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', markerfacecolor = 'white', alpha=0.5)
    
    if adjustBottom:
        plt_fig.subplots_adjust(bottom = adjustBottom)
    if shiftBottomAxis is not None:
        ph.spines["bottom"].set_position(("axes", shiftBottomAxis))
    if xbounds is not None:
        ph.spines['bottom'].set_bounds(xbounds[0], xbounds[1])
    if ybounds is not None:
        ph.spines['left'].set_bounds(ybounds[0], ybounds[1])
    if xlim is not None:
        ph.set_xlim(xlim)
    if ylim is not None:
        ph.set_ylim(ylim)
    if xlabel is not None:
        ph.set_xlabel(xlabel)
    if ylabel is not None:
        ph.set_ylabel(ylabel)
    ph.tick_params('both', length=10, which='major')
    ph.set_xticklabels(xlabels)
    ph.get_xaxis().set_tick_params(which='both', direction='out')
    ph.get_yaxis().set_tick_params(which='both', direction='out')
    ph.spines["top"].set_visible(False)
    ph.spines['right'].set_visible(False)
    ph.xaxis.set_ticks_position('bottom')
    ph.yaxis.set_ticks_position('left')
    for item in ([ph.xaxis.label, ph.yaxis.label] + ph.get_xticklabels() + ph.get_yticklabels()):
        item.set_fontsize(fontsize)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def multi_histogram_plot (samples,
                          colors,
                          xlabel = None,
                          ylabel = None,
                          leg = None,
                          edgecolor = 'none',
                          fontsize = 12,
                          bins = 10,
                          alpha = 0.5,
                          xlim = None,
                          ylim = None,
                          show = True,
                          figdir = None,
                          figname = None):
    """Plot multiple histograms in the same figure.

    Args:
        samples (list): data samples to plot.
        colors (list): color for each histogram.
        xlabel (str): label for x-axis.
        ylabel (str): label for y-axis.
        xlabels (list): label for each bar on the x-axis.
        ylabels (list(numeric)): label for each tick on the y-axis.
        leg (list): label for each histogram.
        edgecolor (str): color for bar edges.
        fontsize (numeric): font size for x and y tick labels and axis labels.
        bins (int): number of bins per histogram.
        alpha (numeric): 0.0 transparent through 1.0 opaque.
        xlim (list): limits for x-axis.
        ylim (list): limits for y-axis.
        show (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()
    ph = plt_fig.add_subplot(111)
    for i, sample in enumerate(samples):
        if leg is None:
            ph.hist(sample, color=colors[i], bins=bins, alpha=alpha, edgecolor=edgecolor)
        else:
            ph.hist(sample, color=colors[i], bins=bins, alpha=alpha, label=leg[i], edgecolor=edgecolor)
    if xlabel is not None:
        ph.set_xlabel(xlabel)
    if ylabel is not None:
        ph.set_ylabel(ylabel)
    if leg is not None:
        ph.legend(loc="upper right")
    if xlim is not None:
        ph.set_xlim(xlim)
    if ylim is not None:
        ph.set_ylim(ylim)
    ph.tick_params('both', length=10, which='major')
    ph.get_xaxis().set_tick_params(which='both', top=False, direction='out')
#     if xlabels is not None:
#         ph.set_xticks(ind)
#         ph.xaxis.set_ticks_position('bottom')
#         ph.get_xaxis().set_tick_params(which='both', direction='out')
#         ph.set_xticklabels(xlabels)
    ph.get_yaxis().set_tick_params(which='both', right=False, direction='out')
    ph.yaxis.set_ticks_position('left')
#     if ylabels is not None:
#         ph.set_yticks(ylabels)
#         ph.set_yticklabels(ylabels)
    if xlabel is not None:
        ph.set_xlabel(xlabel)
    if ylabel is not None:
        ph.set_ylabel(ylabel)
    ph.spines["top"].set_visible(False)
    ph.spines['right'].set_visible(False)
    ph.spines['bottom'].set_position('zero')
    for item in ([ph.xaxis.label, ph.yaxis.label] + ph.get_xticklabels() + ph.get_yticklabels()):
        item.set_fontsize(fontsize)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def heatmap_plot(mat,
                 cmap = 'Blues',
                 interpolation = 'nearest',
                 vmin = None,
                 vmax = None,
                 barPos = None,
                 orientation = 'vertical',
                 barTicks = None,
                 barLabels = None,
                 showBar = True,
                 show = True,
                 figdir = None,
                 figname = None):
    """Two-dimensional heatmap plot.

    Args:
        mat (array): two-dimensional array of data to plot.
        cmap (str): colormap for heatmap.
        interpolation (str): colormap interpolation.
        vmin (numeric): minimum value on color scale.
        vmax (numeric): maximum value on color scale.
        barPos (list): position for colorbar [left, bottom, width, height].
        orientation (str): colorbar orientation (vertical or horizontal).
        barTicks (list): tick positions on colorbar.
        barLabels (list): tick labels on colorbar.
        showBar (boolean): True to show colorbar, otherwise bar is not shown.
        show (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()
    ph = plt_fig.add_subplot(111)
    img = ph.imshow(mat,
                    vmin = vmin,
                    vmax = vmax,
                    cmap = cmap,
                    interpolation = interpolation)
    
    ph.spines["left"].set_visible(False)
    ph.spines['bottom'].set_visible(False)
    ph.get_xaxis().set_tick_params(which='both', bottom=False, top=False, labelbottom=False)
    ph.get_yaxis().set_tick_params(which='both', left=False, right=False, labelleft=False)
    
    if showBar:
        barAxis = None
        if barPos is not None:
            barAxis = plt_fig.add_axes(barPos)
        cbar = plt.colorbar(img,
                            cax = barAxis,
                            orientation = orientation,
                            ticks = barTicks)
        cbar.ax.tick_params(length = 0)
        if barLabels is not None:
            cbar.ax.set_yticklabels(barLabels)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def pie_plot (data,
              labels = None,
              labeldist = 1.1,
              angle = None,
              pct = None,
              pctdist = 0.6,
              colors = None,
              show = True,
              figdir = None,
              figname = None):
    """Pie plot.

    Args:
        data (list): data to plot.
        labels (list): label for each slice in the plot.
        labeldist (float): radial distance at which pie chart labels are drawn.
        angle (numeric): starting angle for first slice in the plot.
        pct (str): format for percentage labels.
        pctdist (numeric): distance of pct labels from center.
        colors (list): colors for pie slices.
        show (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()
    ph = plt_fig.add_subplot(111)
    pie_return = ph.pie(data,
                        labels = labels,
                        labeldistance = labeldist,
                        startangle = angle,
                        colors = colors,
                        autopct = pct,
                        pctdistance = pctdist)
    ph.axis('equal')
    for w in pie_return[0]:
        w.set_linewidth(0)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def venn2_plot (data,
                labels = None,
                colors = None,
                hideNum = False,
                show = True,
                figdir = None,
                figname = None):
    """Venn diagram for two groups.

    Args:
        data (list): groups to plot, each group is a list.
        labels (list): label for each group.
        colors (list): color for each group.
        hideNum (boolean): True to hide numbers on venn diagram.
        show (boolean): 'yes' to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()           
    v = venn2(data, set_labels = labels)
    if colors is not None:
        v.get_patch_by_id('10').set_color(colors[0])
        v.get_patch_by_id('01').set_color(colors[1])
    if hideNum:
        v.get_label_by_id('10').set_text('')
        v.get_label_by_id('01').set_text('')
        if v.get_label_by_id('11') is not None:
            v.get_label_by_id('11').set_text('')
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def venn3_plot (data,
                labels = None,
                colors = None,
                hideNum = False,
                show = True,
                figdir = None,
                figname = None):
    """Venn diagram for two groups.

    Args:
        data (list): groups to plot, each group is a list.
        labels (list): label for each group.
        colors (list): color for each group.
        hideNum (boolean): True to hide numbers on venn diagram.
        show (boolean): 'yes' to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    plt_fig = plt.figure()           
    v = venn3(data, set_labels = labels)
    if colors is not None:
        v.get_patch_by_id('100').set_color(colors[0])
        v.get_patch_by_id('010').set_color(colors[1])
        v.get_patch_by_id('001').set_color(colors[2])
    if hideNum:
        v.get_label_by_id('100').set_text('')
        v.get_label_by_id('010').set_text('')
        v.get_label_by_id('001').set_text('')
        if v.get_label_by_id('110') is not None:
            v.get_label_by_id('110').set_text('')
        if v.get_label_by_id('011') is not None:
            v.get_label_by_id('011').set_text('')
        if v.get_label_by_id('101') is not None:
            v.get_label_by_id('101').set_text('')
        if v.get_label_by_id('111') is not None:
            v.get_label_by_id('111').set_text('')
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)

def network_plot (edges,
                  nodes = None,
                  nodeSizes = None,
                  edgeWidth = 1,
                  nodeColors = None,
                  edgeColors = None,
                  node_cmap = None,
                  node_vmin = None,
                  node_vmax = None,
                  edge_cmap = None,
                  edge_vmin = None,
                  edge_vmax = None,
                  layout = 'graphviz',
                  show = True,
                  figdir = None,
                  figname = None):
    """Network plot.

    Args:
        edges (list): network edges in the form of node tuples.
        nodes (list): network nodes.
        nodeSizes (numeric): node size in network plot.
        edgeWidth (numeric): edge width in network plot.
        nodeColors (list): color for each node in network plot.
        edgeColors (list): color for each edge in network plot.
        node_cmap (str or Colormap instance): colormap for nodes, if node colors are numbers.
        node_vmin (numeric): minimum value on node color scale.
        node_vmax (numeric): maximum value on node color scale.
        edge_cmap (str or Colormap instance): colormap for edges, if edge colors are numbers.
        edge_vmin (numeric): minimum value on edge color scale.
        edge_vmax (numeric): maximum value on edge color scale.
        layout (str): network layout.
        show (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in.
    
    """
    g = nx.Graph()
    if nodes is None:
        nodes = set()
        for edge in edges:
            nodes.update(edge[:2])
    nodes = list(nodes)
    if nodeSizes is None:
        nodeSizes = [20] * len(nodes)
    if edgeColors is None:
        edgeColors = ['black'] * len(edges)
    if nodeColors is None:
        nodeColors = ['blue'] * len(nodes)
    for node, size, color in zip(nodes, nodeSizes, nodeColors):
        g.add_node(node, size=size, color=color)
    for (node1, node2,), color in zip(edges, edgeColors):
        g.add_edge(node1, node2, color=color)
    edges = g.edges()
    nodes = g.nodes()
    edgeColors = [g[u][v]['color'] for u,v in edges]
    nodeSizes = [g.node[u]['size'] for u in nodes]
    nodeColors = [g.node[u]['color'] for u in nodes]
    pos = {'graphviz': graphviz_layout(g),
           'spring': nx.spring_layout(g),
           'shell': nx.shell_layout(g),
           'circular': nx.circular_layout(g),
           'spectral': nx.spectral_layout(g),
           'random': nx.random_layout(g)}
    
    plt_fig = plt.figure()
    nx.draw(g,
            nodes = nodes,
            edges = edges,
            node_size = nodeSizes,
            width = edgeWidth,
            node_color = nodeColors,
            edge_color = edgeColors,
            pos = pos[layout],
            cmap = node_cmap,
            vmin = node_vmin,
            vmax = node_vmax,
            edge_cmap = edge_cmap,
            edge_vmin = edge_vmin,
            edge_vmax = edge_vmax)
    if (figdir is not None) and (figname is not None):
        figdir_eps = figdir / 'eps_figures'
        figdir_pdf = figdir / 'pdf_figures'
        figdir_png = figdir / 'png_figures'
        if not figdir_eps.exists():
            os.makedirs(figdir_eps)
        if not figdir_pdf.exists():
            os.makedirs(figdir_pdf)
        if not figdir_png.exists():
            os.makedirs(figdir_png)
        plt.savefig(os.path.join(figdir_eps, figname+".eps"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_pdf, figname+".pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(figdir_png, figname+".png"), bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(plt_fig)
