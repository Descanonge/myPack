"""Various convenience functions related to plots."""

import numpy as np

from scipy.ndimage import label

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import mpl_toolkits.axes_grid1 as axes_grid1
import matplotlib.ticker as mt

import myPack.analysis as mpa
import myPack._plot_boxes as _plot_boxes


def LatFormatter(fmt):
    """Formatter for latitude values."""
    def func(x, pos):
        return mpa.latlon2str(lat=x, fmt='%lat', fmtF=fmt)
    return mt.FuncFormatter(func)


def LonFormatter(fmt):
    """Formatter for longitude values."""
    def func(x, pos):
        return mpa.latlon2str(lon=x, fmt='%lon', fmtF=fmt)
    return mt.FuncFormatter(func)


def add_cax(ax, loc="right", size=0.1, pad=0.05, **kwargs):
    """Add an ax for colorbar next to a existing axis.

    Wrapper around make_axes_locatable
    """
    divider = axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes(loc, size, pad, **kwargs)
    return cax


def make_imageGrid(fig, rect=111, **grid_kw):
    """Create a grid using axesgrid1.AxesGrid.

    Returns grid
    """
    return axes_grid1.AxesGrid(fig, rect=111, **grid_kw)


def add_box(fig, axes, idx1, idx2, hpad, vpad, offset=None, **kwargs):
    """Add box around axes.

    idx1: lower left axes index in the axes
    idx2: upper right axes index in the axes
    hpad/vpad: horizontal and vertical padding between the axes in points
    offset: additional offset [left, right, bottom, top]

    kwargs transmitted to the rectangle artist

    Returns xy (bottom left position) and wh (width, height)
    """

    if not offset:
        offset = [0., 0., 0., 0.]
    offset_f = [0.]*4

    hpad_f, vpad_f = _plot_boxes.convert_pad(fig, hpad, vpad)
    offset_f[:2] = _plot_boxes.convert_pad(fig, *offset[:2])
    offset_f[2:] = _plot_boxes.convert_pad(fig, *offset[2:])
    xy, wh = _plot_boxes.get_rect(idx1, idx2, axes, hpad_f, vpad_f, offset_f)
    r = plt.Rectangle(xy, *wh, transform=fig.transFigure, **kwargs)
    fig.add_artist(r)

    return xy, wh


def plot_line_mask(x, y, boxes):
    """Plot a line masking inside of boxes.

    By discretisation. Not that elegant but easy to write

    Everything in figure coordinates
    """

    N = 1000
    x = np.linspace(*x, N)
    y = np.linspace(*y, N)

    mask = np.zeros(N, bool)

    for box in boxes:
        mask ^= ((x > box[0][0]) * (x < box[1][0])
                 * (y > box[0][1]) * (y < box[1][1]))

    L, n = label(~mask)
    segx = []
    segy = []
    print(1*mask)
    for i in range(n):
        seg = np.where(L == i)[0]
        segx.append([x[seg[0]], x[seg[-1]]])
        segy.append([y[seg[0]], y[seg[-1]]])

    return segx, segy


def plot_line(ax, slope, intercept=None, length=None, **kwargs):
    """Plot a straight line extending over axes.

    Plot a straight line from a slope and intercept point.
    kwargs is passed to ax.plot

    Whatever the scale is (lin or log) plot a variant of Y=slope*X + intercept
    Intercept can be the y coordinate, or both coordinates of the point.
    If intercept is a float, the intercept point is placed at X=0.
        (so eventually log(x)=0)

    If length is specified, only a portion of the line is drawn
        (in axes fraction units)

    If lin-lin:
        If no intercept is supplied, intercept=(x=0, y=0)
        y = slope*x + intercept
    If log-log:
        If no intercept is supplied, intercept=(x=1, y=1)
        log10(y) = slope*log10(x) + log10(intercept)
        This corresponds to y = intercept * x**(slope)
    If semi-log (log in y):
        If no intercept is supplied, intercept=(x=0, y=1)
        log10(y) = slope*x + intercept
        This corresponds to y = intercept * 10**(slope*x)
    If semi-log (log in x):
        If no intercept is supplid, intercept=(y=0, x=1)
        y = slope*log10(x) + intercept

    All is done by going into axes coordinates to avoid reaching too big
    numbers when using log scales. There is thus a risk of uncorrelating data
    and display.
    """

    # Find axes aspect ratio and intercept depending on scales
    if intercept is None:
        intercept = [None, None]
    if isinstance(intercept, (int, float)):
        intercept = [None, intercept]

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    if ax.get_xscale() == 'linear':
        if intercept[0] is None:
            intercept[0] = 0
    elif ax.get_xscale() == 'log':
        if intercept[0] is None:
            intercept[0] = 1
        xlim = [np.log10(z) for z in xlim]
    else:
        raise ValueError("Unsupported type of axis scale")

    if ax.get_yscale() == 'linear':
        if intercept[1] is None:
            intercept[1] = 0
    elif ax.get_yscale() == 'log':
        if intercept[1] is None:
            intercept[1] = 1
        ylim = [np.log10(z) for z in ylim]
    else:
        raise ValueError("Unsupported type of axis scale")

    xscale = xlim[1] - xlim[0]
    yscale = ylim[1] - ylim[0]

    # Pass intercept and slope into axes coordinates
    slope_ = slope * xscale/yscale
    a = slope_

    tf = ax.transData + ax.transAxes.inverted()
    intercept_ = list(tf.transform(intercept))

    if length is not None:
        points_ = []
        points_.append(intercept_)
        x2 = intercept_[0] + length/np.sqrt(1+a**2)
        y2 = intercept_[1] + (x2-intercept_[0])*a
        points_.append([x2, y2])
    else:
        points_ = mpa.get_intersection_line(intercept_, slope_, [0, 0], [1, 1])

    # Go back to data coordinates
    if len(points_) > 0:
        points_ = np.array(list(points_))[:2]
        points_ = points_[np.argsort(points_[:, 0])]
        points = tf.inverted().transform(points_)

        ax.set_autoscale_on(False)
        ax.plot((points[0, 0], points[1, 0]),
                (points[0, 1], points[1, 1]), **kwargs)
        ax.set_autoscale_on(True)

    return points_


def plot_cursor(ax, xp, yp, **kwargs):
    """Plot lines going from (xp, yp) to view limits."""
    x = (ax.get_xlim()[0], xp, xp)
    y = (yp, yp, ax.get_ylim()[0])

    ax.set_autoscale_on(False)
    ax.plot(x, y, **kwargs)
    ax.set_autoscale_on(True)


def plot_stats(x):
    """Plot some statistics of x."""
    n = len(x)
    mn = np.mean(x)
    st = np.std(x)
    conf = mpa.coef_student(n, 0.95) * st / np.sqrt(n-1)

    fig, ax = plt.subplots(1, 1)

    p = ax.hist(x, rwidth=0.9)

    ax.set_xlim(0, ax.get_xlim()[1])

    ax.axvline(mn, ls='--', color='k', lw='2')
    ax.axvline(mn+st, ls='--', color='k', lw='1')
    ax.axvline(mn-st, ls='--', color='k', lw='1')

    print('Mean: {:f}'.format(mn))
    print('Std dev: {:f}'.format(st))
    print('Confidence interval (95%): [{:f}, {:f}]'.format(mn+conf, mn-conf))

    fig.tight_layout()

    return [fig, ax, p]


def get_discrete_cmap(N, base_cmap=None):
    """Return a discrete version of a cmap."""
    base = plt.cm.get_cmap(base_cmap)

    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)

    try:
        cm = base.from_list(cmap_name, color_list, N)
    except AttributeError:
        cm = plt.cm.get_cmap(base_cmap, N)

    return cm


def _get_seq_cmap(seq):
    """Return cmap from sequence (float in between colors)."""
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])

    return mc.LinearSegmentedColormap('CustomMap', cdict)


def get_lin_cmap(values, colors):
    """Return cmap that linearly evolve between a list of values."""

    # Normalized vmid
    x = [1.0*(v-values[0])/(values[-1]-values[0]) for v in values[1:-1]]

    # Adapt to seq
    c = mc.ColorConverter().to_rgb

    seq = [c(colors[0])]
    for i in range(len(x)):
        seq += [c(colors[i+1]), x[i], c(colors[i+1])]
    seq += [c(colors[-1])]

    return _get_seq_cmap(seq)


def align_twin_ticks(ax1, ax2, change2=True, N_ticks=None):
    """Align ticks of twin axes by changing view limit.

    Deactivate grid

    ax1, ax2: xaxis or yaxis
    change2: change limits of second axis provided
    N_ticks: Number of ticks wanted.
             Defaults to number of ticks on first axis
    """

    if N_ticks is None:
        N_ticks = len(ax1.get_ticklocs())

    ax1.set_major_locator(mt.MaxNLocator(N_ticks, min_n_ticks=N_ticks))
    ax2.set_major_locator(mt.MaxNLocator(N_ticks, min_n_ticks=N_ticks))

    y1_ticks = ax1.get_ticklocs()
    y2_ticks = ax2.get_ticklocs()
    ax1.set_major_locator(mt.FixedLocator(y1_ticks))
    ax2.set_major_locator(mt.FixedLocator(y2_ticks))

    # Left and right positions (1 and 2) to have aligned
    l1, l2 = y1_ticks[[1, -2]]
    r1, r2 = y2_ticks[[1, -2]]
    lmin, lmax = ax1.get_view_interval()
    rmin, rmax = ax2.get_view_interval()

    rmin = r2 - (r2-r1)/(l2-l1)*(l2-lmin)
    rmax = r2 - (r2-r1)/(l2-l1)*(l2-lmax)
    ax2.set_view_interval(rmin, rmax)

    ax2.axes.grid(False)
