"""Box certains axes in a AxesGrid."""


def get_rect(idx1, idx2, axes, hpad, vpad, offset):
    """Get xy and wh for rectangle from axes index in the grid."""
    box1, box2 = [axes[i].get_position().get_points()
                  for i in [idx1, idx2]]
    xy = box1[0]
    xy[0] -= hpad/2.3 - offset[0]
    xy[1] -= vpad/2.3 - offset[2]

    xy2 = box2[1]
    wh = [xy2[i] - xy[i] for i in range(2)]
    wh[0] += hpad/2.3 + offset[1]
    wh[1] += vpad/2.3 + offset[3]
    return xy, wh


def convert_pad(fig, hpad, vpad):
    """Convert point size in fig fraction."""
    hpad_f = hpad / fig.get_figwidth()
    vpad_f = vpad / fig.get_figheight()
    return hpad_f, vpad_f
