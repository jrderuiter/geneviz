from collections import OrderedDict

import seaborn as sns


def build_colormap(data, hue, palette=None, order=None):
    """Builds a colormap, mapping data hue values to colors."""

    if hue is None:
        color_map = None
    else:
        if palette is None:
            palette = sns.color_palette()

        if order is None:
            order = data[hue].unique()

        color_map = OrderedDict(zip(order, palette))

    return color_map
