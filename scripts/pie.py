# %%
import numpy as np
import geopandas
from matplotlib.patches import Patch
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, to_rgb
from .surveillance_region_map import iso_canton_to_region
from .variant_to_pango import variant_to_lineage
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.patches as mpatches
from .plots import save_fig
from .colors import variant_to_color


def draw_pie(ratios, colors, x=0, y=0, size=1):
    fig = plt.gcf()
    ax = plt.gca()
    ratio_acc = 0
    for ratio, color in zip(ratios, colors):
        trans = (fig.dpi_scale_trans +
                 transforms.ScaledTranslation(x, y, ax.transData))
        w = mpatches.Wedge((0, 0), size, 360. * ratio_acc, 360. * (ratio_acc + ratio),
                           clip_on=True, transform=trans, facecolor=color)
        ax.add_patch(w)
        ratio_acc += ratio


def generate_variant_map(df, output_dir):

    xy_coords = {
        "Region 1": {"x": 6.8, "y": 46.40},
        "Region 2": {"x": 7.5, "y": 46.75},
        "Region 3": {"x": 8, "y": 47.42},
        "Region 4": {"x": 8.4, "y": 46.9},
        "Region 5": {"x": 9, "y": 47.40},
        "Region 6": {"x": 9.5, "y": 46.50}
    }

    geo_ch = geopandas.read_file("scripts/swiss_map_files/gadm36_CHE_1.shp")

    geo_merge = geo_ch.assign(region=lambda y: y.HASC_1.apply(
        lambda x: iso_canton_to_region[x[-2:]]))
    ylgn = cm.get_cmap('YlGnBu', 512)
    newcmp = ListedColormap(ylgn(np.linspace(0.25, 0.75, 256)))

    geo_merge.plot(column='region',
                   figsize=(10, 8), cmap=newcmp)

    variant_names = list(variant_to_lineage.keys())
    variant_names.append('others')

    total_cases = df[df.region == 0].cases.sum()
    for reg in range(1, 7):
        data = df[df.region == reg].sum()
        ratios = []
        colors = list(map(cm.get_cmap('Set1'),
                          np.linspace(0, 0.775, len(variant_to_lineage)+1)))
        if len(variant_to_lineage) == 5:
            colors = [(0.819607843137255, 0.4, 0.4), (1, 0.4, 0.4), (1, 0.701960784313725, 0.701960784313725), (0.4, 0.63921568627451, 0.4), (0.956862745098039, 0.862745098039216, 0.505882352941176),(0,0,0)]
        
        for i, variant in enumerate(variant_names):
            if variant in variant_to_color.keys():
                colors[i] = to_rgb(variant_to_color[variant])
        
        for variant in variant_to_lineage.keys():
            ratios.append(data[variant]/data.sequences)
        ratios.append(data['others']/data.sequences)
        x_val = xy_coords[f'Region {reg}']["x"]
        y_val = xy_coords[f'Region {reg}']["y"]

        draw_pie(ratios, colors, x_val, y_val,
                 size=np.sqrt(data.cases/total_cases))


    legend_elements = [Patch(facecolor=color, label=name)
                       for color, name in zip(colors, variant_names)]

    plt.title(
        f'Variants by Region in {df.index[10].strftime("%B %Y")}', size=15)
    plt.figtext(0.6, 0.2, 'Pie area represents number of cases')
    plt.legend(handles=legend_elements, loc="upper left")
    plt.axis("off")
    save_fig(plt.gcf(), output_dir, 'variant_map')
