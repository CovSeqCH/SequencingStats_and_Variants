# %%
from plots import load_data, restrict_dates
import math
import numpy as np
import geopandas
from matplotlib.patches import Patch
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd
from surveillance_region_map import iso_canton_to_region
from variant_to_pango import variant_to_lineage
import matplotlib.pyplot as plt


def draw_pie(ax, ratios=[0.4, 0.3, 0.3], colors=["red", "blue", "green"], X=0, Y=0, size=1000, zorder=10):
    N = len(ratios)
    xy = []
    start = 0.
    for ratio in ratios:
        x = [0] + np.cos(np.linspace(2*math.pi*start, 2 *
                         math.pi*(start+ratio), 30)).tolist() + [0]
        y = [0] + np.sin(np.linspace(2*math.pi*start, 2 *
                         math.pi*(start+ratio), 30)).tolist() + [0]
        xy1 = list(zip(x, y))
        xy.append(xy1)
        start += ratio
    for i, xyi in enumerate(xy):
        ax.scatter([X], [Y], marker=xyi, s=size,
                   facecolor=colors[i], zorder=zorder)


# swiss_divisions = clusters["501YV1"]["Switzerland"]["observed_divisions"]
counts = {}
total = {}
other = {}

times = {
    "February": {"start": "2021-02-01", "end": "2021-02-28"},
    "March": {"start": "2021-03-01", "end": "2021-03-31"},
    "April": {"start": "2021-04-01", "end": "2021-04-30"},
    # "May": {"start": "2021-05-01", "end": "2021-05-31"}
}

xy_coords = {
    "Region 1": {"x": 6.8, "y": 46.40},
    "Region 2": {"x": 7.5, "y": 46.75},
    "Region 3": {"x": 8, "y": 47.42},
    "Region 4": {"x": 8.4, "y": 46.9},
    "Region 5": {"x": 9, "y": 47.40},
    "Region 6": {"x": 9.5, "y": 46.50}
}

text_coords = {
    "Region 1": {"x": 6.3, "y": 46.60},
    "Region 2": {"x": 7.5, "y": 46.95},
    "Region 3": {"x": 8, "y": 47.62},
    "Region 4": {"x": 8.4, "y": 47},
    "Region 5": {"x": 9, "y": 47.60},
    "Region 6": {"x": 9.5, "y": 46.70}
}
# %%
# swiss_reg = pd.read_csv("../regions.tsv", sep="\t", index_col=False)
geo_ch = geopandas.read_file("swiss_map_files/gadm36_CHE_1.shp")
# .merge(swiss_reg, on="NAME_1")
geo_merge = geo_ch.assign(region=lambda y: y.HASC_1.apply(
    lambda x: iso_canton_to_region[x[-2:]]))
# reduce the range of the colormap so not so light or so dark
# https://matplotlib.org/2.0.2/users/colormaps.html
# https://matplotlib.org/3.1.1/tutorials/colors/colormap-manipulation.html#creating-listed-colormaps
ylgn = cm.get_cmap('YlGnBu', 512)
newcmp = ListedColormap(ylgn(np.linspace(0.25, 0.75, 256)))

ax = geo_merge.plot(column='region',
                    figsize=(12, 9), cmap=newcmp)

df = load_data()
df = restrict_dates(df, '20210501', '20210601')
for reg in range(1, 7):
    data = df[df.region == reg].sum()
    ratios = []
    colors = list(map(cm.get_cmap('Accent'),
                  np.linspace(0, 1, len(variant_to_lineage)+1)))
    for variant in variant_to_lineage.keys():
        ratios.append(data[variant]/data.sequences)
    ratios.append(data['others']/data.sequences)
    print(ratios)
    x_val = xy_coords[f'Region {reg}']["x"]
    y_val = xy_coords[f'Region {reg}']["y"]
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import matplotlib as mpl

    paths = [Path.arc(0, 40, is_wedge=True), Path.arc(70, 140, is_wedge=True)]
    for path in paths:
        patch = patches.PathPatch(path, facecolor='orange', lw=0)
        c = ax.add_patch(patch)
        transform = mpl.transforms.Affine2D().translate(8, 47)
        c.set_transform(transform+ax.transData)
    # draw_pie(fig, X=x_val, Y=y_val, ratios=ratios, colors=colors, size=50000)
    # draw_pie(fig, X=x_val, Y=y_val, ratios=[
    #          0, 0.1, 0.2, 0.2, 0.5], colors=colors)

#         plt.text(x=text_coords[reg]['x'], y=text_coords[reg]['y'],
#                  s=f"{reg} ({np.sum(cnts)})", fontsize="x-large", horizontalalignment="center")


#   new_names = [clusters[clus]["display_name"] if clus in clusters else "Other" for clus in region_df.columns]
#    legend_elements = [Patch(facecolor=col, label=c) for c, col in zip(new_names, reg_clus_colors)]

#   plt.title(f'Variants by Region in {time} 2021')
#    plt.legend(handles=legend_elements, loc="upper left")
#     plt.show()
#     plt.axis("off")
#     plt.savefig(f"mapfig_{time}.pdf")

# %%

# %%
