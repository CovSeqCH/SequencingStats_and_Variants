# %%
from plots import load_data, restrict_dates
import numpy as np
import geopandas
from matplotlib.patches import Patch
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from surveillance_region_map import iso_canton_to_region
from variant_to_pango import variant_to_lineage
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.patches as mpatches


def draw_pie(ratios=[0.4, 0.3, 0.3], colors=["red", "blue", "green"], X=0, Y=0, size=1):
    fig = plt.gcf()
    ax = plt.gca()
    ratio_acc = 0
    for ratio, color in zip(ratios, colors):
        trans = (fig.dpi_scale_trans +
                 transforms.ScaledTranslation(X, Y, ax.transData))
        w = mpatches.Wedge((0, 0), size, 360. * ratio_acc, 360. * (ratio_acc + ratio),
                           clip_on=True, transform=trans, facecolor=color)
        ax.add_patch(w)
        ratio_acc += ratio


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
geo_ch = geopandas.read_file("swiss_map_files/gadm36_CHE_1.shp")

geo_merge = geo_ch.assign(region=lambda y: y.HASC_1.apply(
    lambda x: iso_canton_to_region[x[-2:]]))
ylgn = cm.get_cmap('YlGnBu', 512)
newcmp = ListedColormap(ylgn(np.linspace(0.25, 0.75, 256)))

ax = geo_merge.plot(column='region',
                    figsize=(10, 8), cmap=newcmp)

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
    x_val = xy_coords[f'Region {reg}']["x"]
    y_val = xy_coords[f'Region {reg}']["y"]

    draw_pie(ratios, colors, x_val, y_val, size=0.4)

variant_names = list(variant_to_lineage.keys())
variant_names.append('others')
legend_elements = [Patch(facecolor=color, label=name)
                   for color, name in zip(colors, variant_names)]

time = 'May'
plt.title(f'Variants by Region in {time} 2021')
plt.legend(handles=legend_elements, loc="upper left")
plt.axis("off")
plt.show()

# %%
