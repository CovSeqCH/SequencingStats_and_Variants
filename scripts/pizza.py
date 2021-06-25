# %%
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib as mpl

verts = [
    (0., 0.),  # left, bottom
    (0., 1.),  # left, top
    (1., 1.),  # right, top
    (1., 0.),  # right, bottom
    (0., 0.),  # ignored
]

codes = [
    Path.MOVETO,
    Path.LINETO,
    Path.LINETO,
    Path.LINETO,
    Path.CLOSEPOLY,
]

paths = [Path.arc(0, 40, is_wedge=True), Path.arc(70, 140, is_wedge=True)]


fig, ax = plt.subplots()
for path in paths:
    patch = patches.PathPatch(path, facecolor='orange', lw=0)
    c = ax.add_patch(patch)
    transform = mpl.transforms.Affine2D().translate(0.1, 0.1)
    c.set_transform(transform+ax.transData)
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
plt.show()
# %%
