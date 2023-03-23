# This code imports three libraries, matplotlib.pyplot, matplotlib.patches, and shapely.geometry. 
# The first two are used for creating plots and shapes, and the latter is used for creating and 
# manipulating geometric objects.
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from shapely.geometry import Polygon

# Then, the code reads in three sets of polygon coordinates from text files:
# polygonface, decomposedface, and mergedface. Each set of coordinates is read
# from a separate file and converted into a list of tuples, where each tuple
# represents a vertex of the polygon.
with open("Polygon.txt", "r") as infile:
    polygonface = []
    for line in infile:
        coords = line.strip().split(" ")
        polygonface.append([(float(coords[i]), float(coords[i+1])) for i in range(0, len(coords), 2)])

with open("Decomposed.txt", "r") as infile:
    decomposedface = []
    for line in infile:
        coords = line.strip().split(" ")
        decomposedface.append([(float(coords[i]), float(coords[i+1])) for i in range(0, len(coords), 2)])

with open("Merged.txt", "r") as infile:
    mergedface = []
    for line in infile:
        coords = line.strip().split(" ")
        mergedface.append([(float(coords[i]), float(coords[i+1])) for i in range(0, len(coords), 2)])

fig, ax = plt.subplots(1, 3, figsize=(12, 4))

# The code then creates a matplotlib figure with three subplots, each for one set of coordinates. 
# For each set of coordinates, the code creates a Polygon object using the matplotlib.patches.Polygon 
# library and sets the alpha, color, and linewidth properties. The Polygon object is then plotted on the 
# corresponding subplot using the plot function, with the x and y values obtained from the Polygon object's 
# exterior attribute. The number of faces in each set of polygons is also displayed in the subplot title.
for l in polygonface:
    x_coords, y_coords = zip(*l)
    x_vals = [X[0] for X in l]
    y_vals = [Y[1] for Y in l]
    polygon = Polygon(list(zip(x_vals, y_vals)))
    x, y = polygon.exterior.xy
    ax[0].set_xlabel("Faces:" + str(len(polygonface)))
    ax[0].plot(x, y, color='red', alpha=1, linewidth=1, solid_capstyle='round', zorder=2)
    ax[0].set_title('Input Polygon')

for l in decomposedface:
    x_coords, y_coords = zip(*l)
    x_vals = [X[0] for X in l]
    y_vals = [Y[1] for Y in l]
    polygon = Polygon(list(zip(x_vals, y_vals)))
    x, y = polygon.exterior.xy
    ax[1].set_xlabel("Faces:" + str(len(decomposedface)))
    ax[1].plot(x, y, color='green', alpha=1, linewidth=1, solid_capstyle='round', zorder=2)
    ax[1].set_title('Convex Polygons after Decomposition')

for l in mergedface:
    x_coords, y_coords = zip(*l)
    x_vals = [X[0] for X in l]
    y_vals = [Y[1] for Y in l]
    polygon = Polygon(list(zip(x_vals, y_vals)))
    x, y = polygon.exterior.xy
    ax[2].set_xlabel("Faces:" + str(len(mergedface)))
    ax[2].plot(x, y, color='blue', alpha=1, linewidth=1, solid_capstyle='round', zorder=2)
    ax[2].set_title('Convex Polygons after Merging')

# Finally, the figure is displayed using the plt.show() function.
plt.show()