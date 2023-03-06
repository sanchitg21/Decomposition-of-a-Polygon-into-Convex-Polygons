import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LinearRing
import numpy as np

# Define the polygon coordinates in the desired order as a list of tuples

with open("Coordinates.txt", "r") as infile:
    faceCoords = []
    for line in infile:
        coords = line.strip().split(" ")
        faceCoords.append([(float(coords[i]), float(coords[i+1])) for i in range(0, len(coords), 2)])

print(faceCoords)
# faceCoords = np.array(faceCoords, dtype=float)

# Split the coordinates into separate lists of x and y coordinates
fig, ax = plt.subplots()
for l in faceCoords:
    x_coords, y_coords = zip(*l)
    x_vals = [X[0]for X in l]
    y_vals = [Y[1]for Y in l]
    polygon = Polygon(list(zip(x_vals, y_vals)))
    x,y = polygon.exterior.xy
    ax.plot(x, y, color='#000000', alpha=1,
    linewidth=2, solid_capstyle='round', zorder=2)

# Plot the polygon

plt.show()