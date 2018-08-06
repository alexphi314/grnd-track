import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

def plot(inputLat,inputLong,refcoords,name,title):
    """
    Generate a plot of the ground track

    :param inputLat: array of latitudes to plot (deg)
    :param inputLong: array of longitudes to plot (deg)
    :param refcoords: lat/lon pair representing observation location (deg)
    :param name: output file name
    :param line1: first line of plot title
    :param line2: second line of plot title
    :return:
    """

    fig = plt.figure()
    m = Basemap(projection='robin', lon_0=0, resolution='i')

    m.drawcoastlines()
    m.drawparallels(np.arange(-90,90,30), labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60), labels=[0,0,1,1])

    if len(inputLat) > 0 and len(inputLong) > 0:

        split_points = [0]
        for indx, lon in enumerate(inputLong):
            if indx == 0:
                continue

            if inputLong[indx-1] < 180 and inputLong[indx] > 180:
                split_points.append(indx)

        for indx, split_point in enumerate(split_points):
            if indx == 0:
                continue

            lon = inputLong[split_points[indx-1]:split_point]
            lat = inputLat[split_points[indx-1]:split_point]

            x, y = m(lon, lat)
            m.plot(x, y, linewidth=1.5, color='b')

        lon = inputLong[split_points[-1]:]
        lat = inputLat[split_points[-1]:]
        x, y = m(lon, lat)
        m.plot(x, y, linewidth=1.5, color='b')

        m.scatter(inputLong[0], inputLat[0], latlon=True, label='Starting Position')
        m.scatter(inputLong[-1], inputLat[-1], latlon=True, c='r', label='Ending Position')

    if len(refcoords) > 1:
        m.scatter(refcoords[1], refcoords[0], latlon=True, c='r', marker='x', label='Observation Location')

    plt.legend(loc='upper right', bbox_to_anchor=(1.1,1.2))
    plt.title(title)
    fig.savefig(name, dpi=150)


