import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import math as m
from extern_data import data

plt.style.use('dark_background')

datafile = '/Users/Desktop/asc/TARGET__00285.fit'
imdata = fits.open(datafile)
imgdata = np.fliplr(imdata[0].data)
imgdata = imgdata[:1005, 165:1205]
lx, ly = imgdata.shape
print('x: ' + str(ly) + ' y: ' + str(lx))

# Create a circular mask
X, Y = np.ogrid[0:lx, 0:ly]

A = 16  # Y, greater value moves down
B = 0  # X
C = 4.213  # radius  greater value decreases radius
circ_mask = (A + (X - lx / 2)) ** 2 + (B + (Y - ly / 2)) ** 2 > ((lx * ly) / C)
imgdata[circ_mask] = np.min(imgdata)

imgmin = np.min(imgdata)
imgmax = np.max(imgdata)
imgmean = np.mean(imgdata)
imgmed = np.median(imgdata)

print('min: ' + str(imgmin) + ' max: ' + str(imgmax) + ' mean: ' + '%06.02f' % imgmean + ' median: ' + str(imgmed))

nummax = imgdata[np.where(imgdata > (imgmed + 1500.0))]
print('number of nummax: ' + str(len(nummax)))

r1 = 510.0
r2 = 500.0
xc = 520.0
yc = 486.0

theta = np.linspace(0, 2 * np.pi, 100)
cir_x1 = xc + (r1 + 2) * np.cos(theta)
cir_y1 = yc + (r1 + 2) * np.sin(theta)

cir_x2 = xc + r2 * np.cos(theta)
cir_y2 = yc + r2 * np.sin(theta)

x1l, x2l, x3l, x4l, y1l, y2l, y3l, y4l = [], [], [], [], [], [], [], []

offset = 0.5  # Rotation offset applies to every line in each grid. Smaller value moves clockwise.
for angle in range(0, 360):
    x1 = xc + (r1 * m.sin(m.radians(angle - offset)))
    y1 = yc + (r1 * m.cos(m.radians(angle - offset)))

    x2 = xc + (r2 * m.sin(m.radians(angle - offset)))
    y2 = yc + (r2 * m.cos(m.radians(angle - offset)))

    x3 = xc + (0 * m.sin(m.radians(angle - offset)))
    y3 = yc + (0 * m.cos(m.radians(angle - offset)))
    x4 = xc + (r1 * m.sin(m.radians(angle - offset)))
    y4 = yc + (r1 * m.cos(m.radians(angle - offset)))

    x1l.append(x1)
    y1l.append(y1)
    x2l.append(x2)
    y2l.append(y2)
    x3l.append(x3)
    x4l.append(x4)
    y3l.append(y3)
    y4l.append(y4)

x1l = np.asarray(x1l)
x2l = np.asarray(x2l)
y1l = np.asarray(y1l)
y2l = np.asarray(y2l)
x3l = np.asarray(x3l)
x4l = np.asarray(x4l)
y3l = np.asarray(y3l)
y4l = np.asarray(y4l)

# Set up the plot
figure = plt.figure(figsize=(11.04, 11.04), dpi=100)
plt.gca().set_aspect('equal')
ax = plt.gca()
ax.axis('off')
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_aspect(1)
plt.gca().set_aspect('equal')

ax.plot((x3l[::90], x4l[::90]), (y3l[::90], y4l[::90]), linestyle='-', color='w', linewidth=0.9, alpha=0.5)

# Plot compass "spines"
ax.plot((x1l, x2l), (y1l, y2l), linestyle='-', color='orange', linewidth=0.8, alpha=0.9, zorder=1)
ax.plot((x1l[::5], x2l[::5]), (y1l[::5], y2l[::5]), linestyle='-', color='r', linewidth=1, alpha=1, zorder=2)

# Compass circles
ax.plot(cir_x1, cir_y1, linestyle='-', linewidth=0.7, color='orange')
ax.plot(cir_x2, cir_y2, linestyle='-', linewidth=0.7, color='orange')

ax.set_xlim([0, 1040])
ax.set_ylim([-63, 1040])
plt.imshow(imgdata[:, :], cmap='gray', origin='lower', vmin=imgmin + 133, vmax=np.median(nummax))

obsname = data.external_data()[0]
location = data.external_data()[1]
temperature = data.external_data()[2]
rel_humidity = data.external_data()[3]
exposure = data.external_data()[4]
date = data.external_data()[5]
time = data.external_data()[6]
lst = data.external_data()[7]
moon_illum = data.external_data()[8]
moon_elev = data.external_data()[9]
sun_elev = data.external_data()[10]
seqnum = data.external_data()[11]
heatstat = data.external_data()[12]

deg = u"\u00b0"
font = {'weight': 'heavy'}
plt.rc('font', **font)

plt.text(-30, 1030, obsname, color='purple', fontsize=12, alpha=1)
plt.text(-30, 1010, location, color='#696ca5', fontsize=12, alpha=1)
#
plt.text(920, 1030, 'Temp = ' + str(temperature) + ' ' + deg + 'C', color='#696ca5', fontsize=12, alpha=1)
plt.text(971, 1010, 'RH = ' + str(rel_humidity) + '%', color='yellow', fontsize=12, alpha=1)
plt.text(870, 990, 'Exposure = ' + str(exposure) + ' secs', color='#696ca5', fontsize=12, alpha=1)
plt.text(850, 970, 'Date = ' + str(date) + ' UTC', color='#696ca5', fontsize=12, alpha=1)
plt.text(870, 950, 'Time = ' + str(time) + ' UTC', color='#696ca5', fontsize=12, alpha=1)
plt.text(930, 930, 'LST = ' + str(lst), color='#696caf', fontsize=12, alpha=1)
#
plt.text(885, 20, 'Moon illum = ' + str(moon_illum) + '%', color='yellow', fontsize=12, alpha=0.6)
plt.text(900, 0, 'Moon elev = ' + str(moon_elev) + deg, color='yellow', fontsize=12, alpha=0.6)
plt.text(910, -20, 'Sun elev = ' + str(sun_elev) + deg, color='yellow', fontsize=12, alpha=0.6)
plt.text(920, -40, 'Seq num = ' + str(seqnum), color='#696ca5', fontsize=12, alpha=1)
plt.text(900, -60, 'Dew Heater is ' + str(heatstat), color='r', fontsize=12, alpha=1)

plt.tight_layout()
plt.show()
