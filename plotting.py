# Load modules
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath

# Matplotlib colors
npcols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# Define map projections
data_crs = ccrs.PlateCarree()
proj_Global = ccrs.EqualEarth(central_longitude=200)
proj_North = ccrs.NorthPolarStereo()
proj_South = ccrs.NorthPolarStereo()
proj_Global_PC = ccrs.PlateCarree(central_longitude=210)

def setup_figure(projection,nCols=1,nRows=1,size=(15,15),mask=True):

    circ = False
    if projection.lower() == 'north':
        map_proj = proj_North
        circ = True
    elif projection.lower() == 'south':
        map_proj = proj_South
        circ = True
    elif projection.lower() == 'global':
        map_proj = proj_Global
    elif projection.lower() == 'global_pc':
        map_proj = proj_Global_PC
    else:
        raise NotImplementedError("Projection not implemented.")

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)


    land_feature = cfeature.NaturalEarthFeature(
        category='physical',
        name='land',
        scale='50m',
        facecolor=(0.8,0.8,0.8),
        edgecolor='k',
        linewidth=.25,
    )
    

    fig = plt.figure(figsize=size)
    if nCols==1 and nRows==1:
        ax = fig.add_subplot(1, 1, 1, projection=map_proj)

        ax.coastlines('50m', linewidth=0.8)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.gridlines(
            draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
        )

        if circ:
            ax.set_boundary(circle, transform=ax.transAxes)
        if mask:
            ax.add_feature(land_feature)
    elif nCols > 1 and nRows > 1:
        ax = np.empty((nCols,nRows),dtype=object)
        iF = 1
        for jj in range(nRows):
            for ii in range(nCols):
                ax[ii,jj] = fig.add_subplot(nRows, nCols, iF, projection=map_proj)

                ax[ii,jj].coastlines('50m', linewidth=0.8)
                ax[ii,jj].tick_params(axis='both', which='major', labelsize=10)
                ax[ii,jj].gridlines(
                    draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
                )
                if circ:
                    ax[ii,jj].set_boundary(circle, transform=ax.transAxes)
                if mask:
                    ax[ii,jj].add_feature(land_feature)
                iF += 1
    elif nCols > 1:
        ax = np.empty((nCols),dtype=object)
        for ii in range(nCols):
            ax[ii] = fig.add_subplot(1, nCols, ii+1, projection=map_proj)

            ax[ii].coastlines('50m', linewidth=0.8)
            ax[ii].tick_params(axis='both', which='major', labelsize=10)
            ax[ii].gridlines(
                draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
            )
            if circ:
                ax[ii].set_boundary(circle, transform=ax.transAxes)
            if mask:
                ax[ii].add_feature(land_feature)
    elif nRows > 1:
        ax = np.empty((nRows),dtype=object)
        for ii in range(nRows):
            ax[ii] = fig.add_subplot(nRows, 1, ii+1, projection=map_proj)

            ax[ii].coastlines('50m', linewidth=0.8)
            ax[ii].tick_params(axis='both', which='major', labelsize=10)
            ax[ii].gridlines(
                draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--'
            )
            if circ:
                ax[ii].set_boundary(circle, transform=ax.transAxes)
            if mask:
                ax[ii].add_feature(land_feature)

            
    return fig, ax