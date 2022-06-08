
def setnan(var):
    """
    Set the NaN values to the masked array variables of the DPR 2A products 
    according to the file specifications.

    Parameters
    ----------
    var : masked array
        The 2A DPR variables inside the HDF file. Can be from the instruments 
        DPR, Ka or Ku.

    Returns
    -------
    var : masked array
        The same variable with the NaNs set as masked values.
    """
    import numpy as np
    
    if (var.dtype == 'float32') | (var.dtype == 'float64'):
        var[var == -1111.1] = np.nan
        var[var == -9999.9] = np.nan
    
    if (var.dtype == 'int16') | (var.dtype == 'int32') | (var.dtype == 'uint8'):
        var = var.astype('float')
        var[var == -1111.] = np.nan
        var[var == -9999.] = np.nan
        #var[var == 0] = np.nan
        var[var == 255] = np.nan
    var = np.ma.masked_where(np.isnan(var), var)
    return var


def bin_to_height(binRangeNo, localZenithAngle, ellipsoidBinOffset, geoid_Offset,
                  binEllipsoid=176, rangeBinSize=125, band='DPR'):
    """
    Convert the bin range distance to elevation in m.a.s.l.
    It works for Ku, Ka, and DPR level 2A products. The default parameters are set for DPR.

    Parameters
    ----------
    binRangeNo : int
        Bin range number to convert. 0 is close to the ground and 176 is about 20km AGL for DPR.
    localZenithAngle : masked array
        Local zenital angle for each pixel in degree. Provided by the DPR variable with the same name.
    ellipsoidBinOffset : masked array
        The offset of the first bin due to its inclination. Provided by the DPR variable with the same name.
    geoid_Offset : float
        Offset in meters between the geoid and the ellipsoid (supposed WGS84) in the region study. 
        Considered regionally constant here.
    binEllipsoid : int, optional
        The number of bin gates along the radar beam. The default is 176 for 2ADPR.
    rangeBinSize : float, optional
        The size of a bin in meter. The default is 125.
    band : string, optional
        The name of the band. Can be used for 'DPR, 'Ku', or 'Ka'.
        The default is 'DPR'.

    Returns
    -------
    masked array
        A georefferenced vector of values in m.a.s.l of the 2ADPR variable.

    """
    import numpy as np
    
    if band=='DPR' :
        return ((binEllipsoid - binRangeNo)*rangeBinSize + ellipsoidBinOffset)* np.cos(np.deg2rad(localZenithAngle[:,:,0])) + geoid_Offset
    elif band=='Ku' :
        return ((binEllipsoid - binRangeNo)*rangeBinSize + ellipsoidBinOffset)* np.cos(np.deg2rad(localZenithAngle[:,:])) + geoid_Offset
    elif band=='Ka' :
        return ((binEllipsoid - binRangeNo)*rangeBinSize + ellipsoidBinOffset)* np.cos(np.deg2rad(localZenithAngle[:,:])) + geoid_Offset
    else :
        print('Error: band issue')
        

def data_circle(data, latitude, longitude, lon_loc, lat_loc, radius):
    """
    Masked the data of an array outside the circle centered in (lon_loc,lat_loc).
    
    Parameters
    ----------
    data : Array
        The array of variable.
    latitude : Array
        latitude georeferenced.
    longitude : Array
        longitude georeferenced.
    lon_loc : float
        longitude of the circle center in degree.
    lat_loc : float
        Latitude of the circle center in degree.
    radius : float
        Value of the radius in degree.

    Returns
    -------
    masked array
        Masked array of the data where data outside the circle are masked.

    """
    import numpy as np
    import numpy.ma as ma
    
    #mask the values outsides the radius
    return ma.masked_where(np.sqrt((longitude-lon_loc)**2+(latitude-lat_loc)**2)>radius, data)


def single_pt_haversine(lat, lng, degrees=True):
    """
    'Single-point' Haversine: Calculates the great circle distance
    between a point on Earth and the (0, 0) lat-long coordinate
    """
    
    from math import radians, cos, sin, asin, sqrt
    
    r = 6371 # Earth's radius (km). Have r = 3956 if you want miles

    # Convert decimal degrees to radians
    if degrees:
        lat, lng = map(radians, [lat, lng])

    # 'Single-point' Haversine formula
    a = sin(lat/2)**2 + cos(lat) * sin(lng/2)**2
    d = 2 * r * asin(sqrt(a)) 

    return d


# function for plotting variables.
# Must be adapted according to the variables and the region

def proplot_map_var(data1, data2, data_elev, longitude, latitude, Stations,
                    path_dir_out, date, hour, 
                    lon_loc, lat_loc, radius,
                    band='DPR', cat='FS', name_station='Alps'):
    """
    Plot the variables goerefrenced of the 2ADPR products.
    This is coded for the height of the bright band (heightBB)
    and the phase near the surface (phaseNS) above French Alps 
    but can be easily adapted to other variables in other regions.

    Parameters
    ----------
    data1 : masked array
        2-dim 2ADPR variable to show. Initially coded here for heightBB.
    data2 : masked array
        2-dim 2ADPR variable to show. Initially coded here for phaseNS.
    data_elev : masked array
        Elevation of the ground (DEM) given by the HDF File.
    longitude : masked array
        Longitude array of the HDF File.
    latitude : masked array
        Latitude array of the HDF File.
    Stations : DataFrame
        Dataframe that contains the coordinates of the 
        weather stations around Grenoble.
    path_dir_out : string
        The path where to save the figures.
    date : string
        Date of the satellite overpass.
    hour : string
        Time of the satellite overpass.
    lon_loc : float
        Longitude of the circle center.
    lat_loc : float
        Latitude of the circle center.
    radius : float
        Radius in degree of the circle centered in (lat_loc,lon_loc).
    band : string, optional
        3 instruments are accepted: 'DPR', 'Ka', or 'Ku'. 
        The default is 'DPR'.
    cat : string, optional
        The subgroup of the HDF File. It was initially coded for 'FS'
        but can be adapted for 'HS' or 'MS'. The default is 'FS'.
    name_station : string, optional
        Name of the region you are looking at. The default is 'Alps'.

    Returns
    -------
    None.

    """
    import proplot
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    import matplotlib as mpl
    import matplotlib.patches as mpatches
    import cartopy.crs as ccrs
    
    # Plate Carrée map projection
    proplot.rc.reso = 'hi' 

    cmap_phase = proplot.Colormap('blue','green','red', discrete=True)
    cmap = ListedColormap(["light blue", 'light green', "light red", 'blue'])
    cmap_heightBB = proplot.Colormap('viridis')
    cmap_heightBB.set_under(color='light grey', alpha=0.8)
    bounds = [0., 100., 200., 300.]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='max')
    cmap_elev = proplot.Colormap('thermal')
    levels = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
    norm_elev = mpl.colors.BoundaryNorm(levels, cmap_elev.N)
    
    #fig, axs = proplot.subplots(nrows=2, axwidth=5, proj=('cyl', proj))
    fig, axs = proplot.subplots(nrows=1, ncols=2, axwidth=5, proj=('cyl'))
    axs.format(
        land=False, labels=True, borders=True,
        ocean=True, oceancolor='cyan', oceanalpha=0.2, oceanzorder=0.5,
        coast=True, #innerborders=True
        lonlines=1, latlines=1,
        gridminor=True, 
        suptitle='GPM '+name_station+' '+band+'_'+cat+' '+'rad'+str(radius)+' '+
        date+' '+hour+' UTC'
    )
    
    axs[0].add_patch(mpatches.Circle(xy=[lon_loc, lat_loc], 
                                  radius=radius, 
                                  facecolor='none', 
                                  edgecolor = 'k',
                                  # alpha=0.3, 
                                  transform=ccrs.PlateCarree(), 
                                  zorder=30))
    axs[0].plot(lon_loc, lat_loc, 'k+')
    #axs[0].plot(longitude.ravel()[2419], latitude.ravel()[2419], 'r+')
    axs[0].format(
        #lonlim=(67, 71), latlim=(-48, -51), #Kerguelen
        lonlim=(4, 7.5), latlim=(43, 46),  #Alps french
        #lonlim=(-79, -76), latlim=(-7.5, -10), #Huaraz
        #lonlim=(7, 10.5), latlim=(47, 50), #Stuttgart
        labels=True, title= 'BB Bottom'
    )
    axs[1].format(
        lonlim=(4, 7.5), latlim=(43, 46),
        labels=True, title='PhaseNearSurface'
    )

    axs[0].pcolormesh(longitude, latitude, data1, alpha=0.8, vmin=500,
                      cmap=cmap_heightBB,colorbar='b', shading='auto', levels=levels,
                      colorbar_kw={'title': 'Bright Band elevation [m.a.s.l]'})

    phaseplot = axs[1].pcolormesh(longitude, latitude, data2, alpha=0.8,
                      cmap=cmap, #levels=[0,100,200,300]#, norm=norm
                      norm=norm,colorbar='b',
                      colorbar_kw={'title': 'Phase Type',
                                   'ticklabels':['Solid','Mix','Liquid'],
                                   'ticks': [50,150,250]}
                      )
    # fig.colorbar(phaseplot, #ticks=50, 
    #              loc='b', values=[50,150,250],
    #                         #orientation='horizontal', 
    #                         label='phase', col=2, 
    #                         tickloc='bottom')
    #cbar1.axs[1].set_xticklabels(['< -1', '0', '> 1'])
    try:
        axs.contour(longitude, latitude, data_elev, alpha=0.8,
            cmap= cmap_elev, norm=norm_elev,
            labels=True, lw=0.8, ls='--', labels_kw={'size': 7}
        )
        fig.colorbar(cmap_elev, norm=norm_elev, loc='b', 
                     label='Ground elevation [m.a.s.l]')
    except : 
        print('Trouble with the elevation data')
    axs.scatter(Stations.Long, Stations.Lat, c='green5', marker='^',s=100, mec='k')
    fig.show()

    fig.savefig(filename= path_dir_out +date+'_'+hour+'_' + 
                band+'_'+cat+'_rad'+str(radius)+'.png',dpi=300, format= 'png')
    #plt.close(fig)

def proplot_map_BB(data1, data2, data3, data_elev, longitude, latitude, 
                   Stations, path_dir_out, date, hour,
                   lon_loc, lat_loc, radius,
                   band='DPR', cat='FS', name_station='Alps'):
    """
    Plot the bright band elevations (height, top and bottom) of the HDF File
    above the region of Grenoble (location can be modified). Also plot the 
    countour elevation of the ground in the same region.

    Parameters
    ----------
    data1 : masked array
        HeightBB variable, the 2-dim array of height of the bright band 
        in m.a.s.l.
    data2 : masked array
        TopBB. The bright band top. From the bin_BBTop 2ADPR variables
        but previously rescaled to have m.a.s.l values.
    data3 : masked array
        BotBB. The bright band bottom. From the bin_BBBottom 2ADPR variables
        but previously rescaled to have m.a.s.l values.
    data_elev : masked array
        Elevation of the ground (DEM) given by the HDF File.
    longitude : masked array
        Longitude array of the HDF File.
    latitude : masked array
        Latitude array of the HDF File.
    Stations : DataFrame
        Dataframe that contains the coordinates of the 
        weather stations around Grenoble.
    path_dir_out : string
        The path where to save the figures.
    date : string
        Date of the satellite overpass.
    hour : string
        Time of the satellite overpass.
    lon_loc : float
        Longitude of the circle center.
    lat_loc : float
        Latitude of the circle center.
    radius : float
        Radius in degree of the circle centered in (lat_loc,lon_loc).
    band : string, optional
        3 instruments are accepted: 'DPR', 'Ka', or 'Ku'. 
        The default is 'DPR'.
    cat : string, optional
        The subgroup of the HDF File. It was initially coded for 'FS'
        but can be adapted for 'HS' or 'MS'. The default is 'FS'.
    name_station : string, optional
        Name of the region you are looking at. The default is 'Alps'.
        
    Returns
    -------
    None.

    """
    import proplot
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    import matplotlib as mpl
    import matplotlib.patches as mpatches
    import cartopy.crs as ccrs
    import numpy as np
    
    # Plate Carrée map projection
    proplot.rc.reso = 'hi' 

    cmap_phase = proplot.Colormap('blue','green','red', discrete=True)
    cmap = ListedColormap(["light blue", 'light green', "light red", 'blue'])
    cmap_heightBB = proplot.Colormap('viridis')
    cmap_heightBB.set_under(color='light grey', alpha=0.8)
    bounds = [0., 100., 200., 300.]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='max')
    cmap_elev = proplot.Colormap('thermal')
    #levels = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
    levels = np.arange(0,4250,250)
    levels_BB = np.arange(1000,3200,200)
    norm_elev = mpl.colors.BoundaryNorm(levels, cmap_elev.N)
    
    #fig, axs = proplot.subplots(nrows=2, axwidth=5, proj=('cyl', proj))
    fig, axs = proplot.subplots(nrows=1, ncols=3, axwidth=4, proj=('cyl'))
    axs.format(
        land=False, labels=True, borders=True,
        ocean=True, oceancolor='cyan', oceanalpha=0.2, oceanzorder=0.5,
        coast=True, #innerborders=True
        lonlines=1, latlines=1,
        gridminor=True, 
        #suptitle='GPM '+name_station+' '+band+'_'+cat+' '+'rad'+str(thr_dist)+' '+date+' '+hour+' UTC'
    )
    
    for i in range(0,len(axs)):
        axs[i].add_patch(mpatches.Circle(xy=[lon_loc, lat_loc], 
                                      radius=radius, 
                                      facecolor='none', 
                                      edgecolor = 'k',
                                      # alpha=0.3, 
                                      transform=ccrs.PlateCarree(), 
                                      zorder=30))
        axs[i].plot(lon_loc, lat_loc, 'k+')
    #axs[0].plot(longitude.ravel()[2419], latitude.ravel()[2419], 'r+')
    axs[0].format(
        #lonlim=(67, 71), latlim=(-48, -51), #Kerguelen
        lonlim=(4, 7.5), latlim=(43, 46),  #Alps french
        #lonlim=(-79, -76), latlim=(-7.5, -10), #Huaraz
        #lonlim=(7, 10.5), latlim=(47, 50), #Stuttgart
        labels=True, title= 'BB Height'
    )
    axs[1].format(
        lonlim=(4, 7.5), latlim=(43, 46),
        labels=False, title='BB Top'
    )
    axs[2].format(
        lonlim=(4, 7.5), latlim=(43, 46),
        labels=False, title='BB Bottom'
    )
    axs[0].pcolormesh(longitude, latitude, data1, alpha=0.8, vmin=500,
                      cmap=cmap_heightBB, shading='auto', levels=levels_BB)
    axs[1].pcolormesh(longitude, latitude, data2, alpha=0.8, vmin=500,
                      cmap=cmap_heightBB, shading='auto', levels=levels_BB)
    axs[2].pcolormesh(longitude, latitude, data3, alpha=0.8, vmin=500,
                      cmap=cmap_heightBB, shading='auto', levels=levels_BB,
                      colorbar='r',colorbar_kw={'title': 'Bright Band elevation [m.a.s.l]'})

    # phaseplot = axs[1].pcolormesh(longitude, latitude, data2, alpha=0.8,
    #                   cmap=cmap, #levels=[0,100,200,300]#, norm=norm
    #                   norm=norm,colorbar='b',
    #                   colorbar_kw={'title': 'Phase Type',
    #                                'ticklabels':['Solid','Mix','Liquid'],
    #                                'ticks': [50,150,250]}
    #                   )
    # fig.colorbar(phaseplot, #ticks=50, 
    #              loc='b', values=[50,150,250],
    #                         #orientation='horizontal', 
    #                         label='phase', col=2, 
    #                         tickloc='bottom')
    #cbar1.axs[1].set_xticklabels(['< -1', '0', '> 1'])

    try:
        axs.contour(longitude, latitude, data_elev, alpha=0.8,
            cmap= cmap_elev, norm=norm_elev,
            labels=True, lw=0.8, ls='--', labels_kw={'size': 7}
            # labels_kw={'weight': 'bold'}
            #levels = [1000, 2000, 3000]
        )
        fig.colorbar(cmap_elev, norm=norm_elev, loc='b', 
                     label='Ground elevation [m.a.s.l]')
    except : 
        print('Trouble with the elevation data')
    axs.scatter(Stations.Long, Stations.Lat, c='green5', marker='^',s=100, mec='k')
    fig.show()
    fig.savefig(filename= path_dir_out +date+'_'+hour+'_' + 
                band+'_'+cat + 'BB.png',dpi=300, format= 'png')
    # fig.savefig(filename= path_dir_out +date+'_'+hour+'_' + 
    #             band+'_'+cat+'_rad'+str(thr_dist)+'.png',dpi=100, format= 'png')
    # plt.close(fig)
    
    
def proplot_map_snow(data, data_elev, longitude, latitude, Stations,
                    path_dir_out, date, hour,
                    band='DPR', cat='FS', name_station='Alps'):
    """
    Plot the presence of snowfall detected by 2ADPR algorithm during 
    a precipitation event above the French Alps. It can be easily adapted to 
    represent any binary georeferenced variables.

    Parameters
    ----------
    data : masked array
        2-dim 2ADPR variable to show. Initially coded to represent flagSnow,
        indicating the presence of snowfall near the ground (0: no Snow, 1: Snow).
    data_elev : masked array
        Elevation of the ground (DEM) given by the HDF File.
    longitude : masked array
        Longitude array of the HDF File.
    latitude : masked array
        Latitude array of the HDF File.
    Stations : DataFrame
        Dataframe that contains the coordinates of the 
        weather stations around Grenoble.
    path_dir_out : string
        The path where to save the figures.
    date : string
        Date of the satellite overpass.
    hour : string
        Time of the satellite overpass.
    band : string, optional
        3 instruments are accepted: 'DPR', 'Ka', or 'Ku'. 
        The default is 'DPR'.
    cat : string, optional
        The subgroup of the HDF File. It was initially coded for 'FS'
        but can be adapted for 'HS' or 'MS'. The default is 'FS'.
    name_station : string, optional
        Name of the region you are looking at. The default is 'Alps'.

    Returns
    -------
    None.

    """
    import proplot
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    import matplotlib as mpl
    
    # Plate Carrée map projection
    proplot.rc.reso = 'hi' 

    cmap_snow = ListedColormap(["light grey", 'cyan'])
    bounds = [0., 0.5, 1.]
    norm = mpl.colors.BoundaryNorm(bounds, cmap_snow.N, extend='neither')
    cmap_elev = proplot.Colormap('thermal')
    levels = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
    norm_elev = mpl.colors.BoundaryNorm(levels, cmap_elev.N)
    
    #fig, axs = proplot.subplots(nrows=2, axwidth=5, proj=('cyl', proj))
    fig, axs = proplot.subplots(nrows=1, ncols=1, axwidth=5, proj=('cyl'))
    axs.format(
        land=False, labels=True, borders=True,
        ocean=True, oceancolor='cyan', oceanalpha=0.2, oceanzorder=0.5,
        coast=True, #innerborders=True
        lonlines=1, latlines=1,
        gridminor=True, suptitle='GPM '+name_station+' '+band+' '+cat+' '+
        date+' '+hour+' UTC'
    )
    
    axs[0].format(
        lonlim=(4, 7.5), latlim=(43, 46),
        labels=True, title='Snowfall on the ground')

    axs[0].pcolormesh(longitude, latitude, data, alpha=0.8,
                      cmap=cmap_snow, colorbar='b', shading='auto',norm=norm,
                       colorbar_kw={'title': 'Snowfall on the ground',
                                    'ticklabels':['No','Yes'],
                                    'ticks': [0.25,0.75]}
                      )

    try:
        axs.contour(longitude, latitude, data_elev, alpha=0.8,
            cmap= cmap_elev, norm=norm_elev,
            labels=True, lw=0.8, ls='--', labels_kw={'size': 7}
        )
        fig.colorbar(cmap_elev, norm=norm_elev, loc='b', 
                     label='Ground elevation [m.a.s.l]')
        print('contourplot')
    except : 
        print('Trouble with the elevation data')
    axs.scatter(Stations.Long, Stations.Lat, c='green5', marker='^',s=60, mec='k')
    fig.show()
    fig.savefig(filename= path_dir_out +date+'_'+hour+'_' + 
                band+'_'+cat + '_Snow.png', format= 'png')
    # plt.close(fig)