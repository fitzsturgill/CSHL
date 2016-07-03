
% make a tiled lick raster figure for SessionData

function pavlovian_lick_rasters(SessionData, sessionName)

    fig = figure('Name', sessionName);
    lickRaster(SessionData, 1, 2, fig, 'toneA reward')
    lickRaster(SessionData, 1, 3, fig, 'toneA omission')
    lickRaster(SessionData, 1, 4, fig, 'toneA punish')

    lickRaster(SessionData, 2, 2, fig, 'toneB reward')
    lickRaster(SessionData, 2, 3, fig, 'toneB omission')
    lickRaster(SessionData, 2, 4, fig, 'toneB punish')

    splayAxisTile;