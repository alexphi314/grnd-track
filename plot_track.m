function [] = plot_track(inputLat,inputLong)
%Given two vectors of latitudes and longitudes, plot their ground track on
%a map
%   inputLat: vector of latitudes
%   inputLong: vector of longitudes

worldmap world;
load coastlines;
plotm(coastlat,coastlon);

if length(inputLong) > 0 && length(inputLat) > 0
    plotm(inputLat,inputLong,'k');
    p1 = plotm(inputLat(1),inputLong(1),'ob','DisplayName','Starting Position');
    p2 = plotm(inputLat(end),inputLong(end),'or','DisplayName','Ending Position');
    legend([p1,p2]);
end

print('ground_track.jpg','-djpeg');

end
