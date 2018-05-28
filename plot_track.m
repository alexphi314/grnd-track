function [] = plot_track(inputLat,inputLong)
%Given two vectors of latitudes and longitudes, plot their ground track on
%a map
%   inputLat: vector of latitudes
%   inputLong: vector of longitudes

worldmap world;
load coastlines;
plotm(coastlat,coastlon);

if length(inputLong) > 0 && length(inputLat) > 0
    plotm(inputLong,inputLat);
end

print('ground_track.jpg','-djpeg');

end

