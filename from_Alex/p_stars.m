function stars = p_stars(p)
% Returns stars matching sum( p < [ 0.05, 0.1, 0.001 ])
stars = repmat('*',1,sum(p< [0.05 0.01 0.001]));

end