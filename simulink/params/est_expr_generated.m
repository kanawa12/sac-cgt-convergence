function y = calc_y_generated(am,bm,kx,ku)
% Do not edit by hand.

y = -1*[-am(2,3).*kx(1,1) - am(3,3).*kx(1,2) + kx(1,3) am(1,3).*kx(1,1) am(1,3).*kx(1,2); -am(3,3).*kx(1,1) + kx(1,2) -am(3,3).*kx(1,2) + kx(1,3) am(1,3).*kx(1,1) + am(2,3).*kx(1,2); kx(1,1) kx(1,2) kx(1,3); -bm(2,1).*kx(1,1) - bm(3,1).*kx(1,2) + ku(1,1) bm(1,1).*kx(1,1) bm(1,1).*kx(1,2)]*[bm(1,1); bm(2,1); bm(3,1)] + [am(1,3); am(2,3); am(3,3); bm(1,1)];
end
