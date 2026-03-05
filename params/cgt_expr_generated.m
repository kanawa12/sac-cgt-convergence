function y = calc_y_generated(am,bm,ap,bp)
% Do not edit by hand.

y = [am(1,3).*bp(2,1) - am(2,3).*bp(1,1) am(1,3).*bp(3,1) - am(3,3).*bp(1,1) bp(1,1) 0; am(1,3).*bp(3,1) - am(3,3).*bp(1,1) am(2,3).*bp(3,1) - am(3,3).*bp(2,1) + bp(1,1) bp(2,1) 0; bp(1,1) bp(2,1) bp(3,1) 0; bm(1,1).*bp(2,1) - bm(2,1).*bp(1,1) bm(1,1).*bp(3,1) - bm(3,1).*bp(1,1) 0 bp(1,1)] \ [am(1,3) - ap(1,3); am(2,3) - ap(2,3); am(3,3) - ap(3,3); bm(1,1)];
end
