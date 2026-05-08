function model_stop()
% called by stopFcn callback
%{
g = evalin("base", "glogs");

%u_t = g.err.Time;
t = g.err.Time;
u   = g.u_m.Data;
yp = g.y_p.Data;
%k_t = g.kxkulog.Time;
kxku   = g.kxkulog.Data;   % 2501x4
ke = g.K_Ie.Data;
err = g.err.Data;

save('simlog2.mat','t','u','kxku','ke',"err","yp",'-v7')

%}
end

