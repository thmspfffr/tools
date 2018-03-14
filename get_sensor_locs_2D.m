function sa = get_sensor_locs_2D(fn)
%fn: path to MEG data

sens = ft_read_sens(fn);

inds = [];
n=length(sens.label);
for i=1:n;
    x=sens.label{i};
    if x(1)=='M';
        inds=[inds,i];
    end
end

sa.locs_3D_indi=sens.coilpos(inds,:);
para1.rot=90;
locs_2D=mk_sensors_plane(sa.locs_3D_indi,para1);
para2.rot=90;
para2.nin=50;
locs_2D_sparse=mk_sensors_plane(sa.locs_3D_indi,para2);
sa.locs_2D=locs_2D;
sa.locs_2D_sparse=locs_2D_sparse;