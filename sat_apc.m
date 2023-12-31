function [apc_sat] = sat_apc(rsat,rrcv,rsun,s_apc,satnum,opt)

% opt=1 for L1    &    opt=2 for L2

k=(-1).*(rsat./(norm(rsat)));
rs  = rsun - rsat;
e   = rs./(norm(rs));
j   = cross(k,e);
i   = cross(j,k);

R=[i ; j ; k];

apc_sat=R'*s_apc(satnum, : , opt)';

end