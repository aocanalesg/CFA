if save_v_index ==1
    pr_v_vec = reshape(pr_mat,nv*nuA*nuG*nzrs,1);
    cdf_v_vec = reshape(cdf_mat,nv*nuA*nuG*nzrs,1);
    save pr_v_vec.dat pr_v_vec -ascii -double
    save cdf_v_vec.dat cdf_v_vec -ascii -double
    save v_grid.dat v_grid -ascii -double
    vA_grid = Abar*exp(uA_grid);
    vg_grid = gbar*exp(uG_grid);
    save vA_grid.dat vA_grid -ascii -double
    save vg_grid.dat vg_grid -ascii -double
end
if save_bstar_index ==1
    save bstar_mat_N100k_mu1015_zrs40 bstar_mat
end
