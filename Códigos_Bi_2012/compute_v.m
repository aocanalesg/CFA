%------------------------------------------------------------
% Purpose:  kernal estimate the simulated bstar distribution
%           multi-dimension
% Last update: Aug 16, 2011
%------------------------------------------------------------



if compute_vmat_index == 1

%     % with mu=1.015 and z11 = 1-1/20--> needs theta
%     bstar_mat = bstar_mat*0.61;
%     bstar_mat = bstar_mat*0.78;
    
    % initialization
    nv = 80;                        
    bstar_low = 0.2*ybar;       % truncate the simulated data above bstar_low 
    max_v = 3*ybar; max(max(max(max(bstar_mat))));    % needs to be slightly larger than max(bstar_pos);
    min_v = max(min(min(min(min(bstar_mat)))),bstar_low);
    vstep = (max_v-min_v)/(nv-1);
    v_grid = min_v:vstep:max_v;
    
    cdf_mat = zeros(nv,nuA,nuG,nzrs);
    pr_mat = zeros(nv,nuA,nuG,nzrs);    

    for iuA = 1:nuA
        for iuG = 1:nuG
            for izrs = 1:nzrs
                bstar_vec = bstar_mat(:,iuA,iuG,izrs);

                % kernal estimation    
                bstar_pos = bstar_vec(bstar_vec>bstar_low); %     bstar_pos = bstar_vec;    
                cdf_v = ksdensity(bstar_vec,v_grid,'function','cdf');
                pdf_v = ksdensity(bstar_vec,v_grid,'function','pdf');
                cdf_v_1 = [0 cdf_v(1:end-1)];
                pr_mat(:,iuA,iuG,izrs) = cdf_v-cdf_v_1;
                cdf_mat(:,iuA,iuG,izrs) = cdf_v;

            end
        end
    end
else 
    load v_grid.dat
    nv = length(v_grid);
    load pr_v_vec.dat
    load cdf_v_vec.dat
    pr_mat = reshape(pr_v_vec,nv,nuA,nuG,nzrs);
    cdf_mat = reshape(cdf_v_vec,nv,nuA,nuG,nzrs);
end

