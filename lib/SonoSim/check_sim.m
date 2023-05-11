function check_sim(img_p, exc_p, est_p, u_true)
%CHECK_SIM Raises warnings if simulation parameters may lead to errors.

% Check axial kernel length depending on coarse estimation's reference
if est_p.crs_mref

    % Calculate minimal kernel length based on differential displacement
    u_diff = diff(u_true, 1, 3);
    z_min = ceil(max(abs(u_diff), [], 'all'));

    if est_p.z_len < z_min % [wvls]
        warning("Axial kernel size is smaller than max. diff. displacement. " + ...
            "Has to be at least %d wvls.", z_min)
    end

else
    % Calculate minimal kernel length
    z_min = ceil(exc_p.A_n / img_p.c_c * img_p.f_c); % [wvls]

    if est_p.z_len < z_min % [wvls]
        warning("Axial kernel size is smaller than max. displacement. " + ...
            "Has to be at least %d wvls.", z_min)
    end
end

end
