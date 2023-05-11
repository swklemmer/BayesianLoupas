function u_t = create_sin_exc(img_p, exc_p)
%CREATE_SIN_EXC Create sinusoidal amplitude timelines of endogenous and
%induced displacements. The resulting displacements is expressed in
%wavelengths.

u_t = exc_p.A_n * img_p.f_c / img_p.c_c * ...
      sin(2 * pi * 2 * (0:img_p.t_max-1) / img_p.f_f) + ...
      exc_p.A_e * img_p.f_c / img_p.c_c * ...
      sin(2 * pi * exc_p.f_e * (0:img_p.t_max-1) / img_p.f_f);
end
