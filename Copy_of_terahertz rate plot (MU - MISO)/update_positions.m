function para = update_positions(para)
    n = 0:para.M_t-1;
    para.y_bs_t = (n - (para.M_t-1)/2) * para.d;

    m = 0:para.M_r-1;
    para.y_bs_r = (m - (para.M_r-1)/2) * para.d;
end
