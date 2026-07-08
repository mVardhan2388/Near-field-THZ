function para = update_positions(para)
    para.b = zeros(2,para.M_t);
    for i = 0:(para.M_t-1)
        para.b(2,i+1) = (i - (para.M_t - 1)/2)*para.d;
    end
end
