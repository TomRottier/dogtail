# Fixed parameters

function parameters(;fx, fy, fz, ka, kb, ba, bb, la, lb)
    ma = 1;
    mb = 1;
    # la = lb = 0.5;
    lao = 0.5 * la; lbo = 0.5 * lb
    ixa = ixb = 0.01;
    iya = iza = 1 / 12 * ma * lao;
    iyb = izb = 1 / 12 * mb * lbo
    g = -9.81
    # ka = 1; 
    # kb = 10;
    # ba = bb = 0.1

    return (ma, mb, la, lb, lao, lbo, ixa, ixb, iya, iyb, iza, izb, g, ka, kb, ba, bb, fx = fx, fy = fy, fz = fz)
end
