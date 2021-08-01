# File generated automatically from dogtail.f
# Equations of motion for model
function eom(u, p, t)
    @unpack ba, bb, eqx, eqy, eqz, g, ixa, ixb, iya, iyb, iza, izb, ka, kb, la, lao, lb, lbo, ma, mb, ox, oy, oz, px, py, pz, z = p 
    @inbounds q1, q2, q3, q4, q5, q6, u1, u2, u3, u5, u6 = u

    u4 = 0
	u4p = 0

    oxpp = derivative(ox, t, 2)
	oypp = derivative(oy, t, 2)
	ozpp = derivative(oz, t, 2)
	pxpp = derivative(px, t, 2)
	pypp = derivative(py, t, 2)
	pzpp = derivative(pz, t, 2)
    oxp = derivative(ox, t, 1)
    oyp = derivative(oy, t, 1)
    ozp = derivative(oz, t, 1)
    ox = ox(t)
    oy = oy(t)
    oz = oz(t)

    # Torques
    r1 = RotXYZ(q1 - eqx, q2 - eqy, q3 - eqz); r2 = RotXYZ(q4, q5, q6)
    atorx, atory, atorz = -ka * rotation_angle(r1) * rotation_axis(r1) .- ba * [u1,u2,u3]
    btorx, btory, btorz = -kb * rotation_angle(r2) * rotation_axis(r2) .- bb * [u4,u5,u6]

    z[16] = cos(q3)
	z[15] = cos(q2)
	z[46] = z[16]/z[15]
	z[18] = sin(q3)
	z[47] = z[18]/z[15]
	q1p = z[46]*u1 - z[47]*u2
	q2p = z[16]*u2 + z[18]*u1
	z[48] = tan(q2)
	z[50] = z[18]*z[48]
	z[49] = z[16]*z[48]
	q3p = u3 + z[50]*u2 - z[49]*u1
	z[30] = cos(q6)
	z[29] = cos(q5)
	z[51] = z[30]/z[29]
	z[32] = sin(q6)
	z[52] = z[32]/z[29]
	q4p = z[51]*u4 - z[52]*u5
	q5p = z[30]*u5 + z[32]*u4
	z[53] = tan(q5)
	z[55] = z[32]*z[53]
	z[54] = z[30]*z[53]
	q6p = u6 + z[55]*u5 - z[54]*u4
	z[17] = z[15]*z[16]
	z[19] = z[15]*z[18]
	z[20] = sin(q2)
	z[21] = cos(q1)
	z[22] = sin(q1)
	z[23] = z[18]*z[21] + z[16]*z[20]*z[22]
	z[24] = z[16]*z[21] - z[18]*z[20]*z[22]
	z[25] = z[15]*z[22]
	z[26] = z[18]*z[22] - z[16]*z[20]*z[21]
	z[27] = z[16]*z[22] + z[18]*z[20]*z[21]
	z[28] = z[15]*z[21]
	z[31] = z[29]*z[30]
	z[33] = z[29]*z[32]
	z[34] = sin(q5)
	z[35] = cos(q4)
	z[36] = sin(q4)
	z[37] = z[32]*z[35] + z[30]*z[34]*z[36]
	z[38] = z[30]*z[35] - z[32]*z[34]*z[36]
	z[39] = z[29]*z[36]
	z[40] = z[32]*z[36] - z[30]*z[34]*z[35]
	z[41] = z[30]*z[36] + z[32]*z[34]*z[35]
	z[42] = z[29]*z[35]
	z[60] = z[17]*z[31] + z[20]*z[40] - z[19]*z[37]
	z[61] = z[20]*z[41] - z[17]*z[33] - z[19]*z[38]
	z[62] = z[17]*z[34] + z[19]*z[39] + z[20]*z[42]
	z[63] = z[23]*z[31] + z[24]*z[37] - z[25]*z[40]
	z[64] = z[24]*z[38] - z[23]*z[33] - z[25]*z[41]
	z[65] = z[23]*z[34] - z[24]*z[39] - z[25]*z[42]
	z[66] = z[26]*z[31] + z[27]*z[37] + z[28]*z[40]
	z[67] = z[27]*z[38] + z[28]*z[41] - z[26]*z[33]
	z[68] = z[26]*z[34] + z[28]*z[42] - z[27]*z[39]
	z[69] = lbo*z[39]
	z[70] = lbo*z[34]
	z[71] = lbo*z[42]
	z[73] = lbo*z[33]
	z[74] = lbo*z[38]
	z[75] = lbo*z[41]
	z[94] = z[15]*z[21]*q1p - z[20]*z[22]*q2p
	z[95] = -z[15]*z[22]*q1p - z[20]*z[21]*q2p
	z[97] = z[15]*z[16]*q3p - z[18]*z[20]*q2p
	z[98] = -z[16]*z[22]*q1p - z[18]*z[21]*q3p - z[15]*z[18]*z[22]*q2p - z[16]*z[20]*z[22]*q3p - z[18]*z[20]*z[21]*q1p
	z[99] = z[16]*z[21]*q1p + z[15]*z[18]*z[21]*q2p + z[16]*z[20]*z[21]*q3p - z[18]*z[22]*q3p - z[18]*z[20]*z[22]*q1p
	z[113] = -z[15]*z[18]*q3p - z[16]*z[20]*q2p
	z[114] = z[29]*z[35]*q4p - z[34]*z[36]*q5p
	z[115] = -z[29]*z[36]*q4p - z[34]*z[35]*q5p
	z[116] = z[15]*z[42]*q2p + z[17]*z[29]*q5p + z[19]*z[114] + z[20]*z[115] + z[34]*z[113] + z[39]*z[97]
	z[117] = z[16]*z[21]*q3p + z[15]*z[16]*z[22]*q2p + z[16]*z[20]*z[21]*q1p - z[18]*z[22]*q1p - z[18]*z[20]*z[22]*q3p
	z[118] = z[23]*z[29]*q5p + z[34]*z[117] - z[24]*z[114] - z[25]*z[115] - z[39]*z[98] - z[42]*z[94]
	z[119] = z[16]*z[22]*q3p + z[18]*z[21]*q1p + z[16]*z[20]*z[22]*q1p + z[18]*z[20]*z[21]*q3p - z[15]*z[16]*z[21]*q2p
	z[120] = z[26]*z[29]*q5p + z[28]*z[115] + z[34]*z[119] + z[42]*z[95] - z[27]*z[114] - z[39]*z[99]
	z[122] = lbo*z[29]*q5p
	z[124] = z[29]*z[30]*q6p - z[32]*z[34]*q5p
	z[125] = z[30]*z[35]*q4p + z[29]*z[32]*z[35]*q5p + z[30]*z[34]*z[35]*q6p - z[32]*z[36]*q6p - z[32]*z[34]*z[36]*q4p
	z[126] = -z[30]*z[36]*q4p - z[32]*z[35]*q6p - z[29]*z[32]*z[36]*q5p - z[30]*z[34]*z[36]*q6p - z[32]*z[34]*z[35]*q4p
	z[127] = z[15]*z[41]*q2p + z[20]*z[125] - z[17]*z[124] - z[19]*z[126] - z[33]*z[113] - z[38]*z[97]
	z[128] = z[24]*z[126] + z[38]*z[98] - z[23]*z[124] - z[25]*z[125] - z[33]*z[117] - z[41]*z[94]
	z[129] = z[27]*z[126] + z[28]*z[125] + z[38]*z[99] + z[41]*z[95] - z[26]*z[124] - z[33]*z[119]
	z[162] = ixa*u1
	z[166] = iya*u2
	z[170] = iza*u3
	z[181] = ixb*z[31]*u1
	z[182] = ixb*z[37]*u2
	z[183] = ixb*z[40]*u3
	z[187] = iyb*u5
	z[188] = iyb*z[33]*u1
	z[189] = iyb*z[38]*u2
	z[190] = iyb*z[41]*u3
	z[194] = izb*u6
	z[195] = izb*z[34]*u1
	z[196] = izb*z[39]*u2
	z[197] = izb*z[42]*u3
	z[201] = ixb*z[31]
	z[202] = ixb*z[37]
	z[203] = ixb*z[40]
	z[205] = iyb*z[38]
	z[206] = iyb*z[41]
	z[207] = iyb*z[33]
	z[209] = izb*z[34]
	z[210] = izb*z[42]
	z[211] = izb*z[39]
	z[216] = ixa + z[31]*z[201] + z[33]*z[207] + z[34]*z[209] + mb*(z[70]^2+z[73]^2)
	z[217] = z[31]*z[202] - z[33]*z[205] - z[34]*z[211] - mb*(z[69]*z[70]+z[73]*z[74]+la*z[41]*z[70]+la*z[42]*z[73])
	z[218] = z[31]*z[203] + z[34]*z[210] + mb*(z[70]*z[71]+la*z[38]*z[70]-z[73]*z[75]-la*z[39]*z[73]) - z[33]*z[206]
	z[220] = izb*z[34] + z[219]*z[70]
	z[221] = -iyb*z[33] - z[219]*z[73]
	z[225] = z[223] + z[37]*z[202] + z[38]*z[205] + z[39]*z[211] + mb*(z[224]+z[69]^2+z[74]^2+2*la*z[41]*z[69]+2*la*z[42]*z[74])
	z[226] = z[37]*z[201] - z[38]*z[207] - z[39]*z[209] - mb*(z[69]*z[70]+z[73]*z[74]+la*z[41]*z[70]+la*z[42]*z[73])
	z[227] = z[37]*z[203] + z[38]*z[206] - z[39]*z[210] - mb*(z[69]*z[71]+la*z[38]*z[69]+la*z[41]*z[71]-z[74]*z[75]-la*z[39]*z[74]-la*z[42]*z[75])
	z[228] = iyb*z[38] + z[219]*(z[74]+la*z[42])
	z[229] = -izb*z[39] - z[219]*(z[69]+la*z[41])
	z[233] = z[232] + z[40]*z[203] + z[41]*z[206] + z[42]*z[210] + mb*(z[224]+z[71]^2+z[75]^2+2*la*z[38]*z[71]+2*la*z[39]*z[75])
	z[234] = z[40]*z[201] + z[42]*z[209] + mb*(z[70]*z[71]+la*z[38]*z[70]-z[73]*z[75]-la*z[39]*z[73]) - z[41]*z[207]
	z[235] = z[40]*z[202] + z[41]*z[205] - z[42]*z[211] - mb*(z[69]*z[71]+la*z[38]*z[69]+la*z[41]*z[71]-z[74]*z[75]-la*z[39]*z[74]-la*z[42]*z[75])
	z[236] = iyb*z[41] + z[219]*(z[75]+la*z[39])
	z[237] = izb*z[42] + z[219]*(z[71]+la*z[38])
	z[240] = z[205] + z[219]*(z[74]+la*z[42])
	z[241] = z[206] + z[219]*(z[75]+la*z[39])
	z[242] = -z[207] - z[219]*z[73]
	z[245] = z[209] + z[219]*z[70]
	z[246] = z[210] + z[219]*(z[71]+la*z[38])
	z[247] = -z[211] - z[219]*(z[69]+la*z[41])
	z[1] = cos(oy)
	z[2] = cos(oz)
	z[3] = z[1]*z[2]
	z[4] = sin(oz)
	z[5] = z[1]*z[4]
	z[6] = sin(oy)
	z[7] = cos(ox)
	z[8] = sin(ox)
	z[9] = z[4]*z[7] + z[2]*z[6]*z[8]
	z[10] = z[2]*z[7] - z[4]*z[6]*z[8]
	z[11] = z[1]*z[8]
	z[12] = z[4]*z[8] - z[2]*z[6]*z[7]
	z[13] = z[2]*z[8] + z[4]*z[6]*z[7]
	z[14] = z[1]*z[7]
	z[43] = z[4]*oyp + z[1]*z[2]*oxp
	z[44] = z[2]*oyp - z[1]*z[4]*oxp
	z[45] = ozp + z[6]*oxp
	z[56] = lao*(z[44]*z[25]-z[43]*z[20]-z[45]*z[28])
	z[57] = lao*(z[43]*z[19]-z[44]*z[24]-z[45]*z[27])
	z[58] = la*(z[44]*z[25]-z[43]*z[20]-z[45]*z[28])
	z[59] = la*(z[43]*z[19]-z[44]*z[24]-z[45]*z[27])
	z[72] = lbo*(z[43]*z[62]+z[44]*z[65]+z[45]*z[68])
	z[76] = lbo*(z[43]*z[61]+z[44]*z[64]+z[45]*z[67])
	z[77] = z[2]*oyp*ozp + z[4]*oypp + z[1]*z[2]*oxpp - z[1]*z[4]*oxp*ozp - z[2]*z[6]*oxp*oyp
	z[78] = z[4]*z[6]*oxp*oyp + z[2]*oypp - z[4]*oyp*ozp - z[1]*z[2]*oxp*ozp - z[1]*z[4]*oxpp
	z[79] = z[1]*oxp*oyp + ozpp + z[6]*oxpp
	z[80] = u2*(z[44]*z[25]-z[43]*z[20]-z[45]*z[28]-u3) - u3*(z[43]*z[19]-z[44]*z[24]-z[45]*z[27]-u2)
	z[81] = -u3*(z[43]*z[17]+z[44]*z[23]+z[45]*z[26]+u1) - u1*(z[44]*z[25]-z[43]*z[20]-z[45]*z[28]-u3)
	z[82] = u2*(z[43]*z[17]+z[44]*z[23]+z[45]*z[26]+u1) + u1*(z[43]*z[19]-z[44]*z[24]-z[45]*z[27]-u2)
	z[83] = z[77]*z[17] + z[78]*z[23] + z[79]*z[26] + z[80]
	z[84] = z[78]*z[24] + z[79]*z[27] + z[81] - z[77]*z[19]
	z[85] = z[77]*z[20] + z[79]*z[28] + z[82] - z[78]*z[25]
	z[86] = u5*(z[39]*u2-z[43]*z[62]-z[44]*z[65]-z[45]*z[68]-u6-z[34]*u1-z[42]*u3) - u6*(z[33]*u1-z[43]*z[61]-z[44]*z[64]-z[45]*z[67]-u5-z[38]*u2-z[41]*u3)
	z[87] = -u6*(z[43]*z[60]+z[44]*z[63]+z[45]*z[66]+u4+z[31]*u1+z[37]*u2+z[40]*u3) - u4*(z[39]*u2-z[43]*z[62]-z[44]*z[65]-z[45]*z[68]-u6-z[34]*u1-z[42]*u3)
	z[88] = u5*(z[43]*z[60]+z[44]*z[63]+z[45]*z[66]+u4+z[31]*u1+z[37]*u2+z[40]*u3) + u4*(z[33]*u1-z[43]*z[61]-z[44]*z[64]-z[45]*z[67]-u5-z[38]*u2-z[41]*u3)
	z[89] = z[77]*z[60] + z[78]*z[63] + z[79]*z[66] + z[86] + z[31]*z[80] + z[37]*z[81] + z[40]*z[82]
	z[90] = z[77]*z[61] + z[78]*z[64] + z[79]*z[67] + z[87] + z[38]*z[81] + z[41]*z[82] - z[33]*z[80]
	z[91] = z[77]*z[62] + z[78]*z[65] + z[79]*z[68] + z[88] + z[34]*z[80] + z[42]*z[82] - z[39]*z[81]
	z[92] = lao*u3 - z[56]
	z[93] = z[57] - lao*u2
	z[96] = lao*(z[78]*z[25]+z[44]*z[94]-z[77]*z[20]-z[79]*z[28]-z[43]*z[15]*q2p-z[45]*z[95])
	z[100] = lao*(z[77]*z[19]+z[43]*z[97]-z[78]*z[24]-z[79]*z[27]-z[44]*z[98]-z[45]*z[99])
	z[102] = -(z[43]*z[17]+z[44]*z[23]+z[45]*z[26]+u1)*z[93] - z[96]
	z[103] = (z[43]*z[17]+z[44]*z[23]+z[45]*z[26]+u1)*z[92] + z[100]
	z[104] = la*u3 - z[58]
	z[105] = z[59] - la*u2
	z[106] = la*(z[78]*z[25]+z[44]*z[94]-z[77]*z[20]-z[79]*z[28]-z[43]*z[15]*q2p-z[45]*z[95])
	z[107] = la*(z[77]*z[19]+z[43]*z[97]-z[78]*z[24]-z[79]*z[27]-z[44]*z[98]-z[45]*z[99])
	z[108] = (z[44]*z[25]-z[43]*z[20]-z[45]*z[28]-u3)*z[104] - (z[43]*z[19]-z[44]*z[24]-z[45]*z[27]-u2)*z[105]
	z[109] = -(z[43]*z[17]+z[44]*z[23]+z[45]*z[26]+u1)*z[105] - z[106]
	z[110] = (z[43]*z[17]+z[44]*z[23]+z[45]*z[26]+u1)*z[104] + z[107]
	z[111] = z[72] + lbo*u6 + z[70]*u1 + z[71]*u3 - z[69]*u2
	z[112] = z[73]*u1 - z[76] - lbo*u5 - z[74]*u2 - z[75]*u3
	z[121] = lbo*(z[77]*z[62]+z[78]*z[65]+z[79]*z[68]+z[43]*z[116]+z[44]*z[118]+z[45]*z[120])
	z[123] = z[121] + u1*z[122] + lbo*u3*z[115] - lbo*u2*z[114]
	z[130] = lbo*(z[77]*z[61]+z[78]*z[64]+z[79]*z[67]+z[43]*z[127]+z[44]*z[128]+z[45]*z[129])
	z[131] = lbo*u1*z[124] - z[130] - lbo*u2*z[126] - lbo*u3*z[125]
	z[132] = (z[39]*u2-z[43]*z[62]-z[44]*z[65]-z[45]*z[68]-u6-z[34]*u1-z[42]*u3)*z[111] - (z[33]*u1-z[43]*z[61]-z[44]*z[64]-z[45]*z[67]-u5-z[38]*u2-z[41]*u3)*z[112]
	z[133] = z[123] - (z[43]*z[60]+z[44]*z[63]+z[45]*z[66]+u4+z[31]*u1+z[37]*u2+z[40]*u3)*z[112]
	z[134] = (z[43]*z[60]+z[44]*z[63]+z[45]*z[66]+u4+z[31]*u1+z[37]*u2+z[40]*u3)*z[111] + z[131]
	z[137] = z[3]*z[17] + z[6]*z[26] - z[5]*z[23]
	z[138] = z[9]*z[17] + z[10]*z[23] - z[11]*z[26]
	z[139] = z[12]*z[17] + z[13]*z[23] + z[14]*z[26]
	z[140] = z[6]*z[27] - z[3]*z[19] - z[5]*z[24]
	z[141] = z[10]*z[24] - z[9]*z[19] - z[11]*z[27]
	z[142] = z[13]*z[24] + z[14]*z[27] - z[12]*z[19]
	z[143] = z[3]*z[20] + z[5]*z[25] + z[6]*z[28]
	z[144] = z[9]*z[20] - z[10]*z[25] - z[11]*z[28]
	z[145] = z[12]*z[20] + z[14]*z[28] - z[13]*z[25]
	z[149] = z[38]*z[140] + z[41]*z[143] - z[33]*z[137]
	z[150] = z[38]*z[141] + z[41]*z[144] - z[33]*z[138]
	z[151] = z[38]*z[142] + z[41]*z[145] - z[33]*z[139]
	z[152] = z[34]*z[137] + z[42]*z[143] - z[39]*z[140]
	z[153] = z[34]*z[138] + z[42]*z[144] - z[39]*z[141]
	z[154] = z[34]*z[139] + z[42]*z[145] - z[39]*z[142]
	z[155] = atorx + z[136]*(z[70]*z[151]+z[73]*z[154])
	z[157] = atory - z[156]*z[145] - z[136]*(la*z[145]+z[69]*z[151]+z[74]*z[154])
	z[158] = atorz + z[156]*z[142] + z[136]*(la*z[142]+z[71]*z[151]-z[75]*z[154])
	z[160] = btory - z[159]*z[154]
	z[161] = btorz + z[159]*z[151]
	z[163] = ixa*z[43]*z[17]
	z[164] = ixa*z[44]*z[23]
	z[165] = ixa*z[45]*z[26]
	z[167] = iya*z[43]*z[19]
	z[168] = iya*z[44]*z[24]
	z[169] = iya*z[45]*z[27]
	z[171] = iza*z[43]*z[20]
	z[172] = iza*z[44]*z[25]
	z[173] = iza*z[45]*z[28]
	z[174] = ixa*z[83]
	z[175] = iya*z[84]
	z[176] = iza*z[85]
	z[177] = z[43]*z[17]*z[168] + z[43]*z[17]*z[169] + z[43]*z[19]*z[163] + z[43]*z[19]*z[164] + z[43]*z[19]*z[165] + z[44]*z[23]*z[168] + z[44]*z[23]*z[169] + z[45]*z[26]*z[168] + z[45]*z[26]*z[169] + z[168]*u1 + z[169]*u1 + u1*z[166] + z[43]*z[17]*z[166] + z[43]*z[19]*z[162] + z[44]*z[23]*z[166] + z[45]*z[26]*z[166] - z[43]*z[17]*z[167] - z[44]*z[23]*z[167] - z[44]*z[24]*z[163] - z[44]*z[24]*z[164] - z[44]*z[24]*z[165] - z[45]*z[26]*z[167] - z[45]*z[27]*z[163] - z[45]*z[27]*z[164] - z[45]*z[27]*z[165] - z[163]*u2 - z[164]*u2 - z[165]*u2 - z[167]*u1 - u2*z[162] - z[44]*z[24]*z[162] - z[45]*z[27]*z[162]
	z[178] = z[43]*z[17]*z[172] + z[43]*z[20]*z[163] + z[43]*z[20]*z[164] + z[43]*z[20]*z[165] + z[44]*z[23]*z[172] + z[45]*z[26]*z[172] + z[45]*z[28]*z[163] + z[45]*z[28]*z[164] + z[45]*z[28]*z[165] + z[163]*u3 + z[164]*u3 + z[165]*u3 + z[172]*u1 + u3*z[162] + z[43]*z[20]*z[162] + z[45]*z[28]*z[162] - z[43]*z[17]*z[171] - z[43]*z[17]*z[173] - z[44]*z[23]*z[171] - z[44]*z[23]*z[173] - z[44]*z[25]*z[163] - z[44]*z[25]*z[164] - z[44]*z[25]*z[165] - z[45]*z[26]*z[171] - z[45]*z[26]*z[173] - z[171]*u1 - z[173]*u1 - u1*z[170] - z[43]*z[17]*z[170] - z[44]*z[23]*z[170] - z[44]*z[25]*z[162] - z[45]*z[26]*z[170]
	z[179] = z[43]*z[19]*z[172] + z[43]*z[20]*z[167] + z[44]*z[24]*z[171] + z[44]*z[24]*z[173] + z[44]*z[25]*z[168] + z[44]*z[25]*z[169] + z[45]*z[27]*z[171] + z[45]*z[27]*z[173] + z[45]*z[28]*z[167] + z[167]*u3 + z[171]*u2 + z[173]*u2 + u2*z[170] + z[44]*z[24]*z[170] + z[44]*z[25]*z[166] + z[45]*z[27]*z[170] - z[43]*z[19]*z[171] - z[43]*z[19]*z[173] - z[43]*z[20]*z[168] - z[43]*z[20]*z[169] - z[44]*z[24]*z[172] - z[44]*z[25]*z[167] - z[45]*z[27]*z[172] - z[45]*z[28]*z[168] - z[45]*z[28]*z[169] - z[168]*u3 - z[169]*u3 - z[172]*u2 - u3*z[166] - z[43]*z[19]*z[170] - z[43]*z[20]*z[166] - z[45]*z[28]*z[166]
	z[184] = ixb*z[43]*z[60]
	z[185] = ixb*z[44]*z[63]
	z[186] = ixb*z[45]*z[66]
	z[191] = iyb*z[43]*z[61]
	z[192] = iyb*z[44]*z[64]
	z[193] = iyb*z[45]*z[67]
	z[198] = izb*z[43]*z[62]
	z[199] = izb*z[44]*z[65]
	z[200] = izb*z[45]*z[68]
	z[204] = ixb*z[89]
	z[208] = iyb*z[90]
	z[212] = izb*z[91]
	z[213] = z[43]*z[60]*z[191] + z[43]*z[60]*z[192] + z[43]*z[60]*z[193] + z[44]*z[63]*z[191] + z[44]*z[63]*z[192] + z[44]*z[63]*z[193] + z[45]*z[66]*z[191] + z[45]*z[66]*z[192] + z[45]*z[66]*z[193] + z[191]*u4 + z[192]*u4 + z[193]*u4 + z[31]*z[191]*u1 + z[31]*z[192]*u1 + z[31]*z[193]*u1 + z[33]*z[184]*u1 + z[33]*z[185]*u1 + z[33]*z[186]*u1 + z[37]*z[191]*u2 + z[37]*z[192]*u2 + z[37]*z[193]*u2 + z[40]*z[191]*u3 + z[40]*z[192]*u3 + z[40]*z[193]*u3 + u4*z[187] + u4*z[189] + u4*z[190] + z[43]*z[60]*z[187] + z[43]*z[60]*z[189] + z[43]*z[60]*z[190] + z[44]*z[63]*z[187] + z[44]*z[63]*z[189] + z[44]*z[63]*z[190] + z[45]*z[66]*z[187] + z[45]*z[66]*z[189] + z[45]*z[66]*z[190] + z[31]*u1*z[187] + z[31]*u1*z[189] + z[31]*u1*z[190] + z[33]*u1*z[180] + z[33]*u1*z[181] + z[33]*u1*z[182] + z[33]*u1*z[183] + z[37]*u2*z[187] + z[37]*u2*z[189] + z[37]*u2*z[190] + z[40]*u3*z[187] + z[40]*u3*z[189] + z[40]*u3*z[190] - z[43]*z[61]*z[184] - z[43]*z[61]*z[185] - z[43]*z[61]*z[186] - z[44]*z[64]*z[184] - z[44]*z[64]*z[185] - z[44]*z[64]*z[186] - z[45]*z[67]*z[184] - z[45]*z[67]*z[185] - z[45]*z[67]*z[186] - z[184]*u5 - z[185]*u5 - z[186]*u5 - z[38]*z[184]*u2 - z[38]*z[185]*u2 - z[38]*z[186]*u2 - z[41]*z[184]*u3 - z[41]*z[185]*u3 - z[41]*z[186]*u3 - u4*z[188] - u5*z[180] - u5*z[181] - u5*z[182] - u5*z[183] - z[43]*z[60]*z[188] - z[43]*z[61]*z[180] - z[43]*z[61]*z[181] - z[43]*z[61]*z[182] - z[43]*z[61]*z[183] - z[44]*z[63]*z[188] - z[44]*z[64]*z[180] - z[44]*z[64]*z[181] - z[44]*z[64]*z[182] - z[44]*z[64]*z[183] - z[45]*z[66]*z[188] - z[45]*z[67]*z[180] - z[45]*z[67]*z[181] - z[45]*z[67]*z[182] - z[45]*z[67]*z[183] - z[31]*u1*z[188] - z[37]*u2*z[188] - z[38]*u2*z[180] - z[38]*u2*z[181] - z[38]*u2*z[182] - z[38]*u2*z[183] - z[40]*u3*z[188] - z[41]*u3*z[180] - z[41]*u3*z[181] - z[41]*u3*z[182] - z[41]*u3*z[183]
	z[214] = z[43]*z[62]*z[184] + z[43]*z[62]*z[185] + z[43]*z[62]*z[186] + z[44]*z[65]*z[184] + z[44]*z[65]*z[185] + z[44]*z[65]*z[186] + z[45]*z[68]*z[184] + z[45]*z[68]*z[185] + z[45]*z[68]*z[186] + z[184]*u6 + z[185]*u6 + z[186]*u6 + z[34]*z[184]*u1 + z[34]*z[185]*u1 + z[34]*z[186]*u1 + z[42]*z[184]*u3 + z[42]*z[185]*u3 + z[42]*z[186]*u3 + u4*z[196] + u6*z[180] + u6*z[181] + u6*z[182] + u6*z[183] + z[43]*z[60]*z[196] + z[43]*z[62]*z[180] + z[43]*z[62]*z[181] + z[43]*z[62]*z[182] + z[43]*z[62]*z[183] + z[44]*z[63]*z[196] + z[44]*z[65]*z[180] + z[44]*z[65]*z[181] + z[44]*z[65]*z[182] + z[44]*z[65]*z[183] + z[45]*z[66]*z[196] + z[45]*z[68]*z[180] + z[45]*z[68]*z[181] + z[45]*z[68]*z[182] + z[45]*z[68]*z[183] + z[31]*u1*z[196] + z[34]*u1*z[180] + z[34]*u1*z[181] + z[34]*u1*z[182] + z[34]*u1*z[183] + z[37]*u2*z[196] + z[40]*u3*z[196] + z[42]*u3*z[180] + z[42]*u3*z[181] + z[42]*u3*z[182] + z[42]*u3*z[183] - z[43]*z[60]*z[198] - z[43]*z[60]*z[199] - z[43]*z[60]*z[200] - z[44]*z[63]*z[198] - z[44]*z[63]*z[199] - z[44]*z[63]*z[200] - z[45]*z[66]*z[198] - z[45]*z[66]*z[199] - z[45]*z[66]*z[200] - z[198]*u4 - z[199]*u4 - z[200]*u4 - z[31]*z[198]*u1 - z[31]*z[199]*u1 - z[31]*z[200]*u1 - z[37]*z[198]*u2 - z[37]*z[199]*u2 - z[37]*z[200]*u2 - z[39]*z[184]*u2 - z[39]*z[185]*u2 - z[39]*z[186]*u2 - z[40]*z[198]*u3 - z[40]*z[199]*u3 - z[40]*z[200]*u3 - u4*z[194] - u4*z[195] - u4*z[197] - z[43]*z[60]*z[194] - z[43]*z[60]*z[195] - z[43]*z[60]*z[197] - z[44]*z[63]*z[194] - z[44]*z[63]*z[195] - z[44]*z[63]*z[197] - z[45]*z[66]*z[194] - z[45]*z[66]*z[195] - z[45]*z[66]*z[197] - z[31]*u1*z[194] - z[31]*u1*z[195] - z[31]*u1*z[197] - z[37]*u2*z[194] - z[37]*u2*z[195] - z[37]*u2*z[197] - z[39]*u2*z[180] - z[39]*u2*z[181] - z[39]*u2*z[182] - z[39]*u2*z[183] - z[40]*u3*z[194] - z[40]*u3*z[195] - z[40]*u3*z[197]
	z[215] = z[43]*z[61]*z[198] + z[43]*z[61]*z[199] + z[43]*z[61]*z[200] + z[44]*z[64]*z[198] + z[44]*z[64]*z[199] + z[44]*z[64]*z[200] + z[45]*z[67]*z[198] + z[45]*z[67]*z[199] + z[45]*z[67]*z[200] + z[198]*u5 + z[199]*u5 + z[200]*u5 + z[38]*z[198]*u2 + z[38]*z[199]*u2 + z[38]*z[200]*u2 + z[39]*z[191]*u2 + z[39]*z[192]*u2 + z[39]*z[193]*u2 + z[41]*z[198]*u3 + z[41]*z[199]*u3 + z[41]*z[200]*u3 + u5*z[194] + u5*z[195] + u5*z[197] + u6*z[188] + z[43]*z[61]*z[194] + z[43]*z[61]*z[195] + z[43]*z[61]*z[197] + z[43]*z[62]*z[188] + z[44]*z[64]*z[194] + z[44]*z[64]*z[195] + z[44]*z[64]*z[197] + z[44]*z[65]*z[188] + z[45]*z[67]*z[194] + z[45]*z[67]*z[195] + z[45]*z[67]*z[197] + z[45]*z[68]*z[188] + z[33]*u1*z[196] + z[34]*u1*z[188] + z[38]*u2*z[194] + z[38]*u2*z[195] + z[38]*u2*z[197] + z[39]*u2*z[187] + z[39]*u2*z[189] + z[39]*u2*z[190] + z[41]*u3*z[194] + z[41]*u3*z[195] + z[41]*u3*z[197] + z[42]*u3*z[188] - z[43]*z[62]*z[191] - z[43]*z[62]*z[192] - z[43]*z[62]*z[193] - z[44]*z[65]*z[191] - z[44]*z[65]*z[192] - z[44]*z[65]*z[193] - z[45]*z[68]*z[191] - z[45]*z[68]*z[192] - z[45]*z[68]*z[193] - z[191]*u6 - z[192]*u6 - z[193]*u6 - z[33]*z[198]*u1 - z[33]*z[199]*u1 - z[33]*z[200]*u1 - z[34]*z[191]*u1 - z[34]*z[192]*u1 - z[34]*z[193]*u1 - z[42]*z[191]*u3 - z[42]*z[192]*u3 - z[42]*z[193]*u3 - u5*z[196] - u6*z[187] - u6*z[189] - u6*z[190] - z[43]*z[61]*z[196] - z[43]*z[62]*z[187] - z[43]*z[62]*z[189] - z[43]*z[62]*z[190] - z[44]*z[64]*z[196] - z[44]*z[65]*z[187] - z[44]*z[65]*z[189] - z[44]*z[65]*z[190] - z[45]*z[67]*z[196] - z[45]*z[68]*z[187] - z[45]*z[68]*z[189] - z[45]*z[68]*z[190] - z[33]*u1*z[194] - z[33]*u1*z[195] - z[33]*u1*z[197] - z[34]*u1*z[187] - z[34]*u1*z[189] - z[34]*u1*z[190] - z[38]*u2*z[196] - z[39]*u2*z[188] - z[41]*u3*z[196] - z[42]*u3*z[187] - z[42]*u3*z[189] - z[42]*u3*z[190]
	z[222] = z[174] + z[179] + z[31]*z[204] + z[31]*z[215] + z[34]*z[212] + z[34]*z[213] - z[33]*z[208] - z[33]*z[214] - mb*(z[33]*z[70]*z[108]+z[39]*z[73]*z[109]-pxpp*z[70]*z[149]-pxpp*z[73]*z[152]-pypp*z[70]*z[150]-pypp*z[73]*z[153]-pzpp*z[70]*z[151]-pzpp*z[73]*z[154]-z[34]*z[73]*z[108]-z[70]*z[133]-z[73]*z[134]-z[38]*z[70]*z[109]-z[41]*z[70]*z[110]-z[42]*z[73]*z[110])
	z[231] = z[175] + z[178] + z[37]*z[204] + z[37]*z[215] + z[38]*z[208] + z[38]*z[214] + mb*(z[33]*z[69]*z[108]+z[39]*z[74]*z[109]-la*pxpp*z[143]-la*pypp*z[144]-la*pzpp*z[145]-pxpp*z[69]*z[149]-pxpp*z[74]*z[152]-pypp*z[69]*z[150]-pypp*z[74]*z[153]-pzpp*z[69]*z[151]-pzpp*z[74]*z[154]-la*z[40]*z[132]-z[34]*z[74]*z[108]-la*z[110]-z[69]*z[133]-z[74]*z[134]-la*z[41]*z[133]-la*z[42]*z[134]-z[38]*z[69]*z[109]-z[41]*z[69]*z[110]-z[42]*z[74]*z[110]) - z[39]*z[212] - z[39]*z[213] - z[230]*(pxpp*z[143]+pypp*z[144]+pzpp*z[145]+z[103])
	z[238] = z[176] + z[177] + z[40]*z[204] + z[40]*z[215] + z[41]*z[208] + z[41]*z[214] + z[42]*z[212] + z[42]*z[213] + z[230]*(pxpp*z[140]+pypp*z[141]+pzpp*z[142]+z[102]) + mb*(la*pxpp*z[140]+la*pypp*z[141]+la*pzpp*z[142]+pxpp*z[71]*z[149]+pypp*z[71]*z[150]+pzpp*z[71]*z[151]+la*z[37]*z[132]+la*z[109]+z[71]*z[133]+la*z[38]*z[133]+z[38]*z[71]*z[109]+z[39]*z[75]*z[109]+z[41]*z[71]*z[110]-pxpp*z[75]*z[152]-pypp*z[75]*z[153]-pzpp*z[75]*z[154]-z[33]*z[71]*z[108]-z[34]*z[75]*z[108]-z[75]*z[134]-la*z[39]*z[134]-z[42]*z[75]*z[110])
	z[243] = z[208] + z[214] + z[219]*(z[39]*z[109]-pxpp*z[152]-pypp*z[153]-pzpp*z[154]-z[34]*z[108]-z[134]-z[42]*z[110])
	z[248] = z[212] + z[213] - z[219]*(z[33]*z[108]-pxpp*z[149]-pypp*z[150]-pzpp*z[151]-z[133]-z[38]*z[109]-z[41]*z[110])
	z[250] = z[155] - z[222] - z[201]*u4p
	z[251] = z[157] - z[231] - z[202]*u4p
	z[252] = z[158] - z[238] - z[203]*u4p
	z[253] = z[160] - z[243]
	z[254] = z[161] - z[248]

    coef11 = -z[216]
	coef12 = -z[217]
	coef13 = -z[218]
	coef14 = -z[221]
	coef15 = -z[220]
	coef21 = -z[226]
	coef22 = -z[225]
	coef23 = -z[227]
	coef24 = -z[228]
	coef25 = -z[229]
	coef31 = -z[234]
	coef32 = -z[235]
	coef33 = -z[233]
	coef34 = -z[236]
	coef35 = -z[237]
	coef41 = -z[242]
	coef42 = -z[240]
	coef43 = -z[241]
	coef44 = -z[239]
	coef45 = 0
	coef51 = -z[245]
	coef52 = -z[247]
	coef53 = -z[246]
	coef54 = 0
	coef55 = -z[244]
    rhs1 = -z[250]
	rhs2 = -z[251]
	rhs3 = -z[252]
	rhs4 = -z[253]
	rhs5 = -z[254]
    
    coef = @SMatrix [coef11 coef12 coef13 coef14 coef15
					 coef21 coef22 coef23 coef24 coef25
					 coef31 coef32 coef33 coef34 coef35
					 coef41 coef42 coef43 coef44 coef45
					 coef51 coef52 coef53 coef54 coef55]

    rhs = @SVector [rhs1, rhs2, rhs3, rhs4, rhs5]

    @inbounds u1p, u2p, u3p, u5p, u6p = coef \ rhs

    @SVector [q1p, q2p, q3p, q4p, q5p, q6p, u1p, u2p, u3p, u5p, u6p]
end  


# Create model
function createModel(fname)
    
    # Parameters
    px, py, pz, ox, oy, oz,  time, base, mid, tip, initconds, la, lb, tspan = getdata(fname)
	deleteat!(initconds, 10) # Remove u4 
    p = initialise_parameters(px, py, pz, ox, oy, oz, la, lb)
	@unpack g, ma, mb, lao, lbo, ixb, iya, iyb, iza, izb, z = p

    # Initial conditions
    u₀ = SVector{11,Float64}(initconds)

    prob = ODEProblem(eom, u₀, tspan, p)

    return prob, p, time, base, mid, tip

end
