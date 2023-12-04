function [m_cold,m_h,m_fuel,m_cooling_s,m_cooling_r,m_h_corr] = massflow_new(m_a,FAR_real,BPR,split)
%MASSFLOW_NEW Summary of this function goes here
%   Detailed explanation goes here

alfa = 0.6;
m_cold = (m_a * BPR) / (BPR +1);
m_h = m_a/(BPR + 1);
m_cooling = split * m_h;
m_h_corr = (1-split) * m_h;
m_fuel = FAR_real * m_h_corr;
m_cooling_s = alfa * m_cooling;
m_cooling_r = (1 - alfa) * m_cooling;

end

