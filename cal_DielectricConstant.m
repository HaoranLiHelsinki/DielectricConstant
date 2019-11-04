function [ epsm ] = cal_DielectricConstant(Temperature, Frequence)
% Calculate dielectric constant as a function of temperature (Celsius) and Radar frequency 
% By: Haoran Li , haoran.li@helsinki.fi, 20191104
% Reference: Turner et al., JTECH, 2016

v = Frequence;%hz    
% 9.7 35.3 95.0

epsm_s = 8.7914e1  -4.0440e-1.*Temperature+9.5873e-4.*Temperature.^2-1.328e-6.*Temperature.^3;

a1 = 8.111e1;
b1 = 4.434e-3;
c1 = 1.302e-13;
d1 = 6.627e2;

a2 = 2.025;
b2 = 1.073e-2;
c2 = 1.012e-14;
d2 = 6.089e2;
Temperaturec = 1.342e2;

delTemperaturea1 = a1*exp(-b1*Temperature);
Temperatureou1 = c1*exp(d1/(Temperature+Temperaturec));
A1 = (Temperatureou1)^2*delTemperaturea1/(1+(2*pi*v*Temperatureou1)^2);

delTemperaturea2 = a2*exp(-b2*Temperature);
Temperatureou2 = c2*exp(d2/(Temperature+Temperaturec));
A2 = (Temperatureou2)^2*delTemperaturea2/(1+(2*pi*v*Temperatureou2)^2);

epsm_real = epsm_s -(2*pi*v)^2 * (A1+A2);

B1 = (Temperatureou1)*delTemperaturea1/(1+(2*pi*v*Temperatureou1)^2);

B2 = (Temperatureou2)*delTemperaturea2/(1+(2*pi*v*Temperatureou2)^2);


epsm_im = (2*pi*v) * (B1+B2);
epsm = complex(epsm_real, epsm_im);
    
KW2 = abs(( (epsm-1 )/(epsm+2 ) )^2);
 