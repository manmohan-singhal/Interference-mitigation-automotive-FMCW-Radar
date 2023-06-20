function rx_Power = reflection_power(range,txAntennaGain,rxAntennaGain,txPower,centerfreq,RCS_object)

c = 3e8;
rx_Power = txPower*txAntennaGain*RCS_object*(c^2)*rxAntennaGain*(1/(((centerfreq)^2)*((4*pi)^3)*(range^4)));

