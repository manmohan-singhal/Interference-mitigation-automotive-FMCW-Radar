function rx_Power = interference_power(range,txAntennaGain,rxAntennaGain,txPower,centerfreq)

c = 3e8;
rx_Power = txPower*txAntennaGain*(c^2)*rxAntennaGain*(1/(((centerfreq)^2)*((4*pi*range)^4)));

