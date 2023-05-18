clearvars
close all


f = figuresetting('centimeters',21,15,'Helvetica',24,'off',2,'off','off');
% g = figurelayout('langmuir');
% g2 = figurelayout('langmuir2');
% fprl0p5C = figurelayout('c0p5_prl');
% fprl1C = figurelayout('c1_prl');
% fprl1p5C = figurelayout('c1p5_prl');
% fprl2C = figurelayout('c2_prl');
% fprl1p5R = figurelayout('r1p5_prl');

constant=struct('Kb',1.380649*10^(-23),'e',1.602176*10^(-19),...
    'Na',6.02214076*10^23,'epsilon0',8.8541878128*10^(-12),'c',299792458,...
    'ThermalEnergy',4.11*10^(-21),'visc_water_23K',0.00092353,'Plank',6.62607004*10^(-34));

SI=struct('meter',1,'second',1,'minute',60,'hour',3600,'newton',1,'joule',1,'kg',1,'g',10^-3,'liter',10^-3,'hertz',1,...
    'kelvin',1,'watt',1,'coulomb',1,'volt',1,'farad',1,'ohm',1,'pascal',1,'bar',10^5,...
    'angstrom',10^-10,'eV',1.602176565*10^-19,'zepto',10^-21,'atto',10^-18,'femto',10^-15,...
    'pico',10^-12,'nano',10^-9,'micro',10^-6,'milli',10^-3,'centi',10^-2,'deci',10^-1,'kilo',10^3,'mega',10^6,'giga',10^9,...
'percent',10^-2);





