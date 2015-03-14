% Script to test HSBounds.
% the double %% separate sections that can be run separately executed in
% MATLAB.  Click in a section to highlight that section and click on the
% "Run Section" arrow in the menu bar
% Test the examples given in 
%   Watt, P. (1987), POLYXSTAL: a FORTRAN program to calculate average elastic properties 
%   of minerals from single-crystal elasticity data, Computers and Geosciences, 13, 441-462.

% Watt 1987 example: zircon.
zircon=[
  422.9000   71.4000  148.7000         0         0         0
   71.4000  422.9000  148.7000         0         0         0
  148.7000  148.7000  490.3000         0         0         0
         0         0         0  112.7000         0         0
         0         0         0         0  112.7000         0
         0         0         0         0         0   48.7000];
% run H-S bounding function:
[hs_zircon,vrh]=HSBounds(zircon)
% should give these values that agree with Watt to 0.01 GPa
% hs_zircon = 
%   228.0300  112.4900
%   226.6700  107.1500 
% vrh =
%   230.4100  119.3000
%   225.1200   98.1400
%   227.7600  108.7200

%%
% Watt1987  example: Stishovite
stish=[
   453   211   203     0     0     0
   211   453   203     0     0     0
   203   203   776     0     0     0
     0     0     0   252     0     0
     0     0     0     0   252     0
     0     0     0     0     0   302];
 % run H-S bounding function:
 [hs_stish,vrh]=HSBounds(stish)
 % should give these values that agree with Watt 1987 to 0.01 GPa
%  hs_stish =
%   317.0900  223.5200
%   314.0600  218.9000
% vrh =
%   324.0000  232.2000
%   308.2900  208.3000
%   316.1400  220.2500
 
%%
% Watt example: muscovite

musc=[
   181.0000   48.8000   25.6000         0  -14.2000         0
   48.8000  178.4000   21.2000         0    1.1000         0
   25.6000   21.2000   58.6000         0    1.0000         0
         0         0         0   16.5000         0   -5.2000
  -14.2000    1.1000    1.0000         0   19.5000         0
         0         0         0   -5.2000         0   72.0000];
% run H-S bounding function:
 [hs_musc,vrh]=HSBounds(musc)       
% Should give these values that agree with Watt  to 0.01 GPa    
% hs_musc =
%    60.6300   37.2400
%    54.2600   32.0500 
% vrh = 
%    67.6800   43.0900
%    48.6700   27.6100
%    58.1800   35.3500

%%

% Plagioclase elastic moduli used to create Table 1. Each column gives data for a
% different composition. The order of elastic constants in each column is:
% C11 C22 C33 C44 C55 C66 C12 C13 C15 C23 C25 C35 C46 C14 C16 C24 C26 C34 C36 C45 C56 
Plag=[
   68.3350   87.0577   96.1748  104.6397  109.3177  120.4077  132.2194
  184.3432  174.9225  189.4431  201.3603  185.4687  191.6068  200.2049
  180.0076  166.0891  171.9478  172.7897  164.1266  163.6588  163.9226
   24.9772   22.9090   23.6214   22.9499   22.1677   23.2520   24.6175
   26.9008   29.0400   33.0783   33.0057   33.1172   32.7587   36.5852
   33.5526   35.0378   35.5199   35.6241   36.8025   35.0393   35.9933
   32.1813   43.9349   46.0683   51.5000   53.0594   56.5598   63.9587
   30.4224   35.3897   38.4058   43.8532   42.1414   49.8835   55.2934
   -2.2533   -0.4010   -0.1520    0.1319    1.2076    3.2471    5.0678
    4.9698   17.9522   15.4201   14.4585   21.9334   26.2717   31.8699
   -7.7979   -2.8696   -5.0887   -4.8006    0.7258    5.3709    3.5499
    7.4850    4.5962    7.2247    6.9137    2.4989    1.6760    0.4826
   -7.1946   -5.1742   -4.8280   -3.8397    1.3602    0.9270   -2.2449
    4.8676    6.1320    5.8709    6.4965    7.6048    9.0466    9.5249
   -0.9297   -0.6003   -0.4071   -0.8010   -7.7187   -3.0004  -10.7509
   -4.3754   -5.8733   -6.9728   -2.4317   -2.9306    2.1181    7.4737
   -6.3760   -6.4679   -6.7722   -9.8523   -6.8497   -9.8744   -7.2259
   -9.1671   -2.9144    2.2077   -0.3603    0.2440    1.7140    6.6292
   -9.3982  -10.6912   -9.8442   -5.7297    0.6965   -8.0794    1.6277
   -2.4157   -1.3141   -1.1262   -1.0003    0.1938    0.7831    2.9544
    0.6107    0.7947    1.3771    2.0575    2.8167    4.5015    5.1855];
% To create the 6x6 matrix
% chose data set:
    data_num=1;  % change this to calculate vales for other compositions in Table 1
% load upper half matrix then use MATLAB commands to create the symmetric lower part   
    cij=zeros(6,6);
    cij(1,1)=Plag(1,data_num);
    cij(2,2)=Plag(2,data_num);
    cij(3,3)=Plag(3,data_num);
    cij(4,4)=Plag(4,data_num);
    cij(5,5)=Plag(5,data_num);
    cij(6,6)=Plag(6,data_num);
    cij(1,2)=Plag(7,data_num);
    cij(1,3)=Plag(8,data_num);
    cij(1,4)=Plag(14,data_num);     
    cij(1,5)=Plag(9,data_num);    
    cij(1,6)=Plag(15,data_num); 
    cij(2,3)=Plag(10,data_num);
    cij(2,4)=Plag(16,data_num);
    cij(2,5)=Plag(11,data_num);
    cij(2,6)=Plag(17,data_num);
    cij(3,4)=Plag(18,data_num);
    cij(3,5)=Plag(12,data_num);
    cij(3,6)=Plag(19,data_num);
    cij(4,5)=Plag(20,data_num);
    cij(4,6)=Plag(13,data_num);
    cij(5,6)=Plag(21,data_num);
    cij=cij+triu(cij,1)';
% run H-S bounding function:
[hs_plag,vrh,ko_go,ustart]=HSBounds(cij)
% for first plagioclase data set should get values that agree with these:
% hs_plag =
%    60.3300   36.7400
%    57.1200   32.8500
% vrh =
%    63.0900   41.4200
%    54.0600   29.8300
%    58.5700   35.6300
% ko_go =
%    39.7561   16.1553
%    70.5533   90.0950
% ustart =
%    16.7678   89.6816
 
