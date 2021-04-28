% Carey Smith (2020). Bessel Derivative Zeros (https://www.mathworks.com/matlabcentral/fileexchange/28001-bessel-derivative-zeros), MATLAB Central File Exchange. Retrieved July 29, 2020.

function [jprimemk,JofJ]=besselderivzero(MM,KK,tol)

%BessDerivZerosBisect2 Zeros of the first derivative of Bessel function of the first kind.
% function [jprimemk,JofJ]=BessDerivZerosBisect2(MM,KK,tol)
% Calculates the zeros of the first derivatives of Bessel function J'm.
%  MM = a vector of the Bessel function orders (>= 0)
%  KK = a vector of the zero indices (>= 1)
%  tol = tolerance on the value of the derivative, not on the value of the zero
%        (May hang for tol < 1e-12)
%
%Use the bisection algorithm, because it gives the desired roots.
%Estimates for bracketing intervals are taken from the
%asymptotic forms for the roots, in Abramowitz and Stegun.
% Following tradition, only positive zeros are given.  
% So the first zero for J0 is 3.83...
%
% jprimemk = the roots
% JofJ     = Bessel function evaluated there.
%
% Routine written by Larry Forbes, University of Tasmania.
% Modified 2010/06/23  Carey Smith
% Modified 2010/07/09  Carey Smith--Changed the initial guess for large values of m
%                      Carey Smith--Corrected the initial guess for n=17, m=5
% Modified 2011/4/27   Vincent: "It tends to crash for high values of m and k > 1.
%                           (cf. for instance m=44 and k=4)
%                            I fixed it, though, by replacing lines 103 to 110."
if nargin == 2
    tol = 1.e-6;
end
len_m = length(MM);
len_k = length(KK);
jprimemk = zeros(len_m,len_k);
JofJ     = jprimemk;
% Use a table for accuracy and speed. Is is a combination of data from:
%   http://wwwal.kuicr.kyoto-u.ac.jp/www/accelerator/a4/besselroot.htmlx
%   Abramowitz and Stegun (p.411)
%   and this routine
% It is transposed from A&S:
% Rows are values of m (starting at 0); columns are values of k (starting at 1)
%    k=1              k=2              k=3              k=4              k=5              k=6      k=7      k=8      k=9      k=10     k=11     k=12     k=13     k=14     k=15     k=16     k=17     k=18     k=19     k=20
BesselDerivativeZerosT = [...
    3.83170597020751 7.01558666981561 10.1734681350627 13.3236919363142 16.4706300508776 19.6158585105 22.7600843806 25.9036720876 29.0468285349 32.1896799110 35.3323075501 38.4747662348 41.6170942128 44.7593189977 47.9014608872 51.0435351836 54.1855536411 57.3275254379 60.4694578453 63.61136 % m=0
    1.84118378134065 5.33144277352503 8.53631636634628 11.7060049025920 14.8635886339090 18.01553 21.16437 24.31133 27.45705 30.60192 33.74618 36.88999 40.03344 43.17663 46.31960 49.46239 52.60504 55.74757 58.89000 62.03235 % m=1
    3.05423692822714 6.70613319415845 9.96946782308759 13.1703708560161 16.3475223183217 19.51291 22.67158 25.82604 28.97767 32.12733 35.27554 38.42265 41.56893 44.71455 47.85964 51.00430 54.14860 57.29260 60.43635 63.57989 % m=2
    4.20118894121052 8.01523659837595 11.3459243107430 14.5858482861670 17.7887478660664 20.97248 24.14490 27.31006 30.47027 33.62695 36.78102 39.93311 43.08365 46.23297 49.38130 52.52882 55.67567 58.82195 61.96775 65.11315 % m=3
    5.31755312608399 9.28239628524161 12.6819084426388 15.9641070377315 19.1960288000489 22.40103 25.58976 28.76784 31.93854 35.10392 38.26532 41.42367 44.57962 47.73367 50.88616 54.03737 57.18752 60.33677 63.48526 66.63309 % m=4
    6.41561637570024 10.5198608737723 13.9871886301403 17.3128424878846 20.5755145213868 23.80358 27.01031 30.20285 33.38544 36.56078 39.73064 42.89627 46.05857 49.21817 52.37559 55.53120 58.68528 61.83809 64.98980 68.14057 % m=5
    7.50126614468414 11.7349359530427 15.2681814610978 18.6374430096662 21.9317150178022 25.18393 28.40978 31.61788 34.81339 37.99964 41.17885 44.35258 47.52196 50.68782 53.85079 57.01138 60.16995 63.32681 66.48221 69.63635 % m=6
    8.57783648971407 12.9323862370895 16.5293658843669 19.9418533665273 23.2680529264575 26.54503 29.79075 33.01518 36.22438 39.42227 42.61152 45.79400 48.97107 52.14375 55.31282 58.47887 61.64239 64.80374 67.96324 71.12113 % m=7
    9.64742165199721 14.1155189078946 17.7740123669152 21.2290626228531 24.5871974863176 27.88927 31.15533 34.39663 37.62008 40.83018 44.03001 47.22176 50.40702 53.58700 56.76260 59.93454 63.10340 66.26961 69.43356 72.59554 % m=8
    10.7114339706999 15.2867376673329 19.0045935379460 22.5013987267772 25.8912772768391 29.21856 32.50525 35.76379 39.00190 42.22464 45.43548 48.63692 51.83078 55.01844 58.20095 61.37915 64.55368 67.72509 70.89378 74.06014 % m=9
    11.7708766749555 16.4478527484865 20.2230314126817 23.7607158603274 27.1820215271905 30.53451 33.84197 37.11800 40.37107 43.60677 46.82896 50.04043 53.24322 56.43889 59.62863 62.81338 65.99389 69.17075 72.34447 75.51545 % m=10
    ];
  
% n = 17, abs(m)=5:  56.68528
  
[row_bdz, col_bdz] = size(BesselDerivativeZerosT);
for m1=1:len_m
    m = MM(m1);
    for k1=1:len_k
        k = KK(k1);
        if((m+1 <= row_bdz) && (k <= col_bdz))
            % Use the tabble of roots for the asymptotic guess at the kth zero of J'm .
            asymptroot = BesselDerivativeZerosT(m+1,k); % spatial freq.= the nth positive root of d/dr(J_m(r)) = 0
            if(k <= 5)
                dx1 = 1e-9;  % small, because these table values should be very close
            else
                dx1 = 3e-5;  % small, because these table values should be very close
            end
            Aleft = asymptroot-dx1;
            Aright= asymptroot+dx1;
        else
            if(m == 0)
                k = k+1;
            end
            % Give the asymptotic guess at the kth zero of J'm .
            if k == 1
                oneterm = m+0.8086165*m^(1/3)+0.072490*m^(-1/3)-0.05097*m^(-1);
                oneterm = oneterm+0.0094*m^(-5/3);
                asymptroot = oneterm;
            else
                betapr=(m/2+k-3/4)*pi;
                mu=4*m^2;
                asymptroot=betapr-(mu+3)/(8*betapr);
                asymptroot=asymptroot-4*(7*mu^2+82*mu-9)/(3*(8*betapr)^3);
            end
            %bracket the root, and check that the root is bracketed.
            if k == 2
              if(abs(m) <= 1)
                Aleft=asymptroot-pi/1.2;
              else
                Aleft=asymptroot-pi/1.0;
              end
            elseif k == 3
                Aleft=asymptroot-pi/2;
            else
                Aleft=asymptroot-pi/(3+0.02*k);
            end
            Aright=asymptroot+pi/3;
        end
        JprimeL=besselderiv(Aleft,m);
        JprimeR=besselderiv(Aright,m);
        
        %if JprimeL*JprimeR > 0
        %disp(['BessDerivZerosBisect1 error for m=',int2str(m),', k=',int2str(k)])
        %disp(['   Aright=',num2str(Aright),', asymptroot=',num2str(asymptroot),', Aleft=',num2str(Aleft)])
        %warning('BessDerivZerosBisect1:Interval_Error',['Original interval does not bracket root, JprimeL=',num2str(JprimeL),', JprimeR=',num2str(JprimeR)]);
        while JprimeL*JprimeR > 0
          if(JprimeL > 0)
              Aleft = Aleft - (Aright - Aleft)/2;
              JprimeL=besselderiv(Aleft,m);
          else
              Aright = Aright + (Aright - Aleft)/2;
              JprimeR=besselderiv(Aright,m);
          end
        end
        Amiddle=(Aleft+Aright)/2;
        JprimeM=besselderiv(Amiddle,m);
        while abs(JprimeM) > tol
            if JprimeL*JprimeM < 0
                Aright=Amiddle;
                %JprimeR=besselderiv(Aright,m);  % Not used
            else
                Aleft=Amiddle;
                JprimeL=besselderiv(Aleft,m);
            end
            Amiddle=(Aleft+Aright)/2;
            JprimeM=besselderiv(Amiddle,m);
        end
        jprimemk(m1,k1) = Amiddle;
        JofJ(m1,k1) = JprimeM;
    end
end
function derivvalue = besselderiv(x,m)
% Derivative of besselj w.r.t. x
if(0) % Original version
    thederiv=(besselj(m-1,x)-besselj(m+1,x))/2;
    % For unknown reasons, this also seems to return the derivative
else
    dx = 0.0001;
    thederiv = (besselj(m,x+dx) - besselj(m,max(x-dx,0)))/(2*dx);
end
if(imag(thederiv) ~= 0)
    warning('BessDerivZerosBisect1:derivvalue',[' thederiv is complex,not real, m=',int2str(m),', x=',num2str(x)])
    disp([' thederiv = (besselj(',int2str(m),',',num2str(x+dx),')'...
        ' - besselj(',int2str(m),',',num2str(x-dx),'))/(2*',num2str(dx),')'])
    disp(['= (',num2str(besselj(m,x+dx)),' - ',num2str(besselj(m,x-dx)),')/(2*',num2str(dx),')',...
        ' = ',num2str(thederiv),' (complex)'])
    thederiv = real(thederiv);
    disp([' Set thederiv = real(thederiv) = ',num2str(real(thederiv))])
end
derivvalue = thederiv;
