function [Xe,e,A,B]=cross_sampen(x,y,M,r,sflag)
%function [Xe,e,A,B]=cross_sampen(x,y,M,r,sflag);
%
%Input
%
%x,y input data
%M maximum template length
%r matching tolerance
%sflag    flag to standardize signals(default yes/sflag=1) 
%
%Output
%
%Xe sample entropy estimates for m=M-1
%e sample entropy estimates for m=0,1,...,M-1
%A number of matches for m=1,...,M
%B number of matches for m=0,...,M-1 excluding last point
% 
% Error corrected by MAX at row 39  - Not correct standardization
%  -> y=y/sx to  x=x/sx
% 
% Lake, D. E., J. S. Richman, M. P. Griffin, and J. R. Moorman.
% Sample entropy analysis of neonatal heart rate variability. Am J Physiol 2002; 283(3):R789-R797;
% http://ajpregu.physiology.org/content/283/3/R789.abstract
% 
% Richman, J. S. and J. R. Moorman.
% Physiological time series analysis using approximate entropy and sample entropy. Am J Physiol 2000; 278(6):H2039-H2049; http://ajpheart.physiology.org/content/278/6/H2039.abstract

if ~exist('M')|isempty(M),M=5;end
if ~exist('r')|isempty(r),r=.2;end
if ~exist('sflag')|isempty(sflag),sflag=1;end
y=y(:);
x=x(:);
ny=length(y);
nx=length(x);
if sflag>0
   y=y-mean(y);
   sy=sqrt(mean(y.^2));   
   y=y/sy;
   x=x-mean(x);
   sx=sqrt(mean(x.^2));   
   x=x/sx;   
end

lastrun=zeros(nx,1);
run=zeros(nx,1);
A=zeros(M,1);
B=zeros(M,1);
p=zeros(M,1);
e=zeros(M,1);
for i=1:ny
   for j=1:nx
      if abs(x(j)-y(i))<r
         run(j)=lastrun(j)+1;
         M1=min(M,run(j));
         for m=1:M1           
            A(m)=A(m)+1;
            if (i<ny)&(j<nx)
               B(m)=B(m)+1;
            end            
         end
      else
         run(j)=0;
      end      
   end
   for j=1:nx
      lastrun(j)=run(j);
   end
end
N=ny*nx;
B=[N;B(1:(M-1))];
p=A./B;
e=-log(p);
Xe=e(end);