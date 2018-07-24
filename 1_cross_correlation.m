% Cross-correlation of two time series
% Completed during MS Geophysics at the University of Western Ontario
% MATLAB code by Benjamin Consolvo and Gerhard Pratt
% Updated in 12/2015

w=[1 2 3 1 2]; e=[1 3 5 7];
nw=length(w); ne=length(e);
s=zeros(1,nw+ne-1); % zeros creates an 1-by-8(nw+ne-1) matrix of zeros.
e=fliplr(e); %this time reverses the time series e, and so cross correlation can be performed
for p=1: nw+ne-1 % p starts at 1, and continues to 8(length of w + length of e - 1).
    for k=max(p-ne+1,1):min(nw,p) % index k can only exist between the numbers 1 and 5.
        % In this line of code, before the calculations 
        % of s(p) occur, the index is limited to only numbers which are actually defined.  
        % The minimum that the code will even look at indices (k) is
        % min(nw,p), which would be k=[1 2 3 4 5 5 5 5].  The max that the
        % code will look at is max(p-ne+1,1), which is k=[1 1 1 1 2 3 4 5].
        %  The reason that the code implements this index of 1-5 only, is
        %  because it makes the program run more quickly.  The program does
        %  not have to run through as many iterations of the for loop.
        s(p)=s(p)+w(k)*e(p-k+1);  % The '+1' is added to (p-k); otherwise the
        % time series e will begin at the zeroth index.
    end
end
disp(w); disp(e); disp(s);
% The output e is the result of the cross-correlation


