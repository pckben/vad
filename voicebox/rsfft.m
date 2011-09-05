function x=srfft(y,n)
%SRFFT    fft of a real symmetric spectrum X=(Y,N)
% The invere fft is the same function but divided by N
% Y contains FIX(1+N/2) complex samples from the spectrum: if argument N
% is specified then Y will be truncated or padded accordingly
% IMPORTANT: If N is odd, it MUST be specified explicitly.




%      Copyright (C) Mike Brookes 1998
%      Version: $Id: rsfft.m,v 1.4 2007/05/04 07:01:39 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isreal(y) error('SRFFT: Input must be real'); end
fl=size(y,1)==1;
if fl y=y(:); end
[m,k]=size(y);
if nargin<2 n=2*m-2;
else
  mm=1+fix(n/2);
  if mm>m y=[y; zeros(mm-m,k)];
  elseif mm<m y(mm+1:m,:)=[];
  end
  m=mm;
end
   x=real(fft([y;y(n-m+1:-1:2,:)]));
   x(m+1:end,:)=[];

if fl x=x.'; end
