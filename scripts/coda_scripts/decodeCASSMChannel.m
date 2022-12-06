  function [channel]=decodeCASSMChannel(t, enc,threshold,fig)

% function [channel]=decodeCASSMChannel(t, enc,threshold,fig)
%
% Decodes binary CASSM encoding for source identification. 
%
% Input Arguments
% -------------------
% 1. t         :  array of time values for the encoding data (in ms (!))
% 2. enc       :  the CASSM encoding channel data
% 3. threshold :  clips data below threshold (use .5)
% 4. fig       :  plots results on figure(fig) id fig>0
%  
% Output Variables
% -------------------
% 1. channel   :  the decoded Cytek channel number
%
% History
% -------------------
% 1. First version by Todd Wood
% 2. Edits by JBAF : April 26th, 2018
% 
% Modifications
% -------------------
% 2. Clean-up and documentation during debugging process. Took out a
% "round" statement which was always pushing bit_width to zero (right
% move)? Might need to run this by Todd.
%
%
% Notes
% -------------------
% 

 BYTE_SIZE=5;
 
 enc(find(enc<threshold*max(enc))) = 0; 
 enc(find(enc>0))                  = max(enc); 
 
 % calculating 1st derivative to get bit boundaries
 fDiff                             = abs(diff(enc)); 
 
 % finding the locations with a high derivative 
 b                                 = t(find(fDiff>threshold*max(fDiff)));
 
 % calculating the width and start location from the first peak
 %bit_width                         = round(b(2)-b(1))
 bit_width                         = (b(2)-b(1));
 bit_start                         = b(1);
 
 % this must assume a guard bit - spaced afterwards ....
 data_start = bit_start + (1.5*bit_width);
 bits       = [];
 channel    = 0;
 
 % looping over "bits" in signal 
 for k=1:BYTE_SIZE
     % finding bit location
     bits = [bits,(data_start+k*bit_width)];
     if enc(min(find(t>=bits(k)))) > 0
        channel=bitset(channel,k);
     end
 end
 
 % just some test graphics
 if fig
    figure(fig);
    clf;
    plot(t,enc);
    set(gca,'XLim',[floor(data_start-2*bit_width) data_start+8*bit_width]);
    hold;
    plot(t,[abs(diff(enc)); 0],'g-');
    plot(bits,ones(1,length(bits)) *threshold*max(enc),'r+');
    xlabel('Time');
    ylabel('Signal');
    text(bits(BYTE_SIZE)+bit_width,threshold*max(enc), ...
        num2str(channel),'fontsize',20)
 end

