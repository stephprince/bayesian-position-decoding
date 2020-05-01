function  out = gaussSmooth(vector, binrange)
% SP 9.13.18
% taken from function inside of Annabelle's old decoding script 

paddinglength = round(binrange*2.5);
padding = ones(1,paddinglength);

out = smoothvect([padding*vector(1) vector padding*vector(end)],gaussian(binrange,binrange*7));
out = out(paddinglength+1:end-paddinglength);

end