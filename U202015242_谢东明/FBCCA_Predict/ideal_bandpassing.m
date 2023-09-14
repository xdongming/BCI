
function filtered = ideal_bandpassing(input, dim, wl, wh, samplingRate)
 
    if (dim > size(size(input),2))
        error('Exceed maximum dimension');
    end
 
    input_shifted = shiftdim(input,dim-1);
    Dimensions = size(input_shifted);
    
    n = Dimensions(1);
    dn = size(Dimensions,2);
    
    
    Freq = 1:n;
    Freq = (Freq-1)/n*samplingRate;
    mask = Freq > wl & Freq < wh;
    
    Dimensions(1) = 1;
    mask = mask(:);
    mask = repmat(mask, Dimensions);
 
    
    F = fft(input_shifted,[],1);
    
    F(~mask) = 0;
    
    filtered = real(ifft(F,[],1));
    
    filtered = shiftdim(filtered,dn-(dim-1));
end