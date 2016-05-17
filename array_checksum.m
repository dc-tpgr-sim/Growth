function cs = array_checksum(in)
    
    flatsum = sum(in);
    cs = flatsum/sum((in/flatsum).^2);
    
end