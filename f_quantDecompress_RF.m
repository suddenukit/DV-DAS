%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decompression of raw channel RF data based on uniform quantization
% Author(s): Daler Rakhmatov (daler@ece.uvic.ca)
% Date: 19 August 2024, University of Victoria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[RFdata] = f_quantDecompress_RF(RFdata_code,RFdata_pilot,codebook)

[Nt,Nx,Na] = size(RFdata_code);

%-- simulate dequantization
RFdata = zeros(Nt,Nx,Na);
for ia = 1:Na
    for it = 1:Nt
        indQ = squeeze(RFdata_code(it,:,ia)) + 1;
        RFdata(it,:,ia) = single(RFdata_pilot(it,ia))*codebook(indQ);
    end
end

end