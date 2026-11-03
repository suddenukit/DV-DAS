%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compression of raw channel RF data based on uniform quantization
% Author(s): Daler Rakhmatov (daler@ece.uvic.ca)
% Date: 19 August 2024, University of Victoria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[RFdata_code,RFdata_pilot] = f_quantCompress_RF(channel_RFdata,partition,codebook)

[Nt,Nx,Na] = size(channel_RFdata);
dataType = class(channel_RFdata);

RFdata_pilot = zeros(Nt,Na,dataType);   % pilot reference (max abs value)
RFdata_code = zeros(Nt,Nx,Na,'uint8');  % quantized data (codebook indices)
% NOTE: 'uint8' --> no more than 8 bits per compressed sample

%-- simulate quantization
for ia = 1:Na

    for it = 1:Nt
        vec = squeeze(channel_RFdata(it,:,ia));

        maxAbs = max(abs(vec));
        [indQ,~] = quantiz(vec,single(maxAbs)*partition,codebook);

        RFdata_pilot(it,ia) = maxAbs;
        RFdata_code(it,:,ia) = indQ;
    end

end

end