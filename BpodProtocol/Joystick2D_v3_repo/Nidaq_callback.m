function Nidaq_callback(src,event)
% NIDAQ_CALLBACK
% Append each acquired NI-DAQ data chunk to the in-memory buffer while the
% session is running in the background.
%Callback function for nidaq acquisition. This function is used by a
%listener requiered when nidaq is started in background.
global nidaq

nidaq.ai_data = [nidaq.ai_data;event.Data]; % Append the latest acquired chunk.
        
% ExpectedSize=nidaq.duration*nidaq.sample_rate;
% nidaq.ai_data=NaN(ExpectedSize,size(event.Data));
% nidaq.ai_data = event.Data(1:ExpectedSize,:);
     
end