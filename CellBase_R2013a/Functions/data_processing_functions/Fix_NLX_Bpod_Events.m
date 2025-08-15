function [Events_Nttls Events_TimeStamps] = Fix_NLX_Bpod_Events()


% Remove partially-read state-code / timestamp pairs
load Events;
nRealEvents = 0;
RealEvents = [];
RealTimestamps = [];

EventTimestamps=Events_TimeStamps;
EventTTL=Events_Nttls;

Interval = [];
for x = 1:length(EventTTL)-1
    Interval(x) = EventTimestamps(x+1) - EventTimestamps(x);
    if Interval(x) > 50E-6; %50us
        nRealEvents = nRealEvents + 1;
        RealEvents(nRealEvents) = EventTTL(x);
        RealTimestamps(nRealEvents) = EventTimestamps(x);
    end
end
Events_Nttls = RealEvents;
Events_TimeStamps = RealTimestamps;