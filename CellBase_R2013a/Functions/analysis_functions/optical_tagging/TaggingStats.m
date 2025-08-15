function [ID_amp,Lr_amp,ID_PC,Lr_PC,Eff,Lat,Jitt,BaseFR,MaxFR,MinFR,lim_act1,lim_act2,lim_inh1,lim_inh2,Wp_act,Wp_inh,R] = TaggingStats(cellid,TrigName,sigma,dt,win)
% Pass the control to the user in case of error
dbstop if error


%             % Add 'PulseOn' event if missing
%             ST = loadcb(cellid,'STIMSPIKES');
%             if isequal(findcellstr(ST.events(:,1),'PulseOn'),0)
%                 prealignSpikes_hp(cellid,'FUNdefineEventsEpochs',...
%                     @defineEventsEpochs_PhotoStim,'filetype','stim',...
%                     'ifsave',1,'ifappend',1)
%             end
                                             
            % Efficiency, latency and jitter for 'chr2'
            [Eff,Lat,Jitt,BaseFR,MaxFR,MinFR,lim_act1,lim_act2,lim_inh1,lim_inh2,Wp_act,Wp_inh]=reliability_latency_jitter(cellid,'event',TrigName,'event_type','stim','window',win,'dt',dt,'sigma',sigma);
                     
if (MaxFR-BaseFR)>(BaseFR-MinFR)
    lim1=lim_act1;
    lim2=lim_act2;
else
    lim1=lim_inh1;
    lim2=lim_inh2;
end
            % Spike shape correlation
            R= spikeshapecorr(cellid,lim1,lim2);     
            %% plot histogram with highlighted stim period
                PlotStimPeriod(cellid,'limwindow',[lim1 lim2]); 
                
            %% save spon. vs. evoked spike waveform
                plotwaveforms2(cellid,'limwindow',[lim1 lim2]); 
            %% save spon.vs evoked mclust projection                  
                plot_mclust_projections2(cellid,'feature_names',{'feature_Peak','feature_Energy'});
        else
            SaltP_act = NaN;
            D_KL_act = NaN;
            SaltP_inh = NaN;
            D_KL_inh = NaN;
            Eff=NaN;
            Lat=NaN;
            Jitt=NaN;
            BaseFR=NaN;
            MaxFR=NaN;
            MinFR=NaN;
            lim_act1=NaN;
            lim_act2=NaN;
            lim_inh1=NaN;
            lim_inh2=NaN;
            Wp_act=NaN;
            Wp_inh=NaN;
            R = NaN;
               
        
    catch ME
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
    end
% end
% keyboard