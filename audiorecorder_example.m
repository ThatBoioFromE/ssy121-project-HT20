function [audio_recorder] = audiorecorder_example

fs = 2*22050; %sampling frequency
audio_recorder = audiorecorder(fs,24,1);% create the recorder

%ADD USER DATA FOR CALLBACK FUNCTION
audio_recorder.UserData.counter = 1; %initialize a counter in the structure UserData

%attach callback function
time_value = 1; %the timer function will be calles every 1 seconds
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); %use cells to specify function and its arguments

record(audio_recorder); %start recording

end

% CALLBACK FUNCTION
%
%   This function will be called at first once every second (since we set
%   time_value = 1), the counter will increase by one every time the
%   callback function is called. Then, as the counter reaches 10, the trigger
%   time will be set to 0.1 seconds. Now, the function will be called with
%   0.1 s time intervals until the counter reaches 20. When the counter
%   equals 20, the audiorecorder will be stopped.
%
% NOTE: we will only use recObj, i.e., the object that calls this function!
function audioTimerFcn(recObj, event, handles)
    if recObj.UserData.counter < 10
        disp('Slowly!')
        recObj.UserData.counter = recObj.UserData.counter + 1;
        % get current audio data
        rec_data = getaudiodata(recObj);
        % plot the audio data
        figure(1); clf;plot(rec_data);
        %Test the difference between pause, stop, or resume the recObj: "pause(recObj)","stop(recObj)", "resume(recObj)"
%         pause(recObj);
%         resume(recObj);
%         figure(1); clf;plot(rec_data);
%         stop(recObj);
%         resume(recObj);
%         figure(1); clf;plot(rec_data);

    elseif recObj.UserData.counter == 10
        resume(recObj);
        disp('Changing timer value')
        time_value = 0.1;
        set(recObj,'TimerPeriod',time_value); 
        recObj.UserData.counter = recObj.UserData.counter + 1;
        % get current audio data
        rec_data = getaudiodata(recObj);
        % plot the audio data
        figure(1); clf;plot(rec_data);

    elseif recObj.UserData.counter < 20
        disp('Superfast!')
        recObj.UserData.counter = recObj.UserData.counter + 1; 
        % get current audio data
        rec_data = getaudiodata(recObj);
        % plot the audio data
        figure(1); clf;plot(rec_data);
    else
        disp('Example done!') 
        stop(recObj);
    end
end
