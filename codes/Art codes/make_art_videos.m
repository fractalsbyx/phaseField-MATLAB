close all; clc; clear;
initialFrameIndex = 0;
finalFrameIndex = 1000; 
step = 5;
MAX_QUALITY = 100; % The maximum video quality
VIDEO_FRAME_RATE = 40;

fileFormat = 'MPEG-4';
fileName = "art/artvid_1_40fps.mp4";
fileName = convertStringsToChars(fileName);
curVideo = VideoWriter(fileName, fileFormat);
curVideo.Quality = MAX_QUALITY;
curVideo.FrameRate = VIDEO_FRAME_RATE;

open(curVideo); % Open the current video
for curFrameIndex = initialFrameIndex : step: finalFrameIndex
  curFileName = "art/n_t=" + ...
                num2str(curFrameIndex) + ".png";
  curFileName = convertStringsToChars(curFileName);
  curFrame = imread(curFileName);
  writeVideo(curVideo, curFrame);
end
close(curVideo); % Close the current video
