clear

global head_info target_directory extension filename chan imbuffer
global anabuffer projRange
global contrastRange

% stuff that one might often change
%directory='d:\data\brian\chronic_data\';
target_directory='d:\data\brian\chronic_data\';
%base='c083d';

% stuff that one might rarely change
extension='.cfd';
sweep=0;
chan=1;
contrastRange=[1 256];
projRange=[1 1]; %default settings for first last image in projections
					  %uses first and last image if unchanged

imbuffer(1:512, 1:512)=0;

sweepStr=num2str(sweep);
chanStr=num2str(chan);

Fig1Handle=readIM1;

set(Fig1Handle, 'Pointer', 'default')

%SweepHandle=findobj(Fig1Handle, 'Tag', 'SweepText');
%set(SweepHandle,'string',sweepStr);

%BaseHandle=findobj(Fig1Handle, 'Tag', 'BaseText');
%set(BaseHandle,'String',base);

%DirectoryHandle=findobj(Fig1Handle, 'Tag', 'DirectoryText');
%set(DirectoryHandle,'String',directory);

%ChanHandle=findobj(Fig1Handle, 'Tag', 'ChanText');
%set(ChanHandle,'String','1');

%slider
Slider1Handle=findobj(Fig1Handle, 'Tag', 'FrameSlider');
set(Slider1Handle,'Max', 1);
set(Slider1Handle,'Min', 0);
set(Slider1Handle,'Value', 0);
set(Slider1Handle, 'SliderStep',[0.01,0.1]); 

MaxSlider1Handle=findobj(Fig1Handle, 'Tag', 'FrameMax');
set(MaxSlider1Handle, 'String', '1');

MinSlider1Handle=findobj(Fig1Handle, 'Tag', 'FrameMin');
set(MinSlider1Handle, 'String', '1');

CurrentSlider1Handle=findobj(Fig1Handle, 'Tag', 'FrameCurrent');
set(CurrentSlider1Handle, 'String', '1');


%set image display defaults
iptsetpref('ImshowBorder','tight');
iptsetpref('ImshowTruesize','auto');

%set default images and axes 
DatHandle=findobj('Tag', 'WinDat');
axes(DatHandle);

xx(1:100)=0;
plot(xx);

ImageHandle=findobj('Tag', 'WinImg');
axes(ImageHandle);
imshow(imbuffer);

%FigureHandle=findobj('Tag', 'FigImAn');
truesize([512 512]);

%custom pointers
BoxPointer(1:16, 1:16)=NaN;
BoxPointer(1:16, 1)=2;
BoxPointer(1:16, 16)=2;
BoxPointer(1, 1:16)=2;
BoxPointer(16, 1:16)=2;
BoxPointer(3:14, 3)=1;
BoxPointer(3:14, 14)=1;
BoxPointer(3, 3:14)=1;
BoxPointer(14, 3:14)=1;

LinePointer(1:16, 1:16)=NaN;
LinePointer(1, 1:16)=2;
LinePointer(16, 1:16)=2;
LinePointer(3, 1:16)=1;
LinePointer(14, 1:16)=1;
