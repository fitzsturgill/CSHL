
function out1 = CF4header(filename)
% CF4header reads header of confoc4 or cfcNT1.1
%		CF4header returns a structure containing the header
%		Complete filename required including .cfd
%		
%   Karel Svoboda 1/10/00 Matlab 5.3
%	 svoboda@cshl.org
%
%

fid = fopen(filename, 'r');
if fid < 3; error(['### CfdRead: ', filename, ' NOT found.']); end

%head_info version confoc 4. 
% peculiarities: zpos given in units of steps for the stepper motor; this corresponds to 0.5 nm units
% 

head_info.version=fread(fid, 1, 'short');

head_info.name=fread(fid, 14, 'char');
head_info.name=strcat(deblank(char(head_info.name')),'.cfd');

head_info.user=fread(fid, 16, 'char');
head_info.user=deblank(char(head_info.user'));

head_info.time=fread(fid, 16, 'char');
head_info.time=deblank(char(head_info.time'));

head_info.date=fread(fid, 16, 'char');
head_info.date=deblank(char(head_info.date'));

head_info.grab_start=fread(fid, 1, 'long');
head_info.grab_stop=fread(fid, 1, 'long');

head_info.zpos=fread(fid, 1, 'long');

head_info.hdrsize=fread(fid, 1, 'long');
head_info.cfh_resvd=fread(fid, 44, 'long');
head_info.buflen=fread(fid, 1, 'long');
head_info.pixbuflen=fread(fid, 1, 'long');
head_info.str_poletc=fread(fid, 1, 'long');
head_info.darate=fread(fid, 1, 'long');
head_info.dirate=fread(fid, 1, 'long');
head_info.anaoff=fread(fid, 1, 'long');
head_info.pixels_xy=fread(fid, 2, 'short');
head_info.PrePostpixels_xy=fread(fid, 2, 'short');
head_info.ppp_scan=fread(fid, 2, 'short');
head_info.n_images=fread(fid, 1, 'short');
head_info.chans=fread(fid, 1, 'short');
head_info.acv_res=fread(fid, 4, 'long');
head_info.unsign_xy=fread(fid, 2, 'short');
head_info.scandur=fread(fid, 1, 'float');
head_info.retracedur=fread(fid, 1, 'float');
head_info.StepSize=fread(fid, 1, 'long');
head_info.DwellTime=fread(fid, 1, 'long');
head_info.cflags=fread(fid, 1, 'long');
head_info.num=fread(fid, 1, 'short');
head_info.toss=fread(fid, 1, 'short');
head_info.AcqTimeInt=fread(fid, 1, 'long');
head_info.TacqCounter=fread(fid, 1, 'long');
head_info.acv_rsvd=fread(fid, 41, 'long');
%
head_info.rangex=fread(fid, 1, 'uint16');
head_info.rangex=double(head_info.rangex)*20/65536;
%
head_info.rangey=fread(fid, 1, 'uint16');
head_info.rangey=double(head_info.rangey)*20/65536;
%
head_info.offsetx=fread(fid, 1, 'uint16');
head_info.offsetx=double(head_info.offsetx)*20/65536-10.0;
head_info.offsetx=head_info.offsetx+0.5*head_info.rangex;
%
head_info.offsety=fread(fid, 1, 'uint16');
head_info.offsety=double(head_info.offsety)*20/65536-10.0;
head_info.offsety=head_info.offsety+0.5*head_info.rangey;
%
head_info.park=fread(fid, 1, 'uint16');
head_info.park=double(head_info.park)*20/65536;
%
head_info.orientation=fread(fid, 1, 'uint16');
head_info.orientation=double(head_info.orientation)*20/65536;
%
head_info.discrete=fread(fid, 1, 'long');
head_info.zoom_fac=fread(fid, 1, 'float');
head_info.asv_rsvd=fread(fid, 9, 'long');
head_info.ACM0_Gain=fread(fid, 1,'short');
head_info.ACM0_Offset=fread(fid, 1,'short');
head_info.ACM1_Gain=fread(fid, 1,'short');
head_info.ACM1_Offset=fread(fid, 1,'short');
head_info.asv_rsvd=fread(fid, 48,'long');

out1=head_info;
