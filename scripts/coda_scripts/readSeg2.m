function [data,t,notes]=readSeg2(filename,varargin)
%readSeg2: read a seg2 format seismic data file.
% Neill Symons; Sandia National Laboratories; 10/2/03
% Modified by Bruce Engler 2011
% If the SWAMI code is having trouble opening the file, be sure to
%    "addpath" to where the file is located.

%% New 10/3/07, make filename into a cell array.
if ~isa(filename,'cell')
  if size(filename,1)==1
    filename={filename};
  else
    newFilename=cell(1,size(filename,1));
    for i=1:size(filename,1)
      newFilename{i}=filename(i,:);
    end
    filename=newFilename;
    clear newFilename;
  end
end

%% Default values for optional arguments.
verbose=0;

doPlot=1;
params.fontSize=10;

params.printFileStrings=0;
params.printTraceStrings=0;

%% Check varargin for modifiers to the default arguments.
i=1;
while i<=length(varargin)
  currArg=varargin{i};
  i=i+1;
  argType=whos('currArg');
  if ~strcmp(argType.class,'char')
    error('Optional argument %i, type %s must be char',...
      i,argType.class);
  end
  
  switch lower(currArg)      
    case {'verbose' 'v'}
      verbose=verbose+1;
      
    case {'timestamps' 'timestamp' 'tstamp' 'ts'}
      params.timeStamps=varargin{i};i=i+1;
      
    case {'filestrings' 'fstrings' 'strings' 'fstring' 'string'}
      params.printFileStrings=1;
    case {'tracestrings' 'tstrings' 'tstring'}
      params.printTraceStrings=varargin{i};
      i=i+1;
      
    case 'dt'
      params.dt=varargin{i};
      i=i+1;
      
    case {'timemax' 'maxtime' 'tmax' 'maxt'}
      params.taxis=[0 varargin{i}];
      i=i+1;
    case {'timeaxis' 'timelimit' 'taxis'}
      params.taxis=[varargin{i+0} varargin{i+1}];
      i=i+2;
      
    case 'doplot'
      doPlot=1;
    case 'noplot'
      doPlot=0;
      
    case {'tiff' 'save'}
      params.doTiff=1;
    case {'doprint' 'print' 'printer'}
      params.printer=varargin{i};
      i=i+1;
      
    case 'font'
      params.fontSize=varargin{i};
      i=i+1;
    case 'title'
      params.title=varargin{i};
      i=i+1;
    otherwise
      filename=cat(1,filename,currArg);
  end
end

%% Read the files.
for i=1:length(filename)
  if i==1
    [data,notes,params]=readFile(filename{i},params,verbose,doPlot);
    if isfield(params,'timeStamps') && isempty(params.timeStamps{1})
      params.timeStamps{1}=sprintf('%s %s',...
        notes.fileStrings{1}(18:end),notes.fileStrings{2}(18:end));
    end
  else
    % Read the next file.
    [cdata,cnotes,params]=readFile(filename{i},params,verbose,doPlot);
    % Determine the time difference between the first segment and this one.
    if isfield(params,'timeStamps')
      if ~isempty(params.timeStamps{i})
        tdiff=etime(...
          datevec(params.timeStamps{i}),...
          datevec(params.timeStamps{1}));
      else
        tdiff=etime(...
          datevec([cnotes.fileStrings{1}(18:end) cnotes.fileStrings{2}(18:end)]),...
          datevec(params.timeStamps{1}));
      end
    else
      tdiff=etime(...
        datevec([cnotes.fileStrings{1}(18:end) cnotes.fileStrings{2}(18:end)]),...
        datevec([notes.fileStrings{1}(18:end) notes.fileStrings{2}(18:end)]));
    end
    idiff=round(tdiff/params.dt);

    notes.fileStrings=cat(2,notes.fileStrings,cnotes.fileStrings);
    notes.fileNote=cat(2,notes.fileNote,cnotes.fileNote);
    notes.traceStrings=cat(1,notes.traceStrings,cnotes.traceStrings);
    if idiff<0 || idiff>5*size(data,1)
      fprintf('WARNING: Bad date stamp; file %i (%s--%s %s)\n  file 1 (%s--%s %s)\n',...
        i,filename{i},cnotes.fileStrings{1}(18:end),cnotes.fileStrings{2}(18:end),...
        filename{1},notes.fileStrings{1}(18:end),notes.fileStrings{2}(18:end));
     data=cat(1,data,cdata);
    elseif idiff>=size(data,1)
%      data=cat(1,padarray(data,[idiff-size(data,1) 0],'post'),cdata);
     data=cat(1,padZeros(data,[idiff-size(data,1) 0]),cdata);
    else
      newData=zeros([idiff+size(cdata,1) size(data,2)]);
      for ii=1:8:size(data,2)
        if all(data(idiff+1:end,ii)==cdata(1:size(data,1)-idiff,ii))
          if size(data,2)>=ii+7
            newData(:,ii:ii+7)=cat(1,data(1:idiff,ii:ii+7),cdata(:,ii:ii+7));
          else
            newData(:,ii:end)=cat(1,data(1:idiff,ii:end),cdata(:,ii:end));
          end
        else
          jj=1;
          while 1
            if all(data(idiff+1+jj:end,ii)==cdata(1:size(data,1)-idiff-jj,ii))
              if size(data,2)>=ii+7
                newData(:,ii:ii+7)=cat(1,data(1:idiff+jj,ii:ii+7),cdata(1:end-jj,ii:ii+7));
              else
                newData(:,ii:end)=cat(1,data(1:idiff+jj,ii:end),cdata(1:end-jj,ii:end));
              end
              break;
            elseif all(data(idiff+1-jj:end,ii)==cdata(1:size(data,1)-idiff+jj,ii))
              if size(data,2)>=ii+7
                newData(:,ii:ii+7)=cat(1,data(1:idiff-jj,ii:ii+7),cdata(jj:end,ii:ii+7));
              else
                newData(:,ii:end)=cat(1,data(1:idiff-jj,ii:end),cdata(jj:end,ii:end));
              end
              break;
            end
            jj=jj+1;
          end
        end
      end
      data=newData;
    end
  end
end

%% Create the time vector.
dt=params.dt;
t=dt*(0:size(data,1)-1);

%% Do the final stuff.
if doPlot
  plotData=zeros(size(data));
  for i=1:size(data,2)
    if doPlot || nargout>2
      plotData(:,i)=data(:,i)/max(abs(data(:,i)))+i;
    end
  end
  if length(filename)>1
    params.title=strrep(filename{1},'.dat','');
    for i=2:length(filename)
      params.title=sprintf('%s, %s',params.title,strrep(filename{i},'.dat',''));
    end
  end
  doDataPlot(plotData,params)
end



%% LOCAL FUNCTIONS
%

%% function [data,notes,params]=readFile(filename,params,verbose,doPlot)
function [data,notes,params]=readFile(filename,params,verbose,doPlot)
%Open the file.
if isa(filename,'char')
  fid=fopen(filename,'r','ieee-le');
  if fid<0
    error('Unable to open %s',filename);
  end
  %[fid,message]=fopen(filename,'r');
  if ~isfield(params,'title')
    params.title=filename;
  end
else
  t=whos('filename');
  error('Unable to process filename of type %s',t.class);
end
 
%Read the first 2 bytes, should be 3a55 in hex.
A=fread(fid,1,'uint16');
test=sprintf('%0x\n',A);
if ~strncmp(test,'3a55',4)
  error('First 2 bytes should be 3a55 (hex) got %s',test);
end

%Read the first 3 critical parameters from the file.
revNum=fread(fid,1,'uint16');
sizeTpSb=fread(fid,1,'uint16');
nTrace=fread(fid,1,'uint16');

if verbose
  fprintf('File Rev #%i; TracePointerSubblock Size %i; %i Traces\n',revNum,sizeTpSb,nTrace);
end

%String terminator stuff.
fread(fid,6,'uint8');

%bytes 14-31 reserved.
fread(fid,18,'uint8');

%Pointer to trace descripter blocks (1 per trace).
tdbPtr=zeros(1,nTrace);
for i=1:nTrace
  tdbPtr(i)=fread(fid,1,'uint32');
  %if i==1
  %  fprintf('Trace Block Ptr %i: %i\n',i,tdbPtr(i));
  %else
  %  fprintf('Trace Block Ptr %i: %i (interval %i)\n',i,tdbPtr(i),tdbPtr(i)-tdbPtr(i-1));
  %end
end

%Read the strings.
[nBytesRead,fileStrings]=readSeg2Strings(fid,33,sizeTpSb,params.printFileStrings,'File',1);
if nargout>1
  notes.fileStrings=fileStrings;
end

%Read the note.
[nBytesRead,fileNote]=readSeg2Note(fid,nBytesRead,params.printFileStrings);
if nargout>1
  notes.fileNote=fileNote;
end

%Make sure the file is aligned to the correct position.
fseek(fid,tdbPtr(1),'bof');
% %Read any extra required bytes.
% nBytesRead=nBytesRead+4*nTrace;
% if nBytesRead<tdbPtr(1)+1
%   %fprintf('Reading %i extra bytes\n',tdbPtr(i)-nBytesRead+1);
%   fread(fid,tdbPtr(1)-nBytesRead+1,'uint8');
% end

%Read trace descriptor and data blocks.
for i=1:nTrace
  %Check first two bytes
  [curr,nRead]=fread(fid,1,'uint16');
  if nRead~=1
    error('Could not read Trace Descriptor Block ID for trace %i',i);
  end
  if ~strncmp(sprintf('%0x',curr),'4422',4)
    error('Trace Descriptor Block %i: %0x; should be 4422\n',i,curr);
  end
  
  %Read the 4 critical parameters.
  sizeOfBlock=fread(fid,1,'uint16');
  sizeOfData=fread(fid,1,'uint32');
  nSamps=fread(fid,1,'uint32');
  if i==1
    data=zeros(nSamps,nTrace);
  end
  
  
  dataFormat=fread(fid,1,'uint8');
  if verbose
    fprintf('  Trace Block %i: size %i; data size %i; data samps %i; data format %i\n',...
      i,sizeOfBlock,sizeOfData,nSamps,dataFormat);
  end
  switch dataFormat
    case 4
      dataFormatFlag='float32';
    otherwise
      dataFormatFlag='int32';
  end
  
  %Read reserved.
  fread(fid,19,'uint8');
  
  %Read the trace strings.
  [nBytesRead,tStrings]=readSeg2Strings(fid,32,sizeOfBlock,...
    sum(params.printTraceStrings==i),sprintf('Trace #%i',i),1);
  if nargout>1
    notes.traceStrings{i}=tStrings;
  end
  rparams=extractRecordParams(tStrings);
  if ~isfield(params,'dt')
    params.dt=rparams(2);
  elseif params.dt~=rparams(2)
    error('Dt for trace %s.%i (%f); does not match earlier value %f',...
      filename,i,rparams(2),params.dt);
  end
  
  %Read the data.
  [data(:,i),nRead]=fread(fid,nSamps,dataFormatFlag);
  if nRead~=nSamps
    error('Only read %i samples for trace %i; attempted to read %i',nRead,i,nSamps);
  end
  % Nedra added comment 7/27/10 - This converts trace data into mV output.
  data(:,i)=data(:,i)*rparams(1)/rparams(3);
end
fclose(fid);
% END function [data,params]=readFile(filename)

%% START FUNCTION doDataPlot
function doDataPlot(data,params)
if isfield(params,'fftTrace')
  subplot(2,1,1);
end

t=params.dt*(0:size(data,1)-1);
figure;
plot(t,data);
%set(h,'LineWidth',1);
axis tight;
set(gca,'FontSize',params.fontSize,'LineWidth',2);
xlabel('t (s)');
ylabel('Trace #');
if isfield(params,'title')
  title(params.title,'FontSize',params.fontSize,'FontWeight','Bold');
end
if isfield(params,'taxis')
  a=axis;
  axis([params.taxis(1:2) a(3:4)]);
end
if isfield(params,'note')
  a=axis;
  text(0.35*a(1)+0.65*a(2),0.8*a(3)+0.2*a(4),strrep(params.note,'_','\_'),...
    'FontSize',0.75*params.fontSize,'VerticalAlignment','top','BackgroundColor','w');
end

figure(gcf);
drawnow;
if isfield(params,'doTiff') && params.doTiff
  print('-dtiff','-r300',sprintf('%s.tif',strrep(strrep(strrep(params.title,'; ','_'),'.dat',''),'/','_')));
end
if isfield(params,'printer')
  print(sprintf('-P%s',params.printer));
end
%%END FUNCTION doDataPlot

%% START FUNCTION readSeg2Strings
function [startbyte,theStrings]=readSeg2Strings(fid,startbyte,endbyte,doPrint,type,doalign)
if doPrint
  fprintf('%s Strings\n',type);
end
currString=1;
theStrings=cell(0);
while startbyte+1<endbyte
  currSSize=fread(fid,1,'uint16');
  startbyte=startbyte+2;
  if currSSize>2 && currSSize<100
    chars=fread(fid,currSSize-2,'char');
    theStrings{currString}=char(chars)';
    startbyte=startbyte+currSSize-2;
    if doPrint
      fprintf(' %s\n',theStrings{currString});
    end
    currString=currString+1;
  end
end

%Add any required offset.
if doalign && startbyte<endbyte
  %fprintf('Reading %i extra bytes\n',sizeOfBlock-nBytesRead+0);
  fread(fid,endbyte-startbyte,'uint8');
  startbyte=endbyte;
end
%%END FUNCTION readSeg2Strings

%% START FUNCTION readSeg2Note
function [startbyte,theNote]=readSeg2Note(fid,startbyte,doPrint)
currSSize=fread(fid,1,'uint16');
startbyte=startbyte+2;
if currSSize<=2
  theNote='';
else
  chars=fread(fid,currSSize-2,'char');
  startbyte=startbyte+currSSize-2;
  theNote=char(chars);
  if doPrint
    fprintf(' %s\n',theNote);
  end
end
%%END FUNCTION readSeg2Note

%% START FUNCTION extractRecordParams
function [rparams,note]=extractRecordParams(theStrings)

keys={'DESCALING_FACTOR' 'SAMPLE_INTERVAL' 'STACK'};
rparams=zeros(1,length(keys));
for i=1:length(keys)
  key=keys{i};
  index=strmatch(key,theStrings);
  if length(index)~=1
    error('Found %i matches for %s in %i strings (should be exactly 1)',...
      length(index),key,length(theStrings));
  end
  rparams(i)=str2double(theStrings{index}(length(key)+1:end));
  if i>1
    note{i-1}=theStrings{index};
  end
end
%%END FUNCTION extractRecordParams

%% START FUNCTION padZeros
function [d]=padZeros(d,n)
s=size(d);
for i=1:length(s)
  if n(i)>0
    zs=s;
    zs(i)=n(i);
    d=cat(i,d,zeros(zs));
  end
end
%%END FUNCTION padZeros

