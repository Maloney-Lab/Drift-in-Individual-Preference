path="/Users/ryanmaloney/Documents/GitHub/Continuous Turns/ContinuousData"
continfiles=dir("/Users/ryanmaloney/Documents/GitHub/Continuous Turns/ContinuousData")
% fullfile

addpath('/Users/ryanmaloney/Documents/GitHub/Continuous Turns/ContinuousData')
for i=3:length(continfiles)
%     i=3
    if strcmp(continfiles(i).name(end-3:end), ".mat")

        ff=fullfile(path, continfiles(i).name)
        load(ff)
    
        disp(i)
    % % end
        filename=continfiles(i).name
        ncfilename=strcat(filename(1:end-4),".nc");
    %     nccreate(ncfilename,'angles')
    %     ncwrite(ncfilename, 'angles', Anglearray_days)
    %     S = ncinfo(ncfilename);
    %     S
    %
        [t,nflies]=size(Anglearray_days);
        tarray_unix=posixtime(tarray);
        disp(ncfilename)
    % mode = netcdf.getConstant('NETCDF4');
        ncid=netcdf.create(ncfilename,'NETCDF4');
    % S=ncinfo(ncfilename)
    % S.Format='netcdf4'
    % ncwriteschema(ncfilename,S)
    
        dimid = netcdf.defDim(ncid,'t',t);
        dimid2 = netcdf.defDim(ncid,'fly',nflies);
        varid=netcdf.defVar(ncid,'Angles','NC_DOUBLE',[dimid, dimid2]);
        varid2=netcdf.defVar(ncid,'timestamps','NC_DOUBLE',dimid);
        
        netcdf.endDef(ncid)
        netcdf.putVar(ncid,varid,Anglearray_days)
        netcdf.putVar(ncid,varid2,tarray_unix)
        
        netcdf.close(ncid)
    end
 end

%%
ncid2 = netcdf.open('AngleTstamp1.nc','NC_NOWRITE');
x = netcdf.getVar(ncid2,0);
y = netcdf.getVar(ncid2,1);
whos
%%

ncfilename
