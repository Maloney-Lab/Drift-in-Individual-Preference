function exportdatabyfly()

for duration=[24,2]
    for batch=1:8
%         duration=2;
        % batch=1;
        load(string(batch)+"CentroidArray"+string(duration)+"h.mat", "Centroidarray");
        load(string(batch)+"TimeArray"+string(duration)+"h_u.mat", 'timearray');
        
        Centroidarray_shape=size(Centroidarray);
        if length(Centroidarray_shape)<4
            Centroidarray_shape(4)=1;
        end
        for trial=1:Centroidarray_shape(4)
            timestamps=timearray(:,:,1, trial);
            for fly=1:Centroidarray_shape(3)
                filename="CirclingData_"+duration+"h_B"+batch+"_F"+fly+"_T"+trial;
                disp(filename)
                try
                    [~,b]=AngleArrays(Centroidarray(:,:,fly, trial), timestamps, false);
                    btable=struct2table(b);
                    btable=[btable, table(timestamps)];
        
                    converttable2netcdf(btable, filename, fly, batch, trial, duration);
                end
    
            end
        end
    end
end