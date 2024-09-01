function [Circlingarray_days, out]=AngleArrays(Centroidarray_days, Time_array_days, vector)
% vector=true;
centroidsize=size(Centroidarray_days);
centroidsize(length(centroidsize)+1:4)=1;
Circlingarray_days= NaN(centroidsize(1), 1, centroidsize(3), centroidsize(4));
t_diffs=diff(Time_array_days);

for fly = 1:centroidsize(3)
    for trial = 1:centroidsize(4)
%     for batch = 1:size(centroidsize(5))
        if all(isnan(Centroidarray_days(:,1,fly,trial)))
            continue
        end
%         width=max(abs([Centroidarray_days(:,1,fly,trial);Centroidarray_days(:,2,fly,trial)]))*2; 
%         disp('Alt vs Old width')
%         disp(width)
        width=max(Centroidarray_days(:,1,fly,trial)) + min(Centroidarray_days(:,1,fly,trial));                                                         % roiwidth in px
%         disp(width)
        out={};
        out.inx=Centroidarray_days(:,1,fly,trial);
        out.iny=Centroidarray_days(:,2,fly,trial);
        out.theta=atan2(out.iny-width/2,out.inx-width/2);
        out.r=sqrt((out.inx-width/2).^2+(out.iny-width/2).^2)/width;



        out.direction=nan(size(out.inx,1),1);
        out.speed=nan(size(out.inx,1),1);
        out.turning=nan(size(out.inx,1),1);
        %vectorized version because this is slow....

%         whos
        if isduration(t_diffs)
            out.speed(2:end)=sqrt(diff(out.iny)).^2+(diff(out.inx).^2)./seconds(t_diffs(1:end,:,:,trial));
        else
            out.speed(2:end)=sqrt(diff(out.iny).^2+(diff(out.inx).^2))./seconds(seconds(t_diffs(1:end,:,:,trial)));
            disp(min(out.speed))
        end

        out.direction(2:end)=atan2(diff(out.iny),diff(out.inx));
        out.turning(2:end)=atan2(sin(diff(out.direction)),cos(diff(out.direction)));

%         out.turning(out.turning>pi())=out.turning(out.turning>pi())-2*pi();
%         out.turning(out.turning<-pi())=out.turning(out.turning<-pi())+2*pi();
        
%         out.angle=zeros(size(inx,1),1);
%
%         for i=1:size(inx,1)-1                                                % one frame at a time
%             if i > 1a
%                 out.direction(i)=atan2(iny(i)-iny(i-1),inx(i)-inx(i-1));
%                 out.speed(i)=sqrt((iny(i)-iny(i-1))^2+(inx(i)-inx(i-1))^2)/seconds(t_diffs(i));
%                 if i>2
%                     out.turning(i)=out.direction(i)-out.direction(i-1);
%                     if out.turning(i) >pi()
%                         out.turning(i)=out.turning(i)-2*pi();
%                     end
%                     if out.turning(i) < -pi()
%                         out.turning(i)=out.turning(i)+2*pi();
%                     end
%                     
%   
%                 end
%             end
%         end
%         angle2=out.theta
%         angle=out.theta(find(out.speed>.5 & out.r<0.4*width)) - out.direction(find(out.speed>.5 & out.r<0.4*width));
        speedthreshold=4;

        reject_indexes=(out.speed<speedthreshold) | (out.r<0.1);
        out.direction(reject_indexes)=NaN;
        out.turning(reject_indexes)=NaN;
%         out.speed(reject_indexes)=NaN;
%         sreject_indexes=(out.speed<.5) | (out.r<0.1*width);

%         reject_indexes=(out.speed>0);

        angle=nan([length(reject_indexes),1]);
        angle(~reject_indexes)=out.theta(~reject_indexes)-out.direction(~reject_indexes);
%         angle=out.theta.*(out.speed>.1).*(out.r<0.1*width) - out.direction.*(out.speed>.1).*(out.r<0.1*width);
%        angle=out.theta-out.direction;

        if ~isempty(angle)
            angle(angle<0)=angle(angle<0)+2*pi();

%             for j = length(angle)
%                 if angle(j)<0
%                     angle(j)=angle(j)+(2*pi);
%                 end
%             end
            if vector
                angle(~reject_indexes)=sin(angle(~reject_indexes)).*out.speed(~reject_indexes);
            else
                angle=sin(angle);
            end



            %                 ha=histc(angle,bins);
            %                 ha=ha/length(~isnan(angle));
            Circlingarray_days(1:length(angle),1,  fly,trial) = angle;
            %                 histarray(:,fly, trial, batch) = ha;

        elseif isempty(angle)

            Circlingarray_days(:,1, fly,trial) = NaN;
            %                 histarray(:,fly, trial, batch) = NaN;

        end


        out.angle=angle;
%         dbstop
    end
end
end
