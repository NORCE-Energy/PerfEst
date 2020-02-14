function Kout = distance_Loc(Kin,~,~)

% Distance-based localization. Assume observations at every position.  

load inputData.mat;
na=sum(options.actnum);
for I = 1:size(options.staticVar,1)
    
    Kpart = Kin((I-1)*na+1:I*na,:);
    
end

for jdata=1:size(Ke,2)
    fprintf(fidLoc,'%s %d\n','index of column of the Kalman gain ',jdata); % real_obs_index is sparse, fprintf will not accept that
    found = 0;
    
    strType1 = deblank(obsType(jdata,:));
    % check this
    % removing 'W' denoting well
    % removing 'H' denoting historical
    if strType1(1)=='W', strType1=strType1(2:end);  end
    if strType1(end)=='H', strType1=strType1(1:end-1);  end
    strLocation1 = deblank(obsLocation(jdata,:));
    for jj = 1:size(kalmanOptions.alphaDB.location,1)
        strType2 = deblank(kalmanOptions.alphaDB.type(jj,:));
        strLocation2 = deblank(kalmanOptions.alphaDB.location(jj,:));
        % need to make sure that the type string does not have the 'H' in
        % it, i.g. using 'WOPR' instead of 'WOPRH'
        if (~isempty(strfind(strType2,strType1)) && ...
                ~isempty(strfind(strLocation2,strLocation1)) )
            found = 1;
            indexFunc = jj;
            break
        end
    end
    
    if found==1 % if find correct localization function
        fprintf(fidLoc,'%s\t %s\t %s %d\n',obsLocation(jdata,:),obsType(jdata,:), ...
            'index for localization function ',indexFunc);
        % localizing only gridblock variables
        % so Kalman gain for free parameters will not be changed
        if (isfield(kalmanOptions,'iterES') && kalmanOptions.iterES==1) % for iES (Ke = kalman gain)
            for j=1:size(options.staticVar,1)
                idx=0;
                for jj=1:size(kalmanOptions.variable,1) % "kalmanOptions.variable" is specified in defaultKalmanOptions.m
                    if length(deblank(options.staticVar(j,:)))==length(deblank(kalmanOptions.variable(jj,:)))
                        if deblank(options.staticVar(j,:))==deblank(kalmanOptions.variable(jj,:))
                            idx=jj;
                            break
                        end
                    end
                end
                if idx==0
                    error('rfmatlab:preProcess: options should contain staticVar')
                else % "kalmanOptions.index" is specified in defaultKalmanOptions.m
                    Ke(kalmanOptions.index(idx,1):kalmanOptions.index(idx,2),jdata)=...
                        kalmanOptions.alphaDB.func(:,indexFunc).*...
                        Ke(kalmanOptions.index(idx,1):kalmanOptions.index(idx,2),jdata);
                end
            end
        else % for EnKF (Ke = kalman gain)
            for j=1:size(kalmanOptions.variable,1)
                Ke(kalmanOptions.index(j,1):kalmanOptions.index(j,2),jdata)=...
                    kalmanOptions.alphaDB.func(:,indexFunc).*Ke(kalmanOptions.index(j,1):kalmanOptions.index(j,2),jdata);
            end
        end
        if isfield(kalmanOptions,'HEnKFNse')
            error('do not use use distance-based localization in combination with HEnKF')
        end
    end
end
fclose(fidLoc);
