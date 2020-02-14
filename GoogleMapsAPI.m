
%[Longitude, Latitude] = readvars('Dual Litter Bins_Tempe_LatLong.xlsx','Sheet','Sheet2','Range','B6:C243');
%BINDATA=[Longitude, Latitude];

D=zeros(238);
T=zeros(238);
for i=1:238
    for j=1:238
        if i==j
            D(i,j)=0;
            T(i,j)=0;
        elseif i<j
            url1='https://maps.googleapis.com/maps/api/distancematrix/json?units=imperial&origins=%.17g';
            url2='%2c';
            url3='%.17g&destinations=%.17g';
            url4='%2C';
            url5='%.17g';
            url6='&key=AIzaSyDBtFlE8EOUwhsOMEYLHY6mjkH-sQYbHnc';
            % origin(latitude, longitude) destination(lat, long)
            Olat=BINDATA(i,2);
            %rows(i) for origin, columns(j) for destination
            Olong=BINDATA(i,1);
            Dlat=BINDATA(j,2);
            Dlong=BINDATA(j,1);
            url=append(sprintf(url1,Olat),url2,sprintf(url3,Olong,Dlat),url4,sprintf(url5,Dlong),url6);
            options=weboptions('ContentType','auto');
            data=webread(url,options);
            D(i,j)=data.rows.elements.distnace.value; %stores distance in meters
            T(i,j)=data.rows.elements.duration.value; %stores times in seconds
        else
            D(i,j)=D(j,i);
            T(i,j)=T(j,i);
        end
    end
end
%url='https://maps.googleapis.com/maps/api/distancematrix/json?units=imperial&origins=33.418915585739%2c-111.93472391278&destinations=33.4214189471555%2C-111.929466806895&key=';
%options=weboptions('ContentType','auto');
%data=webread(url,options);


