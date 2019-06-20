function ret = handleBackground2(data)

for i=1:512 %handle background
    line=data(:,i);
    line(line==-1)=100001;
    line=sort(line);
    data(:,i)=data(:,i)-median(line);
end

ret = data;

end