fid=fopen('100_100_r5.txt','r');
R = 1;
G = 1;
B = 0;
RGB=[0   0   1;
    1   0   0;
    0   1   0;
    0   0   1;
    1   1   0;
    1   0   1;
    0   1   1;
    0.67 0   1;
    1 0.5 0;
    0.5 0   0;
    0.5 0.5 0.5 ];
index = 1;
while( ~feof(fid))
    str=fgetl(fid);
  if ~isempty(str);
     S = regexp(str,'\t','split');
     %fprintf('%s\n',char(S(1)));
     %fprintf('%s\n',char(S(2)));
	 x = str2num(char(S(1)));
	 y = str2num(char(S(2)));
     %fprintf('%f\n',x);
     %fprintf('%f\n',y);
    %if(x == 0) && (y == 0)
     %fprintf('change color\n');
     %R = mod((R + 0.5),0.9);
     %G = mod((G + 0.5),0.9);
     %B = mod((B + 0.5),0.9);
     %index = mod((index+1),11) + 1;
     %fprintf('%d\n',index);
    %end
    plot(x,y,'.','color',[R G B]);
    %plot(x,y,'.','color',RGB(index,:));
    hold on;
  end
end
hold on;
fclose(fid);
    
	
