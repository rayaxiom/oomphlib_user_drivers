


clear all


NP='NP1';
R='R0';
for i=0:5
  for j=0:5
    % Create the name
    matstr = strcat(['j_' NP R '_'], sprintf('%02d',i),sprintf('%02d',j));
    
    % Load the raw matrix
    str = ['load ./rawdata/' matstr ';'];
    eval(str);
    
    % Convert to matlab sparse matrix.
    str = [matstr '=spconvert([ ' matstr '(:,1)+1, ' matstr '(:,2)+1,' matstr '(:,3)]);'];
    eval(str);
    
  end
end

NP='NP2';
R='R0';
for i=0:5
  for j=0:5
    % Create the name
    matstr = strcat(['j_' NP R '_'], sprintf('%02d',i),sprintf('%02d',j));
    
    % Load the raw matrix
    str = ['load ./rawdata/' matstr ';'];
    eval(str);
    
    % Convert to matlab sparse matrix.
    str = [matstr '=spconvert([ ' matstr '(:,1)+1, ' matstr '(:,2)+1,' matstr '(:,3)]);'];
    eval(str);
    
  end
end

NP='NP2';
R='R1';
for i=0:5
  for j=0:5
    % Create the name
    matstr = strcat(['j_' NP R '_'], sprintf('%02d',i),sprintf('%02d',j));
    
    % Load the raw matrix
    str = ['load ./rawdata/' matstr ';'];
    eval(str);
    
    % Convert to matlab sparse matrix.
    str = [matstr '=spconvert([ ' matstr '(:,1)+1, ' matstr '(:,2)+1,' matstr '(:,3)]);'];
    eval(str);
    
  end
end


Ax = [j_NP1R0_0101,j_NP1R0_0103;
      j_NP1R0_0301,j_NP1R0_0303;];
norm(Ax,inf)