[LAM,EX]=IFA_chuan(file,kinit);


dlmwrite(strcat('Lambda.txt'), LAM, 'delimiter', '\t','precision', 6);
dlmwrite(strcat('EX.txt'), EX, 'delimiter', '\t','precision', 6);
