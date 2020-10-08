$heredoc1 = <<END;

#\$ -S /bin/bash -cwd
#\$ -o error.out -j y
#\$ -pe high 8
#\$ -l highprio
#\$ -l mem_free=8G 
END

$heredoc2 = <<END;

cd /nfs/igsp/labs/engelhardtlab/cg148/IFA/result/\$result
touch test.txt
cp /nfs/igsp/labs/engelhardtlab/cg148/IFA/IFA_chuan.m ./
cp /nfs/igsp/labs/engelhardtlab/cg148/IFA/chuan_script.m ./

/opt/apps/MATLAB/R2012b/bin/matlab -nojvm -nodisplay -r "file=\'\$y\';kinit=\$kinit;chuan_script;quit"

END

my $fac;
my $n_seed=10;
for(my $isim=1;$isim<2;$isim++){
    my @data=("Y","Y_sparse");
    
    `mkdir result`;
    for(my $i_data=0;$i_data<2;$i_data++){
	my $data_i=$data[${i_data}]."_".$isim.".txt";
	if($data_i =~ /sparse/){
	    $fac=20;
	}else{
	    $fac=20;
	}
	for(my $m=0;$m<1;$m++){

	    my $result="file${data[${i_data}]}_${isim}_fac${fac}_sd${m}";
	    `mkdir result/$result/`;

	    my $file="${result}.sh";
	    print "$file\n";
	    open FILE,">/nfs/igsp/labs/engelhardtlab/cg148/IFA/result/$result/$file" or die;
	    print FILE "#!/bin/bash\n\n";
	    print FILE $heredoc1;
	    print FILE "result=$result\n";
	    print FILE "y=/home/igsp1/cg148/SFA/simulation/data/$data_i\n";
	    print FILE "kinit=$fac\n";
	    #print FILE "dir=/nfs/igsp/labs/engelhardtlab/cg148/IFA/result/$result/\n";

	    print FILE $heredoc2;
	    close FILE;
	    #`cd /home/igsp1/cg148/SFA/simulation/result_bfrm/$result/`;
	    #open FILE_TEST, ">test" or die;
	    #print FILE_TEST "test\n";
	    #close FILE_TEST;
	    `qsub /nfs/igsp/labs/engelhardtlab/cg148/IFA/result/$result/$file`;
	}
    }
}
