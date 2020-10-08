$heredoc1 = <<END;

#\$ -S /bin/bash -cwd
#\$ -o error.out -j y
#\$ -pe high 8
#\$ -l highprio
#\$ -l mem_free=8G 
END

$heredoc2 = <<END;

cd ~/results_bicmix/BFRM/
touch test.txt
/nfs/igsp/labs/engelhardtlab/cg148/BFRM/bfrm64 /nfs/igsp/labs/engelhardtlab/cg148/BFRM/result_bfrm/\$result/param_\$result

END

my $fac;
my $n_seed=10;
for(my $isim=1;$isim<11;$isim++){
    my @data=("Y","Y_sparse");
    
    `mkdir result_bfrm`;
    for(my $i_data=1;$i_data<2;$i_data++){
	my $data_i=$data[${i_data}]."_".$isim.".txt";
	if($data_i =~ /sparse/){
	    $fac=20;
	}else{
	    $fac=15;
	}
	for(my $m=0;$m<1;$m++){

	    my $result="file${data[${i_data}]}_${isim}_fac${fac}_sd${m}";
	    `mkdir result_bfrm/$result/`;
	    open FILE_p1,"</nfs/igsp/labs/engelhardtlab/cg148/BFRM/parameters.txt" or die;
	    open FILE_p2,">/nfs/igsp/labs/engelhardtlab/cg148/BFRM/result_bfrm/$result/param_$result" or die;
	    while (my $line = <FILE_p1>) {
		$line =~ s/Y\.txt/\/home\/igsp1\/cg148\/SFA\/simulation\/data\/$data_i/;
		$line =~ s/NLatentFactors \= 20/NLatentFactors \= $fac/;
		$line =~ s/Evol \= 0/Evol \= 1/;
		$line =~ s/EvolVarIn \= 0/EvolVarIn \= 25/;
        $line =~ s/EvolVarInFile \= varin.txt/EvolVarInFile \= \/nfs\/igsp\/labs\/engelhardtlab\/cg148\/BFRM\/varin.txt/;
        $line =~ s/EvolMaximumFactors \= 5/EvolMaximumFactors \= 20/;

		print FILE_p2 $line;
	    }
	    close FILE_p1;
	    close FILE_p2;

	    my $file="${result}.sh";
	    print "$file\n";
	    open FILE,">/nfs/igsp/labs/engelhardtlab/cg148/BFRM/result_bfrm/$result/$file" or die;
	    print FILE "#!/bin/bash\n\n";
	    print FILE $heredoc1;
	    print FILE "result=$result";
	    print FILE $heredoc2;
	    close FILE;
	    `cd /nfs/igsp/labs/engelhardtlab/cg148/BFRM/result_bfrm/$result/`;
	    open FILE_TEST, ">test" or die;
	    print FILE_TEST "test\n";
	    close FILE_TEST;
	    `qsub /nfs/igsp/labs/engelhardtlab/cg148/BFRM/result_bfrm/$result/$file`;
	}
    }
}
