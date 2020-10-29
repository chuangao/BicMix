$heredoc1 = <<END;

#\$ -S /bin/bash -cwd
#\$ -o error.out -j y
#\$ -pe high 8
#\$ -l highprio
#\$ -l mem_free=8G 
END

$heredoc2 = <<END;

/opt/apps/MATLAB/R2012b/bin/matlab -nojvm -nodisplay -r "file=\$y;kinit=\$kinit;chuan_script;"

END

# The argument are respectively: Yfile,start_factor #, end_factor #,step_of_factor #, header, which_seed, iteration #, separators, pc #, nu, output directory   

#my @data=("HCE_mergedResiduals.20121009.0.txt","HPC_mergedResiduals.20121009.0.txt","HVC_mergedResiduals.20121009.0.txt");

# change result file, input file, iteration number, pc number, 

#my @nu=(1,1e-4,1e-20);
my $data=("Y","Y_sparse");
my $out_dir='/nfs/igsp/labs/engelhardtlab/cg148/IFA/result/';

`mkdir $out_dir/`;

for(my $i=0;$i<2;$i++){
    my $dataSparse=$data[$i];
    #`mkdir $out_dir/$dataSparse`;
    for(my $j=1;$j<11;$j++){
	my $datai=$dataSparse_$j;
	#`mkdir $out_dir/$datai`;
	my $y="~/SFA/simulation/data/".$datai.".txt";
	for(my $m=0;$m<$n_seed;$m++){
	    `mkdir ${out_dir}/${datai}_sd$m`;
	    my $file="IFA_sd${m}.sh";
	    print "$file\n";
	    open FILE,">$out_dir/$datai/sd$m/$file" or die;
	    print FILE "#!/bin/bash\n\n";
	    print FILE $heredoc1;
	    
	    print FILE "out_dir=$out_dir\n";
	    
	    print FILE $heredoc2;
	    close FILE;
	    `cd $out_dir/$dataSparse/$datai/sd$m/`;
	    #`qsub $out_dir/genesfull_sparse_dense/$file`;
	}
    }
}


=pod
for(my $p=0;$p<1;$p++){
for(my $k=3;$k<4;$k++){
    for(my $i=0;$i<1;$i++){
	`mkdir $data[$i]`;
	for(my $l=0;$l<1;$l++){
	    my $result="result_lam_no_bak_filePC0_PC$pc[$p]_itr100_fac2000_nu_".$nu_v[$k]."_a_".$a."_b_".$b[$l];
	    #my $result="result_filePC0_PC10_itr100_fac500_nu_".$nu_v[$k]."_a_".$a."_b_".$b[$l];
	    `mkdir $data[$i]/$result`;
	    for(my $j=0;$j<$n_seed;$j++){
		my $file="SFA_no_lam_bak".$data[$i]."_pc".$pc[$p]."_nu_".$nu_v[$k]."_a_".$a."_b_".$b[$l]."_seed_".$j.".sh";
		print "$file\n";
		open FILE,">$data[$i]/$file" or die;
		print FILE "#!/bin/bash\n\n";
		print FILE $heredoc1;
		
		print FILE "data_i=$data[$i]\n";
        print FILE "npc=$pc[$p]\n";
		print FILE "nu=$nu_v[$k]\n";
		print FILE "seed_i=$j\n";
		print FILE "result=$result\n";
		print FILE "a=$a\n";
		print FILE "b=$b[$l]";
		print FILE $heredoc2;
		close FILE;
		#`cd $data[$i]`;
		`qsub $data[$i]/$file`;
		#`cd ..`;
	    }
	}
    }
}
}

=cut






































=pod

$heredoc2 = <<END;


# copy the binary and data files to a local directory
cp \$HOME/structure/geno_as_factor/simulations/\$PROGRAM.out \$TMPDIR/\$PROGRAM.out
cp -r \$HOME/structure/geno_as_factor/sim_data/gen_data/\$POP/\$TYPE \$TMPDIR/
cd \$TMPDIR

mkdir \$PROGRAM
mkdir \$PROGRAM/\$POP
mkdir \$PROGRAM/\$POP/\$TYPE
mkdir \$PROGRAM/\$POP/\$TYPE/temp.\$i
mkdir \$PROGRAM/\$POP/\$TYPE/output.\$i

./\$PROGRAM.out 1 3 \$TMPDIR/\$TYPE/y.\$i \$TMPDIR/\$TYPE/geno.\$i \$TMPDIR/\$DIR_DEST/temp.\$i/ \$TMPDIR/\$DIR_DEST/output.\$i/

cp -fr \$TMPDIR/\$PROGRAM/\$POP/\$TYPE/* \$HOME/structure/geno_as_factor/simulations/\$PROGRAM/\$POP/\$TYPE/

END


my @program=("two_step","eHF");
#my @program=("eHF");
my @pop=("pop","no_pop","mix");
#my @pop=("mix");
my @type=("null","no_pleit","pleit");
#my @type=("pleit");

for(my $i=0;$i<2;$i++){
  my $program=$program[$i];
  `mkdir $program`;
  for(my $j=0;$j<3;$j++){
	my $pop=$pop[$j];
	`mkdir $program/$pop`;
	for(my $k=0;$k<3;$k++){
      my $type=$type[$k];
      `mkdir $program/$pop/$type`;
      for(my $m=1;$m<11;$m++){
		my $file="sim_".$program."_".$pop."_".$type."_".$m.".sh";
		print "$file\n";
		open FILE,">$file" or die;
		print FILE "#!/bin/sh\n\n";
		print FILE $heredoc1;
		print FILE "i=$m\nPROGRAM=$program\nPOP=$pop\nTYPE=$type\nDIR_DEST=\$PROGRAM/\$POP/\$TYPE\n";
		print FILE $heredoc2;
		close FILE;
		`nsub $file`;
      }
	}
  }
}


=cut


=pod



$heredoc_vBay = <<END;
# copy the binary and data files to a local directory
cp \$HOME/structure/geno_as_factor/simulations/\$PROGRAM.out \$TMPDIR/\$PROGRAM.out
cp -r \$HOME/structure/geno_as_factor/sim_data/gen_data/\$POP/\$TYPE \$TMPDIR/
cd \$TMPDIR

mkdir \$PROGRAM
mkdir \$PROGRAM/\$POP
mkdir \$PROGRAM/\$POP/\$TYPE
mkdir \$PROGRAM/\$POP/\$TYPE/temp.\$i
mkdir \$PROGRAM/\$POP/\$TYPE/output.\$i

./\$PROGRAM.out 1 3 \$TMPDIR/\$TYPE/y_vBay.\$i \$TMPDIR/\$TYPE/geno.\$i \$TMPDIR/\$DIR_DEST/temp.\$i/ \$TMPDIR/\$DIR_DEST/output.\$i/

cp -fr \$TMPDIR/\$PROGRAM/\$POP/\$TYPE/* \$HOME/structure/geno_as_factor/simulations/\$PROGRAM/\$POP/\$TYPE/

END

my @program=("vBay");
my @pop=("pop","no_pop","mix");
my @type=("null","no_pleit","pleit");

for(my $i=0;$i<1;$i++){
    my $program=$program[$i];
    `mkdir $program`;
    for(my $j=0;$j<3;$j++){
        my $pop=$pop[$j];
        `mkdir $program/$pop`;
        for(my $k=0;$k<3;$k++){
            my $type=$type[$k];
            `mkdir $program/$pop/$type`;
	    for(my $m=1;$m<11;$m++){
		my $file="sim_".$program."_".$pop."_".$type."_".$m.".sh";
		print "$file\n";
		open FILE,">$file" or die;
		print FILE "#!/bin/sh\n\n";
		print FILE $heredoc1;
		print FILE "i=$m\nPROGRAM=$program\nPOP=$pop\nTYPE=$type\nDIR_DEST=\$PROGRAM/\$POP/\$TYPE\n";
		print FILE $heredoc_vBay;
		close FILE;
		`nsub $file`;
	    }
        }
    }
}


=cut
