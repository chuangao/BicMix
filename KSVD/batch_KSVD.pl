my @data=("Y_sparse","Y");
my @factor=(10,20,50);
my @atom=(3,5,10);

for(my $i=0;$i<2;$i++){
    for(my $j=0;$j<3;$j++){
        for(my $k=0;$k<3;$k++){
            my $data=$data[i];
            my $fac=$factor[$j];
            my $atom=$atom[$k];
            print "data $data, fac $fac, atom $atom\n";
            `matlab -nojvm -nodisplay -r "param.L=$atom,param.K=$fac,file=$data[$i];mine_sim"`;
        }
    }
}



#my @nu=(1,0.1,0.01,0.001,0.0001);

#for(my $i=0;$i<5;$i++){
#print "data $i\n";
#`matlab -r "nu=$nu[$i];mine"`;
#}
