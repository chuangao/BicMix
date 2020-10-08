my @nu=(1,0.1,0.01,0.001,0.0001);

for(my $i=0;$i<5;$i++){
print "data $i\n";
`matlab -r "nu=$nu[$i];mine"`;
}
