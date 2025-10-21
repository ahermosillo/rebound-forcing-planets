
#!/usr/bin/perl -w
use List::Util qw(max min sum);

use POSIX;

$node = 128.38;
$inc = 22.82;

$i=0;
$infile = 'list-to-run.txt';
open(DATA,"<$infile");
while(<DATA>)
	{
	$line = $_;
	chomp($line);
	if(!($line=~"#"))
	    {
	    @stuff = split(/\s+/,$line);
	    $p[$i] = int($stuff[0]);
	    $q[$i] = int($stuff[1]);
	    $peri[$i] = int($stuff[3]);
	    #$a0res[$i] = $stuff[4];
	    $i=$i+1;
	    }
	}
close(DATA);
$imax = $i;

resloop: for($i=0;$i<$imax;$i++)
    {
	if(-e("$p[$i]-$q[$i]-sections.txt")){next resloop;}
	if($q[$i] == 1){$ntp = 340;}
	elsif($q[$i] == 2){$ntp = 180;}
	elsif($q[$i] == 3){$ntp = 120;}
	elsif($q[$i] == 4){$ntp = 90;}
	elsif($q[$i] == 5){$ntp = 72;}
	elsif($q[$i] == 6){$ntp = 60;}
	else{$ntp = 180;}
	open(RES,">resonance-parameters.txt");
	print RES "$p[$i].\n";
	print RES "$q[$i].\n";
	#print RES "$a0res[$i]\n";
	print RES "$peri[$i].\n";
	print RES "$inc\n";
	print RES "$node\n";
	print RES "1500.\n";
	print RES "$ntp.\n";
	$nonrestp = $ntp;
	print RES "$nonrestp.\n";
	close(RES);
	
	if(-e("sections.txt"))
	    {
	    $junk = `rm sections.txt`;
	    }

	$junk = `./rebound >> log-$p[$i]-$q[$i].txt`;

	$junk = `perl 04-15-21-calc-lib-tp.pl`;

	$com = "mv sections.txt $p[$i]-$q[$i]-sections.txt";
	system($com);

	$com = "mv lib-parameters.txt $p[$i]-$q[$i]-lib-parameters.txt";
	system($com);
	}

