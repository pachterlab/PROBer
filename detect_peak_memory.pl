#!/usr/bin/perl

use strict;
use Class::Struct;
use FileHandle;
use Switch;

struct( mem_rec => [ value => '$', mem => '$', line => '$' ]);

# redefine special variable's content
$" = "\n";


my $out_mem = new FileHandle;

if (scalar(@ARGV) != 1) {
    print "Usage: detect_peak_memory.pl configF\n";
    print "  Format of configF:\n";
    print "    line 1: output_name\n";
    print "    line 2: the whole program name e.g. the script name\n";
    print "    next lines: the component name that we want to measure memory use\n";
    exit(-1);
}

my $line;
my ($outFN, $n, @names, @mems) = ();

my $empty = new mem_rec;
$empty->value(0); $empty->mem("0"); $empty->line("");

open(INPUT, $ARGV[0]);
$outFN = <INPUT>; chomp($outFN);
while ($line = <INPUT>) {
    chomp($line);
    push(@names, $line);
    push(@mems, $empty);
}
close(INPUT);

$n = @names; # $names[0] is the global program

$out_mem->open(">$outFN.raw.mem");

my $command = "top -n 1 -u bli -b | grep bli";

my $active = 1;

do {
    $line = `$command`;
    my @arr = split(/\n/, $line);
    for (my $i = 0; $i <= $#arr; $i++) {
	$arr[$i] =~ s/^[ ]+//;
	$arr[$i] =~ s/[ ]+$//;
    }

    my $cur_item;
    my %list = ();

    $active = 0;
    
    for (my $i = 0; $i <= $#arr; $i++) {
	my @tmp = split(/[ \t]+/, $arr[$i]);
	
	my $id = -1;
	for (my $j = 0; $j < $n; $j++) {
	    if (index($tmp[11], $names[$j]) >= 0) { 
		$id = $j;
		last;
	    }
	}
	
	if ($id >= 0) {
	    if ($id == 0) { $active = 1; }
	
	    $cur_item = new mem_rec;
	    $cur_item->value(&convert($tmp[5]));
	    $cur_item->mem($tmp[5]);
	    $cur_item->line($arr[$i]);

	    $list{$id} = $cur_item;
	}
    }

    if ($active) {
	foreach my $prog (keys %list) {
	    print $out_mem $list{$prog}->line."\n";
	    if ($mems[$prog]->value < $list{$prog}->value) {
		$mems[$prog] = $list{$prog};
	    }
	}
	print $out_mem "----------------\n";
	sleep 3;
    }

} while ($active);

$out_mem->close();

open(OUTPUT, ">$outFN.mem");
for (my $i = 0; $i < $n; $i++) {
    print OUTPUT "$names[$i] uses ".$mems[$i]->mem." memory.\n";
}
print OUTPUT "\nEvidence:\n";
for (my $i = 0; $i < $n; $i++) {
    print OUTPUT $mems[$i]->line."\n";
}
close(OUTPUT);

sub convert {
    my $tmp = $_[0];
    my $res;

    switch(chop($tmp)) {
	case 'k' { $res = $tmp * 1024; }
	case 'm' { $res = $tmp * 1024 * 1024; }
	case 'g' { $res = $tmp * 1024 * 1024 * 1024; }
	else { $res = $_[0]; }
    }

    return $res;
}
