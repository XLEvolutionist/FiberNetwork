#!usr/bin/perl

##########################################################################################
# A script to trim and parse edge data from a square adjacency matrix to igraph input format
#	Simon Renny-Byfield, Iowa State University, September 2014
##########################################################################################
#usage perl script.pl <input> <output>
use strict;
use warnings;

#open up the data matrix
open ( IN , $ARGV[0] ) || die "Could not open file $ARGV[0]:$!\n";


my $LineNum=0;
my @geneNames;
#the lower limit of adjacency to report the edge 
my $lower = 0.4;
#a hash check to see if an edge has already been given for a specific gene pair
my %check;
#open an output file
open ( OUT , ">$ARGV[1]" ) ;

while ( <IN> ) {
	chomp;
	#the first line is a list of gene names
	if ( m/Gorai/ ) {
		@geneNames = split / /;
	}#if
	else {
		#split the data
		my @adj = split / /;
		for my $i ( 0 .. $#geneNames ) {
			#skip if the adj is 1 or below the limit
			next if ( $adj[$i] == 1 or $adj[$i] < $lower );
			#sort the names, helps with checking of we have already examined this pair
			my @names2sort = ( $geneNames[$LineNum-1] , $geneNames[$i] );
			my @sortedNames = sort ( @names2sort );
			#move on if we have already examined the pair
			next if ( exists $check{$sortedNames[0]}{$sortedNames[1]} );
			#uncomment the next line to print to file and STDOUT
			#print join( "\t" , @sortedNames ) , "\t" ,  $adj[$i] , "\n";
			print OUT join( "\t" , @sortedNames ) , "\t" ,  $adj[$i] , "\n";
			#add the pair to the hash, so it won't be printed again later.
			$check{$sortedNames[0]}{$sortedNames[1]}="yes";
		}
	}#else
	$LineNum++;
}#while

#close the files
close ( IN );
close ( OUT );

exit;
