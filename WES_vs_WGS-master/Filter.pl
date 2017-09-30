#!/usr/bin/perl -w
# Aziz Belkadi July 2014
# This script exclude from WES and WGS separated vcf files variant with DP<8X GQ<20 and MRR<20% and generates new vcf files
# Input vcf files should be in WES/ and WGS/ repertories.
# Works on UG and HC vcf files

use strict;
no warnings;

mkdir 'Filter';
mkdir'Filter/WES';
mkdir'Filter/WGS';



my @files = glob("WES/*.vcf");
my $nb_vcf = scalar(@files);


for (my $i = 0; $i < $nb_vcf; $i++) {

	my(@name1) = split ('/',$files[$i]);
	my(@name) = split (/\./,$name1[1]);
	open(F, "$files[$i]") || die "Problème à l\'ouverture : $!";
	my($name3)=">Filter/WES/$name[0]_filtered.vcf" ;
	open(FILE,"$name3") || die ("Erreur d'ouverture de Resultat") ;
	
	while (my $ligne1 = <F>){
		
		my (@var) = split ('\t',$ligne1);

		if($ligne1 =~ m/#/ ) {
			
			$" = "\t";
			print FILE "@var";
		
		} else{
		
			my (@var) = split ('\t',$ligne1);
			my $exp = $var[9];
			my @exp1 = split (':',$exp);
			my @exp2 = split (',',$exp1[1]);
			my $size = @exp1;
			if($size == 5){
				if(((($exp1[0]eq"0/1") && (($exp2[0] + $exp2[1]) >= 8) && ($exp2[0]/($exp2[0] + $exp2[1])>=0.2) && ($exp2[1]/($exp2[0] + $exp2[1])>=0.2) && ($exp1[3]>= 20))) || (($exp1[0]eq"1/1") && ((($exp2[0] + $exp2[1]) >=8) && ($exp1[3]>=20) ))){
			
					$" = "\t";
					print FILE "@var";
				
				}
			}elsif($size == 4){
				if(((($exp1[0]eq"0/1") && (($exp2[0] + $exp2[1]) >= 8) && ($exp2[0]/($exp2[0] + $exp2[1])>=0.2) && ($exp2[1]/($exp2[0] + $exp2[1])>=0.2) && ($exp1[2]>= 20))) || (($exp1[0]eq"1/1") && ((($exp2[0] + $exp2[1]) >=8) && ($exp1[2]>=20) ))){
			
					$" = "\t";
					print FILE "@var";
				
				}
			}
		}

	}


close FILE || die "Problème à la fermeture : $!";
close F || die "Problème à la fermeture : $!";


}








my @files = glob("WGS/*.vcf");
my $nb_vcf = scalar(@files);


for (my $i = 0; $i < $nb_vcf; $i++) {

	my(@name1) = split ('/',$files[$i]);
	my(@name) = split (/\./,$name1[1]);
	open(F, "$files[$i]") || die "Problème à l\'ouverture : $!";
	my($name3)=">Filter/WGS/$name[0]_filtered.vcf" ;
	open(FILE,"$name3") || die ("Erreur d'ouverture de Resultat") ;
	
	while (my $ligne1 = <F>){
		
		my (@var) = split ('\t',$ligne1);

		if($ligne1 =~ m/#/ ) {
			
			$" = "\t";
			print FILE "@var";
		
		} else{
		
			my (@var) = split ('\t',$ligne1);
			my $exp = $var[9];
			my @exp1 = split (':',$exp);
			my @exp2 = split (',',$exp1[1]);
			my $size = @exp1;
			if($size == 5){
				if(((($exp1[0]eq"0/1") && (($exp2[0] + $exp2[1]) >= 8) && ($exp2[0]/($exp2[0] + $exp2[1])>=0.2) && ($exp2[1]/($exp2[0] + $exp2[1])>=0.2) && ($exp1[3]>= 20))) || (($exp1[0]eq"1/1") && ((($exp2[0] + $exp2[1]) >=8) && ($exp1[3]>=20) ))){
			
					$" = "\t";
					print FILE "@var";
				}
			}elsif($size == 4){
				if(((($exp1[0]eq"0/1") && (($exp2[0] + $exp2[1]) >= 8) && ($exp2[0]/($exp2[0] + $exp2[1])>=0.2) && ($exp2[1]/($exp2[0] + $exp2[1])>=0.2) && ($exp1[2]>= 20))) || (($exp1[0]eq"1/1") && ((($exp2[0] + $exp2[1]) >=8) && ($exp1[2]>=20) ))){
			
				$" = "\t";
				print FILE "@var";

				}
			}
		}

	}


close FILE || die "Problème à la fermeture : $!";
close F || die "Problème à la fermeture : $!";


}








