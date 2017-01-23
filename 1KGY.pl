#!/usr/bin/perl
$|++;

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

use strict;
use FileHandle;

my $prefix = "d:/perl64/programs/1kg/";
&Main();

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Main(){

   print "\n\n\n\n!*!*!*!*!*!*!*!*!*!*!*!*!*!*!\n\n\n\n";
   print "1KGY Data Extraction and Manipulation Program\n\n";
   print "   (1) Internal Nodes\n";
   print "   (2) Distance Calcs\n";
   print "   (3) Consensii\n";
   print "   (x) Exit\n";
   print "\nEnter Selection >";
   my $selection = <STDIN>; chop $selection; print "\n";
   if ($selection eq "1"){&InternalNodes;}
   if ($selection eq "2"){&ThermodynamicDistanceToNodes;}
   if ($selection eq "3"){&Consensii;}
   if ($selection eq "99"){
      #&ConsensusAllelesToList;
      #&RecodeMega;
      #&CreateFlier;
      #&TestYLocs;
      #&Expand1KG;
      #&Reduce1KGtoCG;
      #&CGPrivates;
      #&TestCGVCF;
      #&TestCGTSV;
      #&Combine;
      #&Nearest;
      #&mutate;
      #&TypeSNPs;
      #&Scozzari;
      #&ISOGGType;
      #&CheckYAlleles;
      #&AnalyzeAncestorFile;
      #&ReduceFASTA;
      #&TestNeanderthalVCF;
      #&DatatoFASTA;
      #&condenseclashes;
      #&Percents;
      #&NandD;
      #&Clashes;
      #&Distance;
      #&CountCategories;
      #&MakeFasta;
      #&FindNearestNeighbor;
      #&PruneData;
      #&RDistanceTable;
   }
   elsif ($selection eq "x"){
      exit;
   }
   my $bell = chr(7); print $bell;
   &Main();
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub InternalNodes(){

   my $suffix = "A0GH.csv";

   print "   Loading Group Designations...\n";
   my $file = "d:/perl64/programs/1KG/megaGroups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $numgroups = @data;
   my @allshapsingroup; my %shapsingroup;
   my %mega;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      $allshapsingroup[$groupnum] = join('', @linedata);
      foreach my $shap (@linedata){
         $mega{$shap} = $groupnum;
      }
   }
   $mega{"F"} = 99;
   $mega{"M"} = 99;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading Scozzari data...\n";
   my $file = "d:/perl64/programs/1KG/Scozzari et al. ancestral data.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data; shift @data;
   my %Scozzari;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $Scozzari{$linedata[0]}{"exists"}    = 1;
      $Scozzari{$linedata[0]}{"data"}      = "$linedata[2]$linedata[3]$linedata[4]$linedata[5]";
      $Scozzari{$linedata[0]}{"ancestral"} = $linedata[4];
   }

   print "   Loading Hallast data...\n";
   my $file = "d:/perl64/programs/1KG/Hallast et al. ancestral data.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data; shift @data;
   my %Hallast;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $Hallast{$linedata[0]}{"exists"}    = 1;
      $Hallast{$linedata[0]}{"data"}      = "$linedata[1]$linedata[2]$linedata[3]";
      $Hallast{$linedata[0]}{"ancestral"} = $linedata[3];
   }

   print "   Loading Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;

   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $mega{$shaps[$i]};
   }

   my @sequences; my @locs; my @majors; my @minors; my $counter = 0; my @missing; my @groupnodes; my @groupalleles; my @allalleles;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      $locs[$counter] = shift @linedata;
      my $major = shift @linedata;
      my $minor = shift @linedata;
      $majors[$counter] = $major;
      $minors[$counter] = $minor;
      $allalleles[$counter][0] = $major;
      $allalleles[$counter][1] = $minor;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         $sequences[$i][$counter] = $linedata[$i];
         my $shap = $shaps[$i];
         $groupalleles[$mega{$shap}][$counter][$linedata[$i]]++;
      }
      $counter++;
   }
   close IN;

#   print "   Calculating Privates...";
#   my %privates; my $private;
#   for (my $loc = 0; $loc <= $counter - 1; $loc++){
#      my $count = 0;
#      for (my $peep = 0; $peep <= $counter - 1; $peep++){
#         if($sequences[$peep][$loc] == 0){ $private = $inds[$peep]; $count++;}
#         if ($count > 1){last;}
#      }
#      if($count == 1){$privates{$private}++;}
#      $count = 0;
#      for (my $peep = 0; $peep <= $counter - 1; $peep++){
#         if($sequences[$peep][$loc] == 1){ $private = $inds[$peep]; $count++;}
#         if ($count > 1){last;}
#      }
#      if($count == 1){$privates{$private}++;}
#   }
#   my $file = "d:/perl64/programs/1KG/YvariantsCP privates.csv";
#   open (OUT, ">$file") || die("Couldn't open >$file");
#   foreach my $ind (sort {$a cmp $b} keys %privates){
#      print OUT "$ind\,$privates{$ind}\n";
#   }
#   close OUT;
#   print "DONE\n";

#   print "   Removing Group Privates...\n";
#   for (my $i = 0; $i <= $counter - 1; $i++){
#      for (my $group = 0; $group <= $numgroups - 1; $group++){
#         if($groupalleles[$group][$i][0] == 1){
#            $groupalleles[$group][$i][0] --;
#            $groupalleles[$group][$i][1] ++;
#         }
#         if($groupalleles[$group][$i][1] == 1){
#            $groupalleles[$group][$i][0] ++;
#            $groupalleles[$group][$i][1] --;
#         }
#      }
#   }

   print "   Calculating Inclusions...";
   my @clashes; my %clashes; my @skips;
   for (my $i = 0; $i <= $counter - 1; $i++){
      for (my $group = 0; $group <= $numgroups - 1; $group++){
         my @alleles;
         # alleles in group
         $alleles[0][0] = $groupalleles[$group][$i][0];
         $alleles[0][1] = $groupalleles[$group][$i][1];
         for (my $outgroup = 0; $outgroup <= $numgroups - 1; $outgroup++){
            if($group != $outgroup){
               # $alleles not in group
               $alleles[1][0] += $groupalleles[$outgroup][$i][0];
               $alleles[1][1] += $groupalleles[$outgroup][$i][1];
            }
         }
         my $incount  = $alleles[0][0] + $alleles[0][1];
         my $outcount = $alleles[1][0] + $alleles[1][1];

         # Strip privates here
#         if($alleles[1][0] == $outcount || $alleles[1][1] == $outcount){ # outgroup fixed
#            if($alleles[0][0] == 1){ # only one seq in in group carries ref
#               $skips[$i] = 1;
#               $alleles[0][0] == 0; $alleles[0][1]++; 
#            }
#            if($alleles[0][1] == 1){ # only one seq in in group carries alt
#               $skips[$i] = 1;
#               $alleles[0][1] == 0; $alleles[0][0]++;
#            }
#         }

         # Strip within-group variability here
 #        if($alleles[1][0] == $outcount || $alleles[1][1] == $outcount){ # outgroup fixed
 #           if($alleles[0][0] > 0 && $alleles[0][1] > 0){ # variability restricted to in group
 #              $skips[$i] = 1;
 #              if($alleles[1][0] == $outcount){
 #                 $alleles[0][0] == $incount; $alleles[0][1] = 0;
 #              }
 #              else{
 #                 $alleles[0][0] == 0; $alleles[0][1] = $incount;
 #              }
 #           }
 #        }

         # is outallele fixed? set to outallele
         if    ($alleles[1][0] == $outcount){$groupnodes[$group][$i] = 0;}# print " outfixed 0\n";}
         elsif ($alleles[1][1] == $outcount){$groupnodes[$group][$i] = 1;}# print " outfixed 1\n";}
         else{
            # is inallele fixed? set to inallele
            if    ($alleles[0][0] == $incount){$groupnodes[$group][$i] = 0;}# print " infixed 0\n";}
            elsif ($alleles[0][1] == $incount){$groupnodes[$group][$i] = 1;}# print " infixed 1\n";}
            else{
               # neither is fixed, apply special cases and tests
               # first, count up groups that are heterozygous at this loc
               my $hetcount; my @vars;
               for(my $j = 0; $j <= $numgroups - 1; $j++){
                  if($groupalleles[$j][$i][0] > 0 && $groupalleles[$j][$i][1] > 0){ # group has variation at this loc
                     $hetcount++;
                  }
                  else{ # what alleles are present within the groups with no variation?
                     $vars[0] += $groupalleles[$j][$i][0]; $vars[1] += $groupalleles[$j][$i][1];  
                  }
               }
               # give up, for now
               $groupnodes[$group][$i] = ".";
               $clashes[$group]++;
               $clashes{$i}++;
            }
         }
      }
   }
   print "\n   Second-pass tests\n";
   my @clashes; my %clashes; my @skips;
   for (my $i = 0; $i <= $counter - 1; $i++){
      my $skipped = 0;
      for (my $group = 0; $group <= $numgroups - 1; $group++){
         if($groupnodes[$group][$i] eq "."){
            $skipped++;
         }
      }
      if($skipped > 0){
         my $ref; my $alt;
         for (my $group = 0; $group <= $numgroups - 1; $group++){
            $ref += $groupalleles[$group][$i][0];
            $alt += $groupalleles[$group][$i][1];
         }
         my $major;
         if($ref >= $alt){$major = 0;}
         else            {$major = 1;}
         for (my $group = 0; $group <= $numgroups - 1; $group++){
            if($groupnodes[$group][$i] eq "."){
#               $groupnodes[$group][$i] = $major;
            }
         }
      }
   }

   print "   Saving stripped FASTA\n";
   my $file = "d:/perl64/programs/1KG/FASTA stripped $suffix";
   open (OUT, ">$file") || die("Couldn't open >$file");
   my @chrs;
   for (my $ind = 0; $ind <= $numinds - 1; $ind++){
      for (my $i = 0; $i <= $counter - 1; $i++){
         if($skips[$i] != 1){
            if($sequences[$ind][$i] eq "."){
               $chrs[$ind] .= ".";
            }
            else{
               $chrs[$ind] .= $allalleles[$i][ $sequences[$ind][$i] ];
            }
         }
      }
   }
   for (my $ind = 0; $ind <= $numinds - 1; $ind++){
      print OUT ">$shaps{$inds[$ind]}-$inds[$ind]\n$chrs[$ind]\n\n";
   }
   close OUT;
   my @skipind;
   for (my $ind = 0; $ind <= $numinds - 1; $ind++){
      for (my $ind2 = $ind + 1; $ind2 <= $numinds - 1; $ind2++){
         if($chrs[$ind] eq $chrs[$ind2]){
            print "      matching chromosomes after removing privates: $ind $ind2\n";
            $skipind[$ind2] = 1;
         }
      }
   }
   for (my $ind = 0; $ind <= $numinds - 1; $ind++){
      if($skipind[$ind] != 1){
         print OUT ">$allshapsingroup[ $shaps{$inds[$ind]} ]-$inds[$ind]\n$chrs[$ind]\n\n";
      }
   }
   close OUT;

   print "   Saving nodes\n";
   my $file = "d:/perl64/programs/1KG/Nodes v2 $suffix";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Loc\,Ref\,Alt\,SAnc\,HAnc";
   for (my $group = 0; $group <= $numgroups - 1; $group++){
      print OUT "\,$allshapsingroup[$group]";
   }
   print OUT "\n";
   for (my $i = 0; $i <= $counter - 1; $i++){
      print OUT "$locs[$i]\,$majors[$i]\,$minors[$i]";
      if($Scozzari{$locs[$i]}{"exists"} == 1){
         print OUT "\,", $Scozzari{$locs[$i]}{"ancestral"};
      }
      else{
         print OUT "\,";
      }
      if($Hallast{$locs[$i]}{"exists"} == 1){
         print OUT "\,", $Hallast{$locs[$i]}{"ancestral"};
      }
      else{
         print OUT "\,";
      }
      for (my $group = 0; $group <= $numgroups - 1; $group++){
         print OUT "\,$groupnodes[$group][$i]";
      }
      print OUT "\,";
      for (my $group = 0; $group <= $numgroups - 1; $group++){
         print OUT "\,$groupalleles[$group][$i][0]\,$groupalleles[$group][$i][1]\,$groupalleles[$group][$i][2]\,";
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub ThermodynamicDistanceToNodes(){

   print "   Loading Group Designations...\n";
   my $file = "d:/perl64/programs/1KG/megaGroups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $counter = 0; my %list; my @groupnames; my %shapsingroup;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      foreach my $group (@linedata){
         $list{$group}{$groupnum}++;
         $groupnames[$counter] .= $group;
         $shapsingroup{$groupnum}{$group}++;
      }
      $counter++;
   }
   my $numgroups = @groupnames;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }
   print "   Loading SNPs Missing in hg19+...\n";
   my %missingSNPs;
   my $file = "d:/perl64/programs/1KG/MissingInHG19.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; close IN;
   foreach my $SNP (@data){
      $missingSNPs{$SNP} = 1;
   }

   print "   Loading Nodes...\n";
   my @nodes;
   my $file = "d:/perl64/programs/1KG/Nodes v2 A0GH.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data; shift @data;
   my $counter = @data;
   for (my $i = 0; $i <= $counter - 1; $i++){
      my @linedata = split(/\,/, $data[$i]);
      my $SNP = shift @linedata; shift @linedata; shift @linedata; shift @linedata; shift @linedata;
      for (my $j = 0; $j <= $numgroups - 1; $j++){
         $nodes[$j][$i] = $linedata[$j];         
      }
   }

   # check each node against each node
   print "   Checking each node against each node...\n";
   my @NNDs; my %NNCs;
   for (my $i = 0; $i <= $counter - 1; $i++){
      for (my $j = 0; $j <= $numgroups - 1; $j++){
         for (my $k = 0; $k <= $numgroups - 1; $k++){
            if($nodes[$j][$i] != $nodes[$k][$i]){
               $NNCs{$j}{$k}{$nodes[$j][$i]}{$nodes[$k][$i]}++;
               $NNDs[$j][$k]++;
            }
         }
      }
   }

   print "   Saving node-to-node distances...\n";
   my $file = "d:/perl64/programs/1KG/Node to Node Distances.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");

   print OUT "Group\,Shaps";
   for (my $i = 0; $i <= $numgroups - 1; $i++){ print OUT "\,$i"; }
   print OUT "\n";
   for (my $i = 0; $i <= $numgroups - 1; $i++){
      print OUT "$i\,$groupnames[$i]";
      for (my $j = 0; $j <= $numgroups - 1; $j++){ print OUT "\,$NNDs[$i][$j]"; }
      print OUT "\n";
   }
   close OUT;

   print "   Calculating...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;
   my %shaps; my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){ $shaps[$i] = $recodes{$inds[$i]} }
      $shaps{$inds[$i]} = $shaps[$i];
   }

   my @NodeD; my %changes; my @counts; my $counter = 0;
   while(<IN>){
      my $line = $_; chomp $line;
      my @linedata = split(/\,/, $line);
      my $SNP = shift @linedata;
      my @alleles; $alleles[0] = shift @linedata; $alleles[1] = shift @linedata;
      # make sure SNP is not missing from HG19
      if($missingSNPs{$SNP} != 1){
         # check each sequence against each node
         for (my $i = 0; $i <= $numinds - 1; $i++){
            if($linedata[$i] < 2){ # makes sure sequence is valid at that position, skip over Ns
               for (my $j = 0; $j <= $numgroups - 1; $j++){
                  if ($nodes[$j][$counter] < 2){ # make sure SNP has not been excluded
                     $counts[$i][$j]++;
                     if($linedata[$i] != $nodes[$j][$counter]){
                        $changes{$i}{$j}{$alleles[$nodes[$j][$counter]]}{$alleles[$linedata[$i]]}++;
                        $NodeD[$i][$j]++;
                     }
                  }
               }
            }
         }
      }
      $counter++;
   }
   close IN;
   
   print "   Saving...\n";
   my $file = "d:/perl64/programs/1KG/Distance to nodes.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for (my $j = 0; $j <= $numgroups - 1; $j++){
      print OUT "\,$groupnames[$j]";
   }
   print OUT "\n";
   print OUT "Ind\,Shap";
   for (my $j = 0; $j <= $numgroups - 1; $j++){
      print OUT "\,$j";
   }
   print OUT "\n";

   for (my $i = 0; $i <= $numinds - 1; $i++){
      print OUT "$inds[$i]\,$shaps[$i]";
      for (my $j = 0; $j <= $numgroups - 1; $j++){
         print OUT "\,$NodeD[$i][$j]";
      }
      print OUT "\n";
   }
   close OUT;

   my $file = "d:/perl64/programs/1KG/Distance to nodes 2.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for (my $j = 0; $j <= $numgroups - 1; $j++){
      print OUT "$j\,$groupnames[$j]\n";
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if($shapsingroup{$j}{$shaps{$inds[$i]}} > 0){
            print OUT "$inds[$i]\,$NodeD[$i][$j]\,$shaps[$i]\n";
         }
      }
      print OUT "\n";
   }
   close OUT;

   exit;

   my $file = "d:/perl64/programs/1KG/Thermodynamic Changes by node final edits 2.0.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "\,";
   my @alleles = ('A', 'C', 'G', 'T');
   my %allowed;
   foreach my $allele (@alleles){
      foreach my $allele2 (@alleles){
         $allowed{$allele}{$allele2}++;
      }
   }
   for (my $i = 0; $i <= $numgroups - 1; $i++){
      foreach my $allele (sort {$a cmp $b} keys %allowed){
         foreach my $allele2 (sort {$a cmp $b} keys %allowed){
            if($allele ne $allele2){
               print OUT "\,$i";
            }
         }
      }
   }
   print OUT "\n";
   print OUT "Group\,Shaps";
   for (my $i = 0; $i <= $numgroups - 1; $i++){
      foreach my $allele (sort {$a cmp $b} keys %allowed){
         foreach my $allele2 (sort {$a cmp $b} keys %allowed){
            if($allele ne $allele2){
               print OUT "\,$allele$allele2";
            }
         }
      }
   }
   print OUT "\n";

   for (my $i = 0; $i <= $numgroups - 1; $i++){
      print OUT "$i\,$groupnames[$i]";
      for (my $j = 0; $j <= $numgroups - 1; $j++){
         foreach my $allele (sort {$a cmp $b} keys %allowed){
            foreach my $allele2 (sort {$a cmp $b} keys %allowed){
               if($allele ne $allele2){
                  print OUT "\,$changes{$i}{$j}{$allele}{$allele2}";
               }
            }
         }
      }
      print OUT "\n";
   }
   close OUT;

   my $file = "d:/perl64/programs/1KG/Thermodynamic Changes to nodes v2 by sequence final edits.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   my @states = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG');
   print OUT "Ratio\,Group\,In/Out\,Results\n";
   foreach my $change1 (@states){
      foreach my $change2 (@states){
         my @alleles1 = split(//, $change1);
         my @alleles2 = split(//, $change2);
         if($alleles1[0] eq $alleles2[1] && $alleles1[1] eq $alleles2[0]){
            for (my $group = 0; $group <= $numgroups - 1; $group++){
               my $out1; my $out2;
               for (my $i = 0; $i <= $numinds - 1; $i++){
                  if($shapsingroup{$group}{$shaps{$inds[$i]}} > 0){
                     if($changes{$i}{$group}{$alleles2[0]}{$alleles2[1]} > 0){
                        $out1 .= $changes{$i}{$group}{$alleles1[0]}{$alleles1[1]} / $changes{$i}{$group}{$alleles2[0]}{$alleles2[1]}; $out1 .= "\,";
                     }
                     else{
                        $out1 .= "\,";
                     }
                  }
                  else{
                     if($changes{$i}{$group}{$alleles2[0]}{$alleles2[1]} > 0){
                        $out2 .= $changes{$i}{$group}{$alleles1[0]}{$alleles1[1]} / $changes{$i}{$group}{$alleles2[0]}{$alleles2[1]}; $out2 .= "\,";
                     }
                     else{
                        $out2 .= "\,";
                     }
                  }
               }
               print OUT "$change1/$change2\,$group\,in\,$out1\n";
               print OUT "$change1/$change2\,$group\,out\,$out2\n";
            }
            print OUT "\n";
         }
      }
   }
   close OUT;

   my $file = "d:/perl64/programs/1KG/Changes from nodes v2 to sequences final edits.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Node\,Groups\,Ind\,Shap\,NodeD\,'+A\,'+C\,'+G\,'+T\,'-A\,'-C\,'-G\,'-T\,AC\,AG\,AT\,CA\,CG\,CT\,GA\,GC\,GT\,TA\,TC\,TG\n";
   for (my $j = 0; $j <= $numgroups - 1; $j++){
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if ($shapsingroup{$j}{$shaps[$i]} == 1){
            print OUT "$j\,$groupnames[$j]\,$inds[$i]\,$shaps[$i]\,$NodeD[$i][$j]";
            my $pA = $changes{$i}{$j}{"C"}{"A"} + $changes{$i}{$j}{"G"}{"A"} + $changes{$i}{$j}{"T"}{"A"};
            my $pC = $changes{$i}{$j}{"A"}{"C"} + $changes{$i}{$j}{"G"}{"C"} + $changes{$i}{$j}{"T"}{"C"};
            my $pG = $changes{$i}{$j}{"A"}{"G"} + $changes{$i}{$j}{"C"}{"G"} + $changes{$i}{$j}{"T"}{"G"};
            my $pT = $changes{$i}{$j}{"A"}{"T"} + $changes{$i}{$j}{"C"}{"T"} + $changes{$i}{$j}{"G"}{"T"};
            my $mA = $changes{$i}{$j}{"A"}{"C"} + $changes{$i}{$j}{"A"}{"G"} + $changes{$i}{$j}{"A"}{"T"};
            my $mC = $changes{$i}{$j}{"C"}{"A"} + $changes{$i}{$j}{"C"}{"G"} + $changes{$i}{$j}{"C"}{"T"};
            my $mG = $changes{$i}{$j}{"G"}{"A"} + $changes{$i}{$j}{"G"}{"C"} + $changes{$i}{$j}{"G"}{"T"};
            my $mT = $changes{$i}{$j}{"T"}{"A"} + $changes{$i}{$j}{"T"}{"C"} + $changes{$i}{$j}{"T"}{"G"};
            print OUT "\,$pA\,$pC\,$pG\,$pT\,$mA\,$mC\,$mG\,$mT";
            print OUT "\,", $changes{$i}{$j}{"A"}{"C"};
            print OUT "\,", $changes{$i}{$j}{"A"}{"G"};
            print OUT "\,", $changes{$i}{$j}{"A"}{"T"};
            print OUT "\,", $changes{$i}{$j}{"C"}{"A"};
            print OUT "\,", $changes{$i}{$j}{"C"}{"G"};
            print OUT "\,", $changes{$i}{$j}{"C"}{"T"};
            print OUT "\,", $changes{$i}{$j}{"G"}{"A"};
            print OUT "\,", $changes{$i}{$j}{"G"}{"C"};
            print OUT "\,", $changes{$i}{$j}{"G"}{"T"};
            print OUT "\,", $changes{$i}{$j}{"T"}{"A"};
            print OUT "\,", $changes{$i}{$j}{"T"}{"C"};
            print OUT "\,", $changes{$i}{$j}{"T"}{"G"};
            print OUT "\n";
         }
      }
      print OUT "\n";
   }
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub ThermodynamicDistanceBetweenNodes(){

   my @shaplist = ("A", "B", "C", "D", "E", "G", "H", "I", "J", "L", "N", "O", "P", "Q", "R", "T", "");
   my @groupnames = ('Con', 'Q', 'R', 'P', 'N', 'O', 'L', 'T', 'I', 'J', 'G', 'H', 'C', 'D', 'E', 'B', 'A', 'RP', 'NO', 'LT', 'IJ', 'DE', 'QRP', 'QRPNOLT', 'IJGH', 'CDE', 'AB');
   my $numnodes = @groupnames;

   print "   Loading Nodes...\n";
   my @nodes;
   my $file = "d:/perl64/programs/1KG/Nodes CP2.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $names = shift @data;
   my @names = split(/\,/, $names); shift @names; shift @names; shift @names; shift @names; 
   my $counter = @data;
   my @alleles; 
   for (my $i = 0; $i <= $counter - 1; $i++){
      my @linedata = split(/\,/, $data[$i]);
      my $loc = shift @linedata;
      $alleles[$i][0] = shift @linedata;;
      $alleles[$i][1] = shift @linedata;
      shift @linedata;
      for (my $j = 0; $j <= $numnodes - 1; $j++){
         if($linedata[$j] eq ""){$linedata[$j] = -1;}
         $nodes[$j][$i] = $linedata[$j];
      }
   }

   my @NodeD;
   my %changes;
   for (my $i = 0; $i <= $counter - 1; $i++){
      for (my $j = 0; $j <= $numnodes - 1; $j++){
         for (my $k = 0; $k <= $numnodes - 1; $k++){
            if($nodes[$j][$i] != $nodes[$k][$i] && $nodes[$j][$i] > -1 && $nodes[$k][$i] > -1){
               $changes{$groupnames[$j]}{$groupnames[$k]}{$alleles[$i][$nodes[$j][$i]]}{$alleles[$i][$nodes[$k][$i]]}++;
               $NodeD[$j][$k]++;
            }
         }
      }
   }

   print "   Saving...\n";
   my $file = "d:/perl64/programs/1KG/Therm Distance between nodes CP2.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "\,\,\,\,";
   for (my $j = 0; $j <= 26; $j++){
      print OUT "\,$groupnames[$j]";
   }
   print OUT "\n";
   print OUT "Node\,MinD\,MinG\,MaxD\,MaxG";
   for (my $j = 0; $j <= 26; $j++){
      print OUT "\,$j";
   }
   print OUT "\n";

   for (my $i = 0; $i <= $numnodes - 1; $i++){
      print OUT "$groupnames[$i]";
      my $minD = 1_000_000; my $minG; my $maxD; my $maxG;
      for (my $j = 0; $j <= $numnodes - 1; $j++){
         if($NodeD[$i][$j] > 0 && $NodeD[$i][$j] < $minD){
            $minD = $NodeD[$i][$j];
            $minG = $j;
         }
         if($NodeD[$i][$j] > 0 && $NodeD[$i][$j] > $maxD){
            $maxD = $NodeD[$i][$j];
            $maxG = $j;
         }         
      }
      print OUT "\,$minD\,$minG($groupnames[$minG])";
      print OUT "\,$maxD\,$maxG($groupnames[$maxG])";
      for (my $j = 0; $j <= $numnodes - 1; $j++){
         print OUT "\,$NodeD[$i][$j]";
      }
      print OUT "\n";
   }
   close OUT;

   my $file = "d:/perl64/programs/1KG/Thermodynamic Changes between nodes CP2.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "\,";
   my @alleles = ('A', 'C', 'G', 'T');
   my %allowed;
   foreach my $allele (@alleles){
      foreach my $allele2 (@alleles){
         $allowed{$allele}{$allele2}++;
      }
   }
   for (my $i = 0; $i <= $numnodes - 1; $i++){
      foreach my $allele (sort {$a cmp $b} keys %allowed){
         foreach my $allele2 (sort {$a cmp $b} keys %allowed){
            if($allele ne $allele2){
               print OUT "\,$i";
            }
         }
      }
   }
   print OUT "\n";
   print OUT "Group";
   for (my $i = 0; $i <= $numnodes - 1; $i++){
      foreach my $allele (sort {$a cmp $b} keys %allowed){
         foreach my $allele2 (sort {$a cmp $b} keys %allowed){
            if($allele ne $allele2){
               print OUT "\,$allele$allele2";
            }
         }
      }
   }
   print OUT "\n";

   foreach my $group (sort {$a cmp $b} keys %changes){
      print OUT "$group";
      for (my $i = 0; $i <= 26; $i++){
         foreach my $allele (sort {$a cmp $b} keys %allowed){
            foreach my $allele2 (sort {$a cmp $b} keys %allowed){
               if($allele ne $allele2){
                  print OUT "\,$changes{$group}{$i}{$allele}{$allele2}";
               }
            }
         }
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub PullData(){

   print "   Pulling data...";

   my %inddata;
   my $file = "d:/perl64/programs/1KG/indpops.txt";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; shift @data;
   foreach my $line (@data){
      my @linedata = split(/\t/, $line);
      $inddata{$linedata[0]}{"pop"}  = $linedata[1];
      $inddata{$linedata[0]}{"spop"} = $linedata[2];
   }
   my $file = "d:/perl64/programs/1KG/Y Chr Haplogroups.txt";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; shift @data;
   foreach my $line (@data){
      my @linedata = split(/\t/, $line);
      $inddata{$linedata[2]}{"macro"} = $linedata[3];
      $inddata{$linedata[2]}{"haplo"} = $linedata[4];
      $inddata{$linedata[2]}{"cont"}  = $linedata[0];
      $inddata{$linedata[2]}{"pop2"}  = $linedata[1];
   }

   my $file = "d:/perl64/programs/1KG/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf";
   open (IN, "<$file") || die("Couldn't open <$file");
   for (my $i = 0; $i <= 124; $i++){
      my $line = <IN>;
   }
   my $headers = <IN>; chomp $headers;
   my @headers = split(/\t/, $headers);
   shift @headers; shift @headers; shift @headers; shift @headers; shift @headers; shift @headers; shift @headers; shift @headers; shift @headers;

   my $file = "d:/perl64/programs/1KG/PopFreqs.csv";
   open (OUT2, ">$file") || die("Couldn't open >$file");
   print OUT2 "Loc\,Maj\,Min\,EAS\,SAS\,EUR\,AFR\,AMR\n";

   my $file = "d:/perl64/programs/1KG/Yvariants.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Loc\,Maj\,Min";
   foreach my $ind (@headers){
      print OUT "\,$ind";
   }
   print OUT "\n";
   print OUT "\,\,";
   foreach my $ind (@headers){
      print OUT "\,", $inddata{$ind}{"pop"};
   }
   print OUT "\n";
   print OUT "\,\,";
   foreach my $ind (@headers){
      print OUT "\,", $inddata{$ind}{"spop"};
   }
   print OUT "\n";
   print OUT "\,\,";
   foreach my $ind (@headers){
      print OUT "\,", $inddata{$ind}{"macro"};
   }
   print OUT "\n";
   print OUT "\,\,";
   foreach my $ind (@headers){
      print OUT "\,", $inddata{$ind}{"haplo"};
   }
   print OUT "\n";

   my %chrsums; my $counter1; my $counter2;
   my $line = <IN>;
   while(<IN>){
      $counter1++;
      my $line = $_;
      chomp $line;
      my @linedata = split(/\t/, $line);
      shift @linedata;
      my $loc = shift (@linedata);
      shift @linedata;
      my $majallele = shift (@linedata);
      my $minallele = shift (@linedata);
      shift @linedata;
      shift @linedata;
      my $info = shift (@linedata);
      my @info = split(/;/, $info);
      my $vt = pop @info; $vt = substr($vt, 3); # SNP, SV, INDEL, MNP
      if($vt eq "TARGET"){ # only happens when loc is an SV
         $vt = pop @info; $vt = substr($vt, 3);
      }
      my $EAS = pop @info; $EAS = substr($EAS, 6);
      my $SAS = pop @info; $SAS = substr($SAS, 6);
      my $EUR = pop @info; $EUR = substr($EUR, 6);
      my $AFR = pop @info; $AFR = substr($AFR, 6);
      my $AMR = pop @info; $AMR = substr($AMR, 6);
      shift @linedata;
      print OUT2 "$loc\,$vt\,$majallele\,$minallele\,$EAS\,$SAS\,$EUR\,$AFR\,$AMR\n";

      print OUT "$loc\,$majallele\,$minallele";
      if($vt eq "SNP"){
         foreach my $ind (@linedata){
            print OUT "\,$ind";
         }
      }
      elsif($vt eq "INDEL"){
         foreach my $ind (@linedata){
            print OUT "\,$ind";
         }
      }
      elsif($vt eq "MNP"){
         foreach my $ind (@linedata){
            print OUT "\,$ind";
         }
      }
      elsif($vt eq "SV"){
         foreach my $ind (@linedata){
            my @data = split(/:/, $ind);
            print OUT "\,$data[0]-$data[1]";
         }
      }
      else{
         print "   New vartype: $vt";
      }
      print OUT "\n";
   }
   close OUT; close OUT2;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub PruneData(){

   print "   Pruning data...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsC.csv";
   open (IN, "<$file") || die("Couldn't open <$file");

   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");

   my $inds  = <IN>; print OUT $inds;
   my $pops  = <IN>; print OUT $pops;
   my $spops = <IN>; print OUT $spops;
   my $shaps = <IN>; print OUT $shaps;
   my $haps  = <IN>; print OUT $haps;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata;
      my $maj = shift @linedata;
      my $min = shift @linedata;
      if(length($maj) == 1 && length($min) == 1){
         print OUT "$line\n";
      }
   }
   close IN;
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub SplitData(){

   print "   Loading data...\n";
   my $file = "d:/perl64/programs/1KG/Yvariants.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds = <IN>; chomp $inds;
   my @inds = split(/\,/, $inds);
   shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;
   my $pops = <IN>; chomp $pops;
   my @pops = split(/\,/, $pops);
   shift @pops; shift @pops; shift @pops;
   my $spops = <IN>; chomp $spops;
   my @spops = split(/\,/, $spops);
   shift @spops; shift @spops; shift @spops;
   my $shaps = <IN>; chomp $shaps;
   my @shaps = split(/\,/, $shaps);
   shift @shaps; shift @shaps; shift @shaps;
   my $haps = <IN>; chomp $haps;
   my @haps = split(/\,/, $haps);
   shift @haps; shift @haps; shift @haps;
   print "   $numinds individuals in database\n";

   my @inddata; my $counter = 0; my $majors; my $minors;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata;
      my $maj = shift @linedata;
      my $min = shift @linedata;
      $majors .= "\,$maj";
      $minors .= "\,$min";
      for (my $i = 0; $i <= $numinds - 1; $i++){
         $inddata[$i] .= "\,$linedata[$i]";
      }
      $counter++;
   }
   close IN;

   print "\n   Saving...\n";
   for (my $i = 0; $i <= $numinds - 1; $i++){
      my $file = "d:/perl64/programs/1KG/Individuals/$inds[$i].csv";
      open (OUT, ">$file") || die("Couldn't open >$file");
      my @seq = split(/\,/, $inddata[$i]);
      print OUT "$inds[$i]";
      foreach my $thing (@seq){
         print OUT "$thing\n";
      }
      close OUT;
   }
   my $file = "d:/perl64/programs/1KG/inddata.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for (my $i = 0; $i <= $numinds - 1; $i++){
      print OUT "$inds[$i]\,$pops[$i]\,$spops[$i]\,$haps[$i]\,$shaps[$i]\n"
   }
   close OUT;
   my $file = "d:/perl64/programs/1KG/MajorAlleles.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   my @seq = split(/\,/, $majors);
   print OUT "Majors\n";
   foreach my $thing (@seq){
      print OUT "$thing\n";
   }
   close OUT;
   my $file = "d:/perl64/programs/1KG/MinorAlleles.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   my @seq = split(/\,/, $minors);
   print OUT "Minors\n";
   foreach my $thing (@seq){
      print OUT "$thing\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Distance(){

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; close IN;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[4];
   }

   print "   Loading Data Clashes...\n";
   my %clashes;
   my $file = "d:/perl64/programs/1KG/Clashes by loc CP2NSH.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; close IN;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[1];
   }

   print "   Loading Data...\n";

   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds = <IN>; chomp $inds;
   my @inds = split(/\,/, $inds);
   shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;
   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $shaps[$i];
   }

   my @inddata; my $counter = 0;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata;
      my $maj = shift @linedata;
      my $min = shift @linedata;
      my $consensus = pop @linedata;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         $inddata[$i][$counter] = $linedata[$i];
      }
      $counter++;
   }
   close IN;

   print "   Calculating:";

   my @differences;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      my @results;
      print " $i";
      for (my $j = $i + 1; $j <= $numinds - 1; $j++){
         my $difference;
         for (my $k = 0; $k <= $counter - 1; $k++){
            if($clashes{$counter} <= 0){
               if($inddata[$i][$k] ne "." && $inddata[$j][$k] ne "." && $inddata[$i][$k] ne $inddata[$j][$k]){
                  $difference++;
               }
            }
         }
         $differences[$i][$j] = $difference;
      }
   }

   print "   Saving:";

   my $file = "d:/perl64/programs/1KG/Y NN Table CP2C.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "\,\,$inds\n\,\,$shaps\n";
   for (my $i = 0; $i <= $numinds - 1; $i++){
      print OUT "\#$shaps[$i]-$inds[$i]\,";
      for (my $j = 0; $j <= $numinds - 1; $j++){
         print OUT "\,$differences[$i][$j]";
      }
      print OUT "\n";
   }
   close OUT;

   my $file = "d:/perl64/programs/1KG/NNTable2C.meg";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "#MEGA\n";
   print OUT "!Title 1KGY Node NN data;\n";
   print OUT "!Description 1KGY Node NN data;\n";
   for (my $i = 0; $i <= $numinds - 1; $i++){
      print OUT "\#$shaps[$i]-$inds[$i]\n";
   }
   for (my $i = 0; $i <= $numinds - 1; $i++){
      my $j = "\t" x $i;
      print OUT $j;
      for (my $j = $i + 1; $j <= $numinds - 1; $j++){
         print OUT "\t$differences[$i][$j]";
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub GCbias(){

   # Open up each child and each parent
   # scan down column 4 for singletons
   # column 3 has allele: $results{$allele}
   # column 2 has allele + other allele: $results{$allele}{$altallele}
   
   my %results;
   print "   Testing GC bias: Child";
   for(my $child = 1; $child <= 11; $child++){
      print " $child";
       for(my $p = 0; $p <= 1; $p++){
         print ".$p";
         my $file = "d:/perl64/programs/1463/testvariants HQ4ABL $child-$p.csv";
         open (IN, "<$file") || die("Couldn't open <$file");
         my @data = <IN>; close IN;
         my @inddata;
         foreach my $line (@data){
            chomp $line;
            my @linedata = split(/\,/, $line);
            push (@inddata, $linedata[4]);
         }
         my $count = @inddata; my $laststrand = "A"; my $counter = 1;
         for(my $line = 0; $line <= $count - 1; $line++){
            my $strand = $inddata[$line];
            if($laststrand eq $strand){$counter++;}
            else{
               if($counter == 1){
                  my @linedata = split(/\,/, $data[$line]);
                  my $allele = $linedata[3];
                  my $alt1 = substr($linedata[2], 0, 1);
                  my $alt2 = substr($linedata[2], 1, 1);
                  if($allele eq $alt1){
                     $results{$child}{$p}{$alt2}{$allele}++;
                     $results{12}{0}{$alt2}{$allele}++;
                  }
                 else{
                     $results{$child}{$p}{$alt1}{$allele}++;
                     $results{12}{0}{$alt1}{$allele}++;
                  }
               }
               $counter = 1;
               $laststrand = $strand;
            }
         }
      }
   }
   my $file = "d:/perl64/programs/1463/GC Bias.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for(my $child = 1; $child <= 12; $child++){
      for(my $p = 0; $p <= 1; $p++){
         my $parent;
         if($p == 0){
            $parent = "F";
         }
         else{
            $parent = "M";
         }
         print OUT "Child $child-$parent\,A\,C\,G\,T\n";
         print OUT "A\,", $results{$child}{$p}{"A"}{"A"};
         print OUT "\,",  $results{$child}{$p}{"A"}{"C"};
         print OUT "\,",  $results{$child}{$p}{"A"}{"G"};
         print OUT "\,",  $results{$child}{$p}{"A"}{"T"};
         print OUT "\n";
         print OUT "C\,", $results{$child}{$p}{"C"}{"A"};
         print OUT "\,",  $results{$child}{$p}{"C"}{"C"};
         print OUT "\,",  $results{$child}{$p}{"C"}{"G"};
         print OUT "\,",  $results{$child}{$p}{"C"}{"T"};
         print OUT "\n";
         print OUT "G\,", $results{$child}{$p}{"G"}{"A"};
         print OUT "\,",  $results{$child}{$p}{"G"}{"C"};
         print OUT "\,",  $results{$child}{$p}{"G"}{"G"};
         print OUT "\,",  $results{$child}{$p}{"G"}{"T"};
         print OUT "\n";
         print OUT "T\,", $results{$child}{$p}{"T"}{"A"};
         print OUT "\,",  $results{$child}{$p}{"T"}{"C"};
         print OUT "\,",  $results{$child}{$p}{"T"}{"G"};
         print OUT "\,",  $results{$child}{$p}{"T"}{"T"};
         print OUT "\n";
      }
      print OUT "\n";
   }
   close OUT;   
}


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Skittlize(){

   # traditional: green = A   blue = C   black = G   red = T
   # Skittle:     black = A   red  = C   green = G   blue = T

   my @centromeres;
   my $file = "d:/perl64/programs/1463/centromerelocations.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; close IN;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $centromeres[$linedata[0]] = $linedata[1]; 
   }

   my $centromerecode = "G" x 1000;

   my $file = "d:/perl64/programs/1463/testvariants HQ4AB.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $line = <IN>; chomp $line;
   my @chroms;
   my $counter;
   print "   Parsing data: 1";
   my $chr = 1; my $lastchr = 1; my $centromere = $centromeres[$chr];
   while(<IN>){
      $counter++;
      my $linedata = $_;
      chomp $linedata;
      my @linedata = split(/\,/, $linedata);
      my @fatdat = split(//, $linedata[6]);
      my @motdat = split(//, $linedata[7]);
      for (my $i = 0; $i <= 10; $i++){
         $chr = $linedata[1];
         if($chr != $lastchr){
            print " $chr";
            for (my $i = 0; $i <= 10; $i++){
               my $child = $i + 1;
               my $file = "d:/perl64/programs/1463/haplotypes/$child $lastchr F.fa";
               open (OUT, ">$file") || die("Couldn't open >$file");
               print OUT ">child $i chr $lastchr F\n";
               print OUT $chroms[$i][0], "\n";
               close OUT;
               $file = "d:/perl64/programs/1463/haplotypes/$child $lastchr M.fa";
               open (OUT, ">$file") || die("Couldn't open >$file");
               print OUT ">child $i chr $lastchr M\n";
               print OUT $chroms[$i][1], "\n";
               close OUT;
               $chroms[$i][0] = "";
               $chroms[$i][1] = "";
               if($fatdat[$i] eq "A"){$chroms[$i][0] = "C";}
               if($fatdat[$i] eq "B"){$chroms[$i][0] = "T";}
               if($motdat[$i] eq "A"){$chroms[$i][1] = "C";}
               if($motdat[$i] eq "B"){$chroms[$i][1] = "T";}
            }
            $lastchr = $chr;
            $centromere = $centromeres[$chr];
         }
         else{
            if($linedata[2] > $centromere){
               for (my $i = 0; $i <= 10; $i++){
                  $chroms[$i][0] .= $centromerecode;
                  $chroms[$i][1] .= $centromerecode;
                  $centromere = 1_000_000_000_000;
               }
            }
            if($fatdat[$i] eq "A"){$chroms[$i][0] .= "C";}
            if($fatdat[$i] eq "B"){$chroms[$i][0] .= "T";}
            if($motdat[$i] eq "A"){$chroms[$i][1] .= "C";}
            if($motdat[$i] eq "B"){$chroms[$i][1] .= "T";}
         }
      }
   }
   print "   NumSNPs = $counter\n";
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub SuperSkittle(){

   # traditional: green = A   blue = C   black = G   red = T
   # Skittle:     black = A   red  = C   green = G   blue = T

   print "   Creating SuperSkittle maps\n";

   my %allowed; $allowed{"A"} = 1; $allowed{"B"} = 1; $allowed{"M"} = 0; $allowed{"X"} = 0; $allowed{"E"} = 0;

   my @centromeres;
   my $file = "d:/perl64/programs/hapmap/3generations/centromerelocations.csv"; 
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; close IN;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $centromeres[$linedata[0]] = $linedata[1]; 
   }

   my $centromerecode = "G" x 400;
   my $band           = "A" x 400;

   my $file = "d:/perl64/programs/1463/testvariants HQ4AB.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $line = <IN>;
   print "   Parsing data:";
   my $chr = 0; my $lastchr = 0; my @chroms; my $centromereF = $centromeres[$chr]; my $centromereM = $centromeres[$chr];
   my $counterF = 0; my $counterM = 0;
   
   while(<IN>){
      my $linedata = $_;
      chomp $linedata;
      my @linedata = split(/\,/, $linedata);
      my @fatdat   = split(//, $linedata[6]);
      my @motdat   = split(//, $linedata[7]);
      $chr = $linedata[1];

      # If new chromosome
      if($chr != $lastchr){
         # First, save data
         for (my $i = 0; $i <= 10; $i++){
            print OUTF "$chroms[$i][0]";
            print OUTM "$chroms[$i][1]";
            $chroms[$i][0] = "";
            $chroms[$i][1] = "";
         }
         close OUTF; close OUTM;
         # Next, start new string with this line of data;
         my $countF  = $linedata[6] =~ tr/A//;
            $countF += $linedata[6] =~ tr/B//;
         my $countM += $linedata[7] =~ tr/A//;
            $countM += $linedata[7] =~ tr/B//;
         for (my $i = 0; $i <= 10; $i++){
            if($countF >= 1){
               if($allowed{$fatdat[$i]} == 1){
                  if($fatdat[$i] eq "A"){$chroms[$i][0] = "C";}
                  if($fatdat[$i] eq "B"){$chroms[$i][0] = "T";}
               }
               else{$chroms[$i][0] = "A";}
            }
            else{$chroms[$i][0] = "";}
            if($countM >= 1){
               if($allowed{$motdat[$i]} == 1){
                  if($motdat[$i] eq "A"){$chroms[$i][1] = "C";}
                  if($motdat[$i] eq "B"){$chroms[$i][1] = "T";}
               }
               else{$chroms[$i][1] = "A";}
            }
            else{$chroms[$i][0] = "";}
         }
         $centromereF = $centromeres[$chr];
         $centromereM = $centromeres[$chr];
         # Next, open files for next chromosome
         my $file = "d:/perl64/programs/1463/haplotypes/Family chr$chr F.fa";
         open (OUTF, ">$file") || die("Couldn't open >$file");
         print OUTF ">chr $lastchr F\n";
         my $file = "d:/perl64/programs/1463/haplotypes/Family chr$chr M.fa";
         open (OUTM, ">$file") || die("Couldn't open >$file");
         print OUTM ">chr $lastchr M\n";
         $lastchr = $chr;
         print " $chr";
         $counterF = 0; $counterM = 0;
      }
      # Still same chromosome
      else{
         my $countF  = $linedata[6] =~ tr/A//;
            $countF += $linedata[6] =~ tr/B//;
         my $countM += $linedata[7] =~ tr/A//;
            $countM += $linedata[7] =~ tr/B//;
         if($countF >= 1){
           if($linedata[2] > $centromereF){
               my $num = 400 - $counterF;
               my $centromerecode = "G" x $num;
               for (my $i = 0; $i <= 10; $i++){
                  print OUTF "$chroms[$i][0]$centromerecode";
                  $chroms[$i][0] = "";
               }
               print OUTF "$band$band$band";
               $centromerecode = "G" x 400;
               for (my $i = 0; $i <= 10; $i++){
                  print OUTF $centromerecode;
               }
               print OUTF "$band$band$band";
               $centromereF = 1_000_000_000_000;
               $counterF = 0;
            }
            $counterF++;
            for (my $i = 0; $i <= 10; $i++){
               if($allowed{$fatdat[$i]} == 1){
                  if($fatdat[$i] eq "A"){$chroms[$i][0] .= "C";}
                  if($fatdat[$i] eq "B"){$chroms[$i][0] .= "T";}
               }
               else{$chroms[$i][0] .= "A";}
            }
            if($counterF == 400){
               for (my $i = 0; $i <= 10; $i++){
                  if($i != 3){
                     print OUTF "$chroms[$i][0]";
                  }
                  $chroms[$i][0] = "";
               }
               print OUTF "$band$band$band";
               $counterF = 0;
            }
         }
         if($countM >= 1){
           if($linedata[2] > $centromereM){
               my $num = 400 - $counterM;
               my $centromerecode = "G" x $num;
               for (my $i = 0; $i <= 10; $i++){
                  print OUTM "$chroms[$i][1]$centromerecode";
                  $chroms[$i][1] = "";
               }
               print OUTM "$band$band$band";
               $centromerecode = "G" x 400;
               for (my $i = 0; $i <= 10; $i++){
                  print OUTM $centromerecode;
               }
               print OUTM "$band$band$band";
               $centromereM = 1_000_000_000_000;
                $counterM = 0;
            }
            $counterM++;
            for (my $i = 0; $i <= 10; $i++){
               if($allowed{$motdat[$i]} == 1){
                  if($motdat[$i] eq "A"){$chroms[$i][1] .= "C";}
                  if($motdat[$i] eq "B"){$chroms[$i][1] .= "T";}
               }
               else{$chroms[$i][1] .= "A";}
            }
            if($counterM == 400){
               for (my $i = 0; $i <= 10; $i++){
                  if($i != 3){
                     print OUTM "$chroms[$i][1]";
                  }
                  $chroms[$i][1] = "";
               }
               print OUTM "$band$band$band";
               $counterM = 0;
            }
         }
      }
   }
   close OUTF; close OUTM;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub LongSkittle(){

   print "   Combining Fasta files for child ";
   my $separator  = "A" x 4000;

   for (my $child = 1; $child <= 11; $child++){
      print " $child";
      for (my $p = 0; $p <= 1; $p++){
         my $parent;
         if ($p == 0){$parent = "F";}
         else        {$parent = "M";}
         my $file = "d:/perl64/programs/1463/haplotypes/$child all $parent.fa";
         open (OUT, ">$file") || die("Couldn't open >$file");
         print OUT ">child $child all $parent\n";
         for (my $chr = 1; $chr <= 21; $chr++){
            my $file = "d:/perl64/programs/1463/haplotypes/Child $child chr $chr $parent.fa";
            open (IN, "<$file") || die("Couldn't open <$file");
            my @seq = <IN>; chomp @seq; close IN;
            my $seq = $seq[1]; @seq = "";
            print OUT "$seq$separator\n";
         }
         close OUT;
      }
   }
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Privates(){

   print "   Loading data...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds = <IN>; chomp $inds;
   my @inds = split(/\,/, $inds);
   shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;
   print "   $numinds individuals in database\n";

   my $file = "d:/perl64/programs/1KG/Yvariants NoPrivates.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "$inds\n$pops\n$spops\n$shaps\n$haps\n";

   my %privates;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata;
      my $maj = shift @linedata;
      my $min = shift @linedata;
      my %results;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         $results{$linedata[$i]}++
      }
      if($results{0} == 1){
         for (my $i = 0; $i <= $numinds - 1; $i++){
            if($linedata[$i] == 0){
               $linedata[$i] = 1;
               my $data = join ("\,", @linedata);
               print OUT "$loc\,$maj\,$min\,$data\n";
               $privates{$inds[$i]}++;
               last;
            }
         }
      }
      if($results{1} == 1){
         for (my $i = 0; $i <= $numinds - 1; $i++){
            if($linedata[$i] == 1){
               $linedata[$i] = 0;
               my $data = join ("\,", @linedata);
               print OUT "$loc\,$maj\,$min\,$data\n";
               $privates{$inds[$i]}++;
               last;
            }
         }
      }
   }
   close IN;
   close OUT;
   my $file = "d:/perl64/programs/1KG/Privates.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Ind\,#Privates\n";
   foreach my $ind (@inds){
      print OUT "$ind\,$privates{$ind}\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub DistanceToConsensus(){

   my @shaplist = ("A", "B", "C", "D", "E", "G", "H", "I", "J", "L", "N", "O", "P", "Q", "R", "T", "");

   print "   Calculating...\n";
   my @consensus;
   my $file = "d:/perl64/programs/1KG/Consensus.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; shift @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      push (@consensus, $linedata[3]);
   }

   my $file = "d:/perl64/programs/1KG/YvariantsC.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds = <IN>; chomp $inds;
   my @inds = split(/\,/, $inds);
   shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my %shaps; my %spops;
   my @spops = split(/\,/, $spops); shift @spops; shift @spops; shift @spops;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      $shaps{$inds[$i]} = $shaps[$i];
      $spops{$inds[$i]} = $spops[$i];
   }

   my %distance; my %hdistance; my %pdistance;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata;
      my $maj = shift @linedata;
      my $min = shift @linedata;
      my $consensus = pop @linedata;
      my %results;
      if(substr($min, 0, 3) eq "<CN"){
         my %alleles;
         for (my $i = 0; $i <= $numinds - 1; $i++){
            $alleles{$linedata[$i]}++;
         }
         my $max;
         foreach my $allele (sort {$a cmp $b} keys %alleles){
            if ($alleles{$allele} > $max){
               $max = $alleles{$allele};
               $consensus = $allele;
            }
         }
      }
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if(substr($min, 0, 3) eq "<CN"){
            if($linedata[$i] != $consensus){
               $distance{$inds[$i]}{"consensus"}++;
               $distance{$inds[$i]}{"ancestor"}++;
               $hdistance{$shaps{$inds[$i]}}{"consensus"}++;
               $hdistance{$shaps{$inds[$i]}}{"ancestor"}++;
               $pdistance{$spops{$inds[$i]}}{"consensus"}++;
               $pdistance{$spops{$inds[$i]}}{"ancestor"}++;
            }
         }
         else{
            if($linedata[$i] != $consensus){
               $distance{$inds[$i]}{"consensus"}++;
               $hdistance{$shaps{$inds[$i]}}{"consensus"}++;
               $pdistance{$spops{$inds[$i]}}{"consensus"}++;
            }
            if($linedata[$i] != 0){
               $distance{$inds[$i]}{"ancestor"}++;
               $hdistance{$shaps{$inds[$i]}}{"ancestor"}++;
               $pdistance{$spops{$inds[$i]}}{"ancestor"}++;
            }
         }
      }
   }
   close IN;
   print "   Saving...\n";
   my $file = "d:/perl64/programs/1KG/Distance to Consensus.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Inds\,DConsensus\,DAncestor\n";
   for (my $i = 0; $i <= $numinds - 1; $i++){
      print OUT "$inds[$i]\,", $distance{$inds[$i]}{"consensus"}, "\,", $distance{$inds[$i]}{"ancestor"}, "\n";
   }

   print OUT "\n SHaps to Consensus\n";
   foreach my $shap (@shaplist){
      print OUT "$shap";
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if($shaps{$inds[$i]} eq $shap){
            print OUT "\,", $distance{$inds[$i]}{"consensus"};
         }
      }
      print OUT "\n";
   }

   print OUT "\n SHaps to Ancestor\n";
   foreach my $shap (@shaplist){
      print OUT "$shap";
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if($shaps{$inds[$i]} eq $shap){
            print OUT "\,", $distance{$inds[$i]}{"ancestor"};
         }
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub CountCategories(){

   my $file = "d:/perl64/programs/1KG/Yvariants.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds = <IN>; chomp $inds;
   my @inds = split(/\,/, $inds);
   shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;
   close IN;

   my %spops;
   my @spops = split(/\,/, $spops);
   my $numinds = @spops;
   foreach my $spop (@spops){
      $spops{$spop}++;
   }

   my %pops;
   my @pops = split(/\,/, $pops);
   foreach my $pop (@pops){
      $pops{$pop}++;
   }
   
   my %shaps;
   my @shaps = split(/\,/, $shaps);
   foreach my $shap (@shaps){
      $shaps{$shap}++;
   }

   my %haps;
   my @haps = split(/\,/, $haps);
   foreach my $hap (@haps){
      $haps{$hap}++;
   }
 
   my %spopshaps;
   for(my $i = 0; $i <= $numinds - 1; $i++){
      $spopshaps{$spops[$i]}{$shaps[$i]}++;
   }
   
   my $file = "d:/perl64/programs/1KG/CategoryCounts.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   
   print OUT "Spop\,Count\n";
   foreach my $spop (sort {$a cmp $b} keys %spops){
      print OUT "$spop\,$spops{$spop}\n";
   }
   print OUT "\n";

   print OUT "pop\,Count\n";
   foreach my $pop (sort {$a cmp $b} keys %pops){
      print OUT "$pop\,$pops{$pop}\n";
   }
   print OUT "\n";

   print OUT "SHap\,Count\n";
   foreach my $shap (sort {$a cmp $b} keys %shaps){
      print OUT "$shap\,$shaps{$shap}\n";
   }
   print OUT "\n";
   
   print OUT "Hap\,Count\n";
   foreach my $hap (sort {$a cmp $b} keys %haps){
      print OUT "$hap\,$haps{$hap}\n";
   }
   print OUT "\n";

   my @shaplist = ("A", "B", "C", "D", "E", "G", "H", "I", "J", "L", "N", "O", "P", "Q", "R", "T", "");
   my @spoplist = ("AFR", "AMR", "EAS", "EUR", "SAS", "");

   print OUT "Spop";
   foreach my $shap (@shaplist){
      print OUT "\,$shap";
   }
   print OUT "\n";
   foreach my $spop (@spoplist){
      print OUT "$spop";
      foreach my $shap (@shaplist){
         print OUT "\,$spopshaps{$spop}{$shap}";
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
sub MakeFasta(){

   print "   Loading Sequences...\n";

   my $file = "d:/perl64/programs/1KG/Yvariants.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;
   my @inds  = split(/\,/, $inds);  shift @inds;  shift @inds;  shift @inds;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   my $numinds = @inds;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($shaps[$i] eq ""){$shaps[$i] = "U";}
   }
 
   my %seqs;
   while(<IN>){
      my $line = $_; chomp $line;
      my @linedata = split(/\,/, $line); shift @linedata;
      my @alleles;
      $alleles[0] = shift @linedata;
      $alleles[1] = shift @linedata;
      if(substr($alleles[1], 0, 1) ne "<"){
         if(length($alleles[0]) == 1 && length($alleles[1]) == 1){
            my @counts;
            for (my $i = 0; $i <= $numinds - 1; $i++){
               $counts[$linedata[$i]]++;
            }
            if($counts[0] ne 1 && $counts[1] ne 1){
               for (my $i = 0; $i <= $numinds - 1; $i++){
                  $seqs{$shaps[$i]}{$inds[$i]} .= $alleles[$linedata[$i]];
               }
            }
         }
      }
   }

   print "   Saving...\n";
   my $file = "d:/perl64/programs/1KG/Seqs without privates.fas";
   open (OUT, ">$file") || die("Couldn't open >$file");
   foreach my $shap (sort {$a cmp $b} keys %seqs){
      foreach my $ind (sort {$a cmp $b} keys %{$seqs{$shap}}){
         print OUT ">$shap-$ind\n";
         my @seq = unpack("(A70)*", $seqs{$shap}{$ind});
         foreach my $string(@seq){
            print OUT "$string\n";
         }
      }
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
sub FindNearestNeighbor(){

   print "   Loading NN Table...\n";

   my @NNTable; my @seqnames;
   my $file = "d:/perl64/programs/1KG/Y NN Table.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $numinds = @data;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      my @linedata = split(/\,/, $data[$i]);
      $seqnames[$i] = shift (@linedata);
      for (my $j = $i + 1; $j <= $numinds - 1; $j++){
         $NNTable[$i][$j] = $linedata[$j];
         $NNTable[$j][$i] = $linedata[$j];
      }    
   }

   print "   Loading Shaps...\n";
   my %Shaps;
   my $file = "d:/perl64/programs/1KG/inddata.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $Shaps{$linedata[0]} = $linedata[4];
   }

   my @nearest;
   print "   Calculating...\n";
   for (my $i = 0; $i <= $numinds - 1; $i++){
      my $min = 1_000_000;
      my $nearest;
      for (my $j = 0; $j <= $numinds - 1; $j++){
         if($NNTable[$i][$j] > 0 && $NNTable[$i][$j] < $min){
            $min = $NNTable[$i][$j];
            $nearest = "$seqnames[$j]\,$Shaps{$seqnames[$j]}\,$min";
         }
      }
      push(@nearest, "$seqnames[$i]\,$Shaps{$seqnames[$i]}\,$nearest");
   }
   
   print "   Saving...\n";
   my $file = "d:/perl64/programs/1KG/NNeighbors.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   foreach my $line (@nearest){
      print OUT "$line\n";
   }
   close OUT;
}

#Changed

#HG04158 is closest to an R, not farthest away, changed
#HG04185 is closest to an H, not farthest away, changed
#NA12413 is closest to an I, not farthest away, changed
#HG03890 is closest to a  J, not farthest away, changed

#HG00628 is closest to a  C, but it farther away than any other C-C pair
#HG02040 is closest to a  J, but it farther away than any other J-J pair
#HG02684 is closest to a  J, but it farther away than any other J-J pair
#HG03594 is closest to an H, but it farther away than any other H-H pair

#HG03680 is closest to another unknown (HG03837)
#HG03837 is closest to another unknown (HG03680)

#HG03848 is closest to another unknown (HG03965)
#HG03965 is closest to another unknown (HG03848)

#HG03870 is closest to another unknown (HG03872)
#HG03872 is closest to another unknown (HG03870)
#HG04238 is closest to another unknown (HG03872)
#HG04033 is closest to another unknown (HG03872)
#HG03792 is closest to another unknown (HG04033)

# Consensus is different from ref at 376 SNPs

# B/c P is so recent and b/c it nests within R and b/c R has two main branches, 
# only one of which gave rise to P, rolled haplogroup P in to haplogroup R

# split off R1a1 b/c they were so different from the other Rs, labeled 'V'

# split "A" into 4 groups
# HG02645 = group 1, labeled 'W'
# HG02666 = group 1
# HG02613 = group 2, labeled 'X'
# HG01890 = group 3, labeled 'Y'
# HG02982 = group 4, labeled 'Z'

#Clashes
#0-0
#1-231
#2-561
#3-322
#4-1819
#5-336
#6-585
#7-525
#8-956
#9-321
#10-359
#11-1405
#12-152
#13-594
#14-33
#15-371
#16-127
#17-667
#18-741
#19-1146
#20-21
#21-0
#22-0
#23-0
#24-298
#25-540
#26-2464
#27-1200
#28-1408
#29-1205
#30-2117
#31-453
#32-1541
#33-1829
#34-2327
#35-0

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub RDistanceTable(){

   print "   Loading Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds    = split(/\,/, $inds);  shift @inds;  shift @inds;  shift @inds; pop @inds;
   my @pops    = split(/\,/, $pops);  shift @pops;  shift @pops;  shift @pops;
   my @spops   = split(/\,/, $spops); shift @spops; shift @spops; shift @spops;
   my @shaps   = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   my @haps    = split(/\,/, $haps);  shift @haps;  shift @haps;  shift @haps;
   my $numinds = @inds;

   my @iinds; my @ipops; my @ispops; my @ishaps; my @ihaps; my $icounter = 0;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($shaps[$i] eq "P" || $shaps[$i] eq "R" || $shaps[$i] eq "V"){
         $iinds[$icounter]  = $inds[$i];
         $ipops[$icounter]  = $pops[$i];
         $ispops[$icounter] = $spops[$i];
         $ishaps[$icounter] = $shaps[$i];
         $ihaps[$icounter]  = $haps[$i];
         $icounter++;
      }
   }

   my @sequences; my $counter = 0; my @seqlist;
   while(<IN>){
      my $line = $_; chomp $line;
      my @linedata = split(/\,/, $line); shift @linedata;
      my @alleles; $alleles[0] = shift @linedata; $alleles[1] = shift @linedata;
      my $icounter = 0;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if($shaps[$i] eq "P" || $shaps[$i] eq "R" || $shaps[$i] eq "V"){
            $sequences[$icounter][$counter] = $linedata[$i];
            $icounter++;
         }
      }
      $counter++;
   }
   close IN;

   my @distances;
   print "   Calculating Distances...";
   for (my $i = 0; $i <= $icounter - 1; $i++){
      print " $i";
      for (my $j = $i + 1; $j <= $icounter - 1; $j++){
         for (my $k = 0; $k <= $counter - 1; $k++){
            if($sequences[$i][$k] != $sequences[$j][$k] && $sequences[$i][$k] ne "" && $sequences[$j][$k] ne ""){
               $distances[$i][$j]++;
               $distances[$j][$i]++;
            }
         }
      }
   }
   
   print "\n   Saving...\n";
   my $file = "d:/perl64/programs/1KG/RDistanceTable.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Ind\,Pop\,Spop\,Shap\,Hap";
   for (my $i = 0; $i <= $icounter - 1; $i++){ print OUT "\,", $iinds[$i]; }
   print OUT "\n"; print OUT "\,\,\,\,";
   for (my $i = 0; $i <= $icounter - 1; $i++){ print OUT "\,", $ipops[$i]; }
   print OUT "\n"; print OUT "\,\,\,\,";
   for (my $i = 0; $i <= $icounter - 1; $i++){ print OUT "\,", $ispops[$i]; }
   print OUT "\n"; print OUT "\,\,\,\,";
   for (my $i = 0; $i <= $icounter - 1; $i++){ print OUT "\,", $ishaps[$i]; }
   print OUT "\n"; print OUT "\,\,\,\,";
   for (my $i = 0; $i <= $icounter - 1; $i++){ print OUT "\,", $ihaps[$i]; }
   print OUT "\n";
   for (my $i = 0; $i <= $icounter - 1; $i++){
      print OUT "$iinds[$i]\,$ipops[$i]\,$ispops[$i]\,$ishaps[$i]\,$ihaps[$i]";
      for (my $j = 0; $j <= $icounter - 1; $j++){
         print OUT "\,$distances[$i][$j]";
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Clashes(){

   print "   Loading Group Designations...\n";
   my %list; my @groupnames; my %shapnum;
   my $file = "d:/perl64/programs/1KG/Groups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $counter = 0; my @groupcounts;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      my $shapnames = join('', @linedata);
      $shapnum{$shapnames} = $groupnum;
      foreach my $shap (@linedata){
         $list{$shap}{$groupnum} = 1;
         $groupnames[$counter] .= $shap;
         $groupcounts[$groupnum]++;
      }
      $counter++;
   }
   my $numgroups = @groupnames;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds - 1; # -1 b/c consensus is appended to the end of the list

   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $shaps[$i];
   }

   my @alleles; my $counter = 0;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      shift @linedata; shift @linedata; shift @linedata; pop @linedata;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         foreach my $group (sort {$a<=>$b} keys %{$list{$shaps[$i]}}){
            $alleles[$group][$counter][$linedata[$i]]++;
         }
      }
      $counter++;
   }
   close IN;

   print "   Calculating Hets...\n";   
   my @hets;
   for (my $group = 0; $group <= $numgroups - 1; $group++){
      for (my $i = 0; $i <= $counter - 1; $i++){
         if($alleles[$group][$i][0] > 0 && $alleles[$group][$i][1] > 0){
            $hets[$group][$i] = 1;
         }
         else{
            $hets[$group][$i] = 0;
         }
      }
   }   

   print "   Calculating Clashes...";
   my @clashes;
   for (my $shap1 = 0; $shap1 <= $numgroups - 1; $shap1++){
      print " $groupnames[$shap1]";
      for (my $shap2 = 0; $shap2 <= 24; $shap2++){
         if($shap1 != $shap2 && $list{$groupnames[$shap2]}{$shap1} != 1){
            for (my $i = 0; $i <= $counter - 1; $i++){
               if($hets[$shap1][$i] == 1 && $hets[$shap2][$i] == 1){
                  $clashes[$shap1][$shap2]++;
               }
            }
         }
      }
   }
   print "\n";

   print "Saving...\n";
   my $file = "d:/perl64/programs/1KG/Clashes.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");

   for (my $shap1 = 1; $shap1 <= $numgroups - 1; $shap1++){
      print OUT "\,$groupnames[$shap1]";
   }
   print OUT "\n";
   for (my $shap1 = 1; $shap1 <= $numgroups - 1; $shap1++){
      print OUT "$groupnames[$shap1]";
      for (my $shap2 = 1; $shap2 <= $numgroups - 1; $shap2++){
         print OUT "\,$clashes[$shap1][$shap2]";
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub NandD(){

   print "   Loading Distance Data...\n";
   my $file = "d:/perl64/programs/1KG/Distance to Nodes CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data; close IN;
   my $shaps = shift @data;
   my @shaps = split(/\,/, $shaps);
   my %shapcols;
   for(my $i = 0; $i <= 30; $i++){
      $shapcols{$shaps[$i]} = $i;
   }

   my %Inds;
   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @rdata = <IN>; chomp @rdata; close IN;
   foreach my $line (@rdata){
      my @linedata = split(/\,/, $line);
      $Inds{$linedata[0]}{"shap"} = $linedata[3];
   }

   shift @data;
   foreach my $line(@data){
      my @linedata = split(/\,/, $line);
      if($Inds{$linedata[0]}{"shap"} eq ""){
         $Inds{$linedata[0]}{"shap"} = $linedata[1];
      }
      $Inds{$linedata[0]}{"distance"} = $linedata[ $shapcols{ $Inds{$linedata[0]}{"shap"} } ];
   }

   print "   Calculating Ns...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>;
   my $spops = <IN>;
   my $shaps  = <IN>;
   my $haps  = <IN>;
   my @inds = split(/\,/, $inds);
   shift @inds; shift @inds; shift @inds; pop @inds;
   my $numinds = @inds;

   while(<IN>){
      my $line = $_; chomp $line;
      my @linedata = split(/\,/, $line); shift @linedata; shift @linedata; shift @linedata;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if($linedata[$i] eq "."){
            $Inds{$inds[$i]}{"missing"}++;
         }
      }
   }
   close IN;

   print "Saving...\n";
   my $file = "d:/perl64/programs/1KG/NandDs.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Ind\,Shap\,N\,D\n";
   foreach my $ind (sort {$a cmp $b} keys %Inds){
      print OUT "$ind\,", $Inds{$ind}{"shap"}, "\,", $Inds{$ind}{"missing"}, "\,", $Inds{$ind}{"distance"}, "\n";
   }
   close OUT;

}
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
sub InternalConsensi(){

   print "   Loading Group Designations...\n";
   my %list; my @groupnames;
   my $file = "d:/perl64/programs/1KG/Groups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $counter = 0;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      foreach my $group (@linedata){
         $list{$group}{$groupnum}++;
         $groupnames[$counter] .= $group;
      }
      $counter++;
   }
   my $numgroups = @groupnames;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading and Sorting Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds - 1;

   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]};
      }
      $shaps{$inds[$i]} = $shaps[$i];
   }

   my %consensuscounts; my %outs; my $counter = 0;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata;
      my @alleles;
      my $maj = shift @linedata; $alleles[0] = $maj;
      my $min = shift @linedata; $alleles[1] = $min;
      my $consensus = pop @linedata;
      $outs{"MainConsensus"}{$counter} = $consensus;
      $outs{"Major"}{$counter} = $maj;
      $outs{"Minor"}{$counter} = $min;
      for (my $group = 0; $group <= $numgroups - 1; $group++){
         $consensuscounts{$group}{$counter}{0} = 0;
         $consensuscounts{$group}{$counter}{1} = 0;
      }
      for (my $i = 0; $i <= $numinds - 1; $i++){
         my $shap = $shaps[$i];
         if($linedata[$i] == 0 || $linedata[$i] == 1){
            foreach my $listelement (sort {$a<=>$b} keys %{$list{$shap}}){
               $consensuscounts{$listelement}{$counter}{$linedata[$i]}++;
            }
         }
      }
      $counter++;
   }

   print "   Calculating...\n";
   my %consensus;
   for (my $i = 0; $i <= $counter - 1; $i++){
      for (my $group = 0; $group <= $numgroups - 1; $group++){
         if($consensuscounts{$group}{$i}{0} > $consensuscounts{$group}{$i}{1}){
            $consensus{$group}{$i} = 0;
         }
         else{
            $consensus{$group}{$i} = 1;
         }
         if($consensuscounts{$group}{$i}{0} == $consensuscounts{$group}{$i}{1}){
            $consensus{$group}{$i} = 0;
         }
      }
   }

   print "   Saving Consensi...\n";
   my $file = "d:/perl64/programs/1KG/Consensi CP.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "\,\,\,";
   foreach my $group(sort {$a<=>$b} keys %consensus){
      print OUT "\,";
      foreach my $shap(@shaps){
         if($list{$shap}{$group} > 0){
            print OUT "$shap";
         }
      }
   }
   print OUT "\n";
   print OUT "i\,Main\,Major\,Minor";
   foreach my $group(sort {$a<=>$b} keys %consensus){
      print OUT "\,$group";
   }
   print OUT "\n";
   for (my $i = 0; $i <= $counter - 1; $i++){
      print OUT "$i";
      print OUT "\,", $outs{"MainConsensus"}{$i};
      print OUT "\,", $outs{"Major"}{$i};
      print OUT "\,", $outs{"Minor"}{$i};
      foreach my $group(sort {$a cmp $b} keys %consensus){
         print OUT "\,$consensus{$group}{$i}";
      }
      print OUT "\n";
   }
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Percents(){

   my %allowed; $allowed{0} = 1; $allowed{1} = 1; $allowed{2} = 1; $allowed{3} = 1; $allowed{4} = 1; $allowed{5} = 1;

   print "   Loading Group Designations...\n";
   my $file = "d:/perl64/programs/1KG/Groups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $counter = 0; my %list; my @groupnames; my %shapsingroup;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      foreach my $group (@linedata){
         $list{$group}{$groupnum}++;
         $groupnames[$counter] .= $group;
         $shapsingroup{$groupnum}{$group}++;
      }
      $counter++;
   }
   my $numgroups = @groupnames;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading Sequence Data...\n";

   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;
   my @shaps  = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   my @inds   = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;
   my %shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $shaps[$i];
   }

   my @inddata; my $counter = 0;
   while(<IN>){
      my $line = $_; chomp $line;
      my @linedata = split(/\,/, $line); shift @linedata; shift @linedata; shift @linedata;
      for (my $ind = 0; $ind <= $numinds - 1; $ind++){
         $inddata[$ind][$counter] = $linedata[$ind];
         $counter++;
      }
   }

   print "   Loading Nodes...\n";

   my @shapdata; $counter = 0;
   my $file = "d:/perl64/programs/1KG/nodes cp2.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $shaps = <IN>; chomp $shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   my $numshaps = @shaps;
   while(<IN>){
      my $line = $_; chomp $line;
      my @linedata = split(/\,/, $line); shift @linedata; shift @linedata; shift @linedata;
      for (my $shap = 0; $shap <= $numshaps - 1; $shap++){
         $shapdata[$shap][$counter] = $linedata[$shap];
      }
      $counter++;
   }
   
   print "   Calculating...";
   my %counts;
   for(my $ind = 0; $ind <= $numinds - 1; $ind++){
      print "$inds[$ind] ";
      for(my $shap = 0; $shap <= $numshaps - 1; $shap++){
         for(my $i = 0; $i <= $counter - 1; $i++){
            if($allowed{$inddata[$ind][$i]} == 1 && $allowed{$shapdata[$shap][$i]} == 1){
               $counts{$ind}{$shap}++;
            }
         }
      }
      $inddata[$ind] = "";
   }
   print "\n";

   print "   Saving...\n";
   my $file = "d:/perl64/programs/1KG/Y Percents.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");

   for (my $j = 0; $j <= $numshaps - 1; $j++){
      print OUT "$j\,\,$groupnames[$j]\n";
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if(index ($groupnames[$j], $shaps{$inds[$i]}) >= 0){
            print OUT "$inds[$i]\,$counts{$i}{$j}\,$shaps{$inds[$i]}\n";
         }
      }
      print OUT "\n";
   }
   close OUT;
}


sub condenseclashes(){

   my $file = "d:/perl64/programs/1KG/megaGroups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my %mega; my @allshapsingroup;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      $allshapsingroup[$groupnum] = join('', @linedata);
      foreach my $shap (@linedata){
         $mega{$shap} = $groupnum;
      }
   }
   $mega{"F"} = 99;
   $mega{"M"} = 99;

   my $suffix = "CP2NSH.csv";
   my $file = "d:/perl64/programs/1KG/Clash data mega $suffix";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   shift @data; shift @data; shift @data;
   my $shaps = shift @data;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps; pop @shaps;
   my $numinds = @shaps;
   shift @data;
   my %locdata;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata; shift @linedata; shift @linedata;
      for(my $i = 0; $i <= $numinds - 1; $i++){
         $locdata{$mega{$shaps[$i]}}{$loc}{$linedata[$i]}++;
      }
   }
   my $file = "d:/perl64/programs/1KG/Clash data summary mega $suffix";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "GNum\,Shaps";
   foreach my $loc (sort {$a<=>$b} keys %{$locdata{1}}){
      print OUT "\,$loc";
   }
   print OUT "\n";
   foreach (my $mega = 0; $mega <= 8; $mega++){
      print OUT "$mega\,$allshapsingroup[$mega]";
      foreach my $loc (sort {$a<=>$b} keys %{$locdata{1}}){
         print OUT "\,", $locdata{$mega}{$loc}{"0"}, ".", $locdata{$mega}{$loc}{"1"}, ".", $locdata{$mega}{$loc}{"."};
      }
      print OUT "\n";
   }
   close OUT;
}


# deal with singletons/privates (requires more than one sequence in each group or group will be erased)
#               if    ($alleles[0][0] == 1){$groupnodes[$group][$i] = 1;}
#               elsif ($alleles[0][1] == 1){$groupnodes[$group][$i] = 0;}
#               # test to see if heterozygosity restricted to a single group
#               elsif ($numnamesingroup[$group] > 1){ # this works for combination groups only, so will not be applied to single groups
#                  my $heteros; my $allele;
#                  for (my $ingroup = 1; $ingroup <= 21; $ingroup++){
#                     if($list{$groupnames{$ingroup}}{$group} == 1){ # the group is within the combination group
#                        # mutation must have occurred within that group, node = homozygous allele
#                        if($groupalleles[$ingroup][$i][0] > 0 && $groupalleles[$ingroup][$i][0] > 0){
#                           $heteros++;
#                        }
#                        else{
#                           if($groupalleles[$ingroup][$i][0] == 0){$allele = 1;}
#                           else{$allele = 0;}
#                        }
#                     }
#                  }
#                  if($heteros == 1){
#                     $groupnodes[$group][$i] = $allele;                  
#                  }
#                  else{$groupnodes[$group][$i] = ".";}
#               }
#               else{

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub DatatoFASTA(){

   my $suffix = "mega CP2NSH.csv";

   print "   Loading Group Designations...\n";
   my $file = "d:/perl64/programs/1KG/megaGroups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my %mega;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      foreach my $shap (@linedata){
         $mega{$shap} = $groupnum;
      }
   }
   $mega{"F"} = 99;
   $mega{"M"} = 99;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;

   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $mega{$shaps[$i]};
   }

   my @sequences;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      shift @linedata;
      my %alleles;
      $alleles{0} = shift @linedata;
      $alleles{1} = shift @linedata;
      $alleles{"."} = ".";
      for (my $i = 0; $i <= $numinds - 1; $i++){
         $sequences[$i] .= $alleles{$linedata[$i]};
      }
   }
   close IN;
   my $numseqs = @sequences;

   my $file = "d:/perl64/programs/1KG/YvariantsCP.FAS";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for (my $i = 0; $i <= $numseqs - 1; $i++){
      print OUT ">$shaps[$i]-$inds[$i]\n$sequences[$i]\n\n";
   }
   close OUT;

   my $file = "d:/perl64/programs/1KG/names.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for (my $i = 0; $i <= $numseqs - 1; $i++){
      print OUT "#$shaps[$i]-$inds[$i]\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub TestNeanderthalVCF(){
   my $file = "d:/perl64/programs/1KG/AltaiNea.hg19_1000g.Y.mod.vcf/AltaiNea.hg19_1000g.Y.mod.vcf";
   open (IN, "<$file") || die("Couldn't open <$file");
   print "Pos\tID\tRef\tAlt";
   my $file = "d:/perl64/programs/1KG/AltaiNea.hg19_1000g.Y.mod.vcf/AltaiNea.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Pos\,Ref\,Alt\,Allele1\,Allele2\n";

   while(<IN>){
      my $line = $_;
      chomp $line;
      #print $line;
      my @linedata = split(/\t/, $line);
      if($linedata[0] eq "Y"){
         print OUT "$linedata[1]\,$linedata[3]\,$linedata[4]\,";
         my @tempdata  = split(/\:/, $linedata[8]);
         my $GTpos;
         my $numthings = @tempdata;
         for(my $i = 0; $i <= $numthings - 1; $i++){
            if($tempdata[$i] eq "GT"){
               $GTpos = $i;
               last;
            }
         }
         #print "$linedata[8] $linedata[9] ";
         if($GTpos >= 0){
            my @tempdata  = split(/\:/, $linedata[9]);
            my @tempdata2 = split(/\//, $tempdata[$GTpos]);
            print OUT "$tempdata2[0]\,$tempdata2[1]\n";
            #print "$linedata[1]\,$linedata[3]\,$linedata[4]\,$tempdata2[0] $tempdata2[1]\n";
         }
         #foreach my $datum (@linedata){print "\n   $datum";}
      }
      #print $line;
   }
   close OUT;
   close IN;
}
#A = # A on forward and reverse strands
#C
#G
#T
#IR = # reads with indel starting at this pos
#AD = allelic depth
#DP = read depth
#GQ = genotype quality
#GT = genotype
#PL = phred likelihood
#AF = alt alleles/total alleles in 1KG
#AMR_AF = alt allele freq in AMR
#ASN_AF
#AFR_AF
#EUR_AF
#1000gALT
#TS = within EPO Compara 6 primate block
#TSSeq = sequence of above
#CAnc = ref-chimp/human sequence
#GAnc = ref-gorilla/human sequence
#OAnc = ref-orang/human sequence
#mSC = mamallian conservation score
#pSG = primate conservation score
#GRP = GERP conservation score
#bSC = B score
#Map20 = Duke score something
#CPG = CPG
#RM = repeat masked in EPO Compara 6 primate block
#SysErr = duh
#SysHCB = human-chimp-bonbo
#UR = copy number control region
#AC = allele count in genotypes for each alt allele
#AF = allele frequency for each alt allele
#AN = number of alleles
#BaseQRankSum = ignore
#DP = filtered depth
#DS = were any samples downsampled?
#Dels = deletions present
#FS = Fisher test strand bias
#Hrun = homopolymer run
#HaplotypeScore = 
#InbreedingCoeff = HW
#MQ = RMS mapping quality
#MQ0 = total mapping quality zero reads
#MQRankSum = Wilcoxian Z-score
#QD = quality by depth
#ReadPosRankSum = Wilcoxian something

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub ReduceFASTA(){

   my $suffix = "mega CP2NSH, no privates.csv";

   print "   Opening stripped FASTA\n";
   my $file = "d:/perl64/programs/1KG/FASTA stripped $suffix";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN;
   my $sequences = join ('', @data);
   my @sequencelist = split(/\>/, $sequences); shift @sequencelist;
   my @sequences;
   my $seqnum = 0;
   my @seqnames;
   my $numnucs;
   foreach my $sequence (@sequencelist){
      my @seqdata = split(/\n/, $sequence);
      my $name = shift @seqdata;
      $seqnames[$seqnum] = $name;
      my $sequence = shift @seqdata;
      my @sequence = split(//, $sequence);
      $numnucs = @sequence;
      for(my $i = 0; $i <= $numnucs - 1; $i++){
         $sequences[$seqnum][$i] = $sequence[$i];
      }
      $seqnum++;
   } 
   
   print "   Calculating";
   my $numinds = $seqnum; my $counter = $numnucs;
   for (my $ind = 0; $ind <= $numinds - 1; $ind++){
      for (my $i = 0; $i <= $counter - 1; $i++){
         if($sequences[$ind][$i] eq "."){$sequences[$ind][$i] = "-";}
      }
   }
   print ".";
   for (my $ind = 1; $ind <= $numinds - 1; $ind++){
      for (my $i = 0; $i <= $counter - 1; $i++){
         if($sequences[0][$i] ne "-"){
            if($sequences[$ind][$i] eq $sequences[0][$i]){$sequences[$ind][$i] = ".";}
         }
      }
   }
   print ".";
   my @chrs; 
   for (my $ind = 0; $ind <= $numinds - 1; $ind++){
      for (my $i = 0; $i <= $counter - 1; $i++){
         $chrs[$ind] .= $sequences[$ind][$i];
      }
   }
   print "\n   Saving\n";
   my $file = "d:/perl64/programs/1KG/FASTA stripped and reduced $suffix";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for (my $ind = 0; $ind <= $numinds - 1; $ind++){
      print OUT ">$seqnames[$ind]\n$chrs[$ind]\n\n";
   }
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub AnalyzeAncestorFile(){

   print "   Loading\n";
   my $file = "d:/perl64/programs/1KG/ancestralgroupings.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN;
   my $numseqs = @data;
   my @seqdata;
   my @ancestordata;
   for(my $i = 0; $i <= 84 - 1; $i++){
      # 1. 8-HG00096
      my @linedata = split(/\./, $data[$i]);
      if(substr($linedata[1], 1, 1) eq "("){last;}
      else{$seqdata[$linedata[0]] = substr($linedata[1], 1, 1);}
   }
   for(my $i = 85 - 1; $i <= 164 - 1; $i++){
      # 84. (74 . 82)
      my @linedata = split(/\./, $data[$i]);
      my $first = $seqdata[substr($linedata[1], 2)];
      my $second = $seqdata[substr($linedata[2], 1, -1)];
      $seqdata[$linedata[0]] = $first.$second; #)    1234. (235 . 1172)
   }
   my $file = "d:/perl64/programs/1KG/ancestralgroupingssummary.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for(my $i = 0; $i <= 165 - 1; $i++){
      print OUT "$i\,A$seqdata[$i]\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub CheckYAlleles(){

   my $suffix = ".csv";

   print "   Loading Group Designations...\n";
   my $file = "d:/perl64/programs/1KG/Megagroups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $numgroups = @data;
   my @allshapsingroup; my %shapsingroup;
   my %mega;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      $allshapsingroup[$groupnum] = join('', @linedata);
      foreach my $shap (@linedata){
         $mega{$shap} = $groupnum;
      }
   }

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;

   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $mega{$shaps[$i]};
   }

   my @alleles; my @locs; my @majors; my @minors; my $counter = 0; my @missing; my @groupnodes; my @groupalleles; my @allalleles;
   while(<IN>){
      my $line = $_; chomp $line;
      my @linedata = split(/\,/, $line);
      $locs[$counter] = shift @linedata;
      my $major = shift @linedata;
      my $minor = shift @linedata;
      $majors[$counter] = $major;
      $allalleles[$counter][0] = $major;
      my @minors = split(/\;/, $minor);
      my $numminors = @minors;
      for (my $i = 0; $i <= $numminors - 1; $i++){
         $allalleles[$counter][$i + 1] = $minors[$i];
      }
      for (my $i = 0; $i <= $numinds - 1; $i++){
         my $shap = $shaps[$i];
         $groupalleles[$mega{$shap}][$counter][$linedata[$i]]++;
         $groupalleles[16][$counter][$linedata[$i]]++;
      }
      $counter++;
   }
   close IN;
   my $file = "d:/perl64/programs/1KG/YAlleles.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   for (my $i = 0; $i <= $counter; $i++){
      print OUT "$locs[$i]";
      for (my $j = 0; $j <= 6; $j++){
         print OUT "\,$allalleles[$i][$j]";
      }
      for (my $k = 0; $k <= 16; $k++){
         for (my $j = 0; $j <= 6; $j++){
            print OUT "\,$groupalleles[$k][$i][$j]";
         }
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Consensii(){

   my $suffix = "mega CP2NSH.csv";

   print "   Loading Group Designations...\n";
   my $file = "d:/perl64/programs/1KG/megaGroups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $numgroups = @data;
   my @allshapsingroup; my %shapsingroup;
   my %mega;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      $allshapsingroup[$groupnum] = join('', @linedata);
      foreach my $shap (@linedata){
         $mega{$shap} = $groupnum;
      }
   }
   $mega{"F"} = 99;
   $mega{"M"} = 99;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;

   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   my @shapcounts;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $mega{$shaps[$i]};
      $shapcounts[$mega{$shaps[$i]}]++;
   }

   print "Group sizes:\n";
   for(my $group = 0; $group <= 9; $group++){
      print "   $group: $shapcounts[$group]\n";
   }
   <STDIN>;

   my @sequences; my @locs; my @majors; my @minors; my $counter = 0; my @missing; my @groupnodes; my @groupalleles; my @allalleles;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      $locs[$counter] = shift @linedata;
      my $major = shift @linedata;
      my $minor = shift @linedata;
      $majors[$counter] = $major;
      $minors[$counter] = $minor;
      $allalleles[$counter][0] = $major;
      $allalleles[$counter][1] = $minor;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         $sequences[$i][$counter] = $linedata[$i];
         my $shap = $shaps[$i];
         $groupalleles[$mega{$shap}][$counter][$linedata[$i]]++;
         $groupalleles[9][$counter][$linedata[$i]]++;
      }
      $counter++;
   }
   close IN;
   my $file = "d:/perl64/programs/1KG/Consensii $suffix";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Loc\,Ref\,Alt";
   for(my $group = 0; $group <= 8; $group++){
      print OUT "\,$group";
   }
   print OUT "\n";
   for (my $i = 0; $i <= $counter; $i++){
      print OUT "$locs[$i]\,$majors[$i]\,$minors[$i]";
      for(my $group = 0; $group <= 9; $group++){
         if($groupalleles[$group][$i][0] >= $groupalleles[$group][$i][1] ){
            print OUT "\,0";
         }
         else{print OUT "\,1";}
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub ISOGGType(){

   print "   Loading ISOGG data...\n";
   my $file = "d:/perl64/programs/1KG/ISOGG2015 SNPs.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data; shift @data;
   my %ISOGG; my %groups;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $ISOGG{$linedata[4]}{"exists"} = 1;
      $ISOGG{$linedata[4]}{"name"}   = $linedata[0];
      $ISOGG{$linedata[4]}{"group"}  = $linedata[1];
      $ISOGG{$linedata[4]}{"rsID"}   = $linedata[3];
      $ISOGG{$linedata[4]}{"from"}   = $linedata[5];
      $ISOGG{$linedata[4]}{"to"}     = $linedata[6];
      $groups{$linedata[1]}++;
   }

   print "   Loading Group Designations...\n";
   my $file = "d:/perl64/programs/1KG/megaGroups CP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   my $numgroups = @data;
   my @allshapsingroup; my %shapsingroup;
   my %mega;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $groupnum = shift @linedata;
      $allshapsingroup[$groupnum] = join('', @linedata);
      foreach my $shap (@linedata){
         $mega{$shap} = $groupnum;
      }
   }
   $mega{"F"} = 99;
   $mega{"M"} = 99;

   print "   Loading Recodes...\n";
   my %recodes;
   my $file = "d:/perl64/programs/1KG/Recode.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; chomp @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recodes{$linedata[0]} = $linedata[3];
   }

   print "   Loading Sequences...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @inds = split(/\,/, $inds); shift @inds; shift @inds; shift @inds;
   my $numinds = @inds;

   my %shaps;
   my @shaps = split(/\,/, $shaps); shift @shaps; shift @shaps; shift @shaps;
   for (my $i = 0; $i <= $numinds - 1; $i++){
      if($recodes{$inds[$i]} ne ""){
         $shaps[$i] = $recodes{$inds[$i]}
      }
      $shaps{$inds[$i]} = $mega{$shaps[$i]};
   }

   my %TypedInds; my @locs; my $counter = 0;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      $locs[$counter] = shift @linedata;
      if($ISOGG{$locs[$counter]}{"exists"} == 1){
         my @alleles; $alleles[0] = shift @linedata; $alleles[1] = shift @linedata;
         for (my $i = 0; $i <= $numinds - 1; $i++){
            my $test = $alleles[$linedata[$i]];
            if($test eq $ISOGG{$locs[$counter]}{"to"}){
               $TypedInds{$i}{$ISOGG{$locs[$counter]}{"group"}}++;
            }
         }
      }
   }
   close IN;
   my $file = "d:/perl64/programs/1KG/ISOGGtyped.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Ind\,Mega\,Shap";
   foreach my $group (sort {$a cmp $b} keys %groups){
      print OUT "\,$group($groups{$group})";
   }
   print OUT "\n";
   for (my $i = 0; $i <= $numinds - 1; $i++){
      print OUT "$inds[$i]\,$mega{$inds[$i]}\,$shaps{$inds[$i]}";
      foreach my $group (sort {$a cmp $b} keys %groups){
         print OUT "\,", $TypedInds{$i}{$group};
      }
      print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Scozzari(){

   print "   Loading Scozzari data...\n";
   my $file = "d:/perl64/programs/1KG/Scozzari et al. ancestral data.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data; shift @data;
   my %Scozzari;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $Scozzari{$linedata[0]}{"exists"}    = 1;
      $Scozzari{$linedata[0]}{"data"}      = "$linedata[2]$linedata[3]$linedata[4]$linedata[5]";
      $Scozzari{$linedata[0]}{"ancestral"}  = $linedata[4];
   }

   my $file = "d:/perl64/programs/1KG/Y Nodes Final Edits Scozz.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");

   print "   Loading Nodes...\n";
   my $file = "d:/perl64/programs/1KG/Y Nodes Final Edits.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $headers = shift @data;
   print OUT "ScozzDat\,AncestralAllele\,$headers\n";
   my %Nodes;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      if($Scozzari{$linedata[0]}{"exists"} == 1){
         print OUT $Scozzari{$linedata[0]}{"data"};
         print OUT "\,", $Scozzari{$linedata[0]}{"ancestral"};
         print OUT "\,$line\n";
      }
      else{
         print OUT "\,\,$line\n";
      }
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub TypeSNPs(){

   print "   Loading ISOGG data...\n";
   my $file = "d:/perl64/programs/1KG/ISOGG2015SNPs.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my %ISOGG; # Pos SNP Haplogroup AltSNP rsID Mutation AltMutation AltNuc RefNuc
   my $ISOGGheaders = shift @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $ISOGG{$linedata[0]} = $line;
   }
   
   my $file = "d:/perl64/programs/1KG/ISOGG-Typed 1KG SNPs.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   
   print "   Loading and parsing 1KG data...\n";
   my $file = "d:/perl64/programs/1KG/Nodes v2 mega CP2NSH A0A1Ances.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   # Loc	Ref	Alt	SAnc	HAnc	YZ	WX	B	DE	C	GH	IJ	LT	NO	QR	Num1s	YZ	WX	B	DE	C	GH	IJ	LT	NO	QR
   my $TKGheaders = shift @data;
   print OUT "$TKGheaders";
   print OUT "\,YZ\,WX\,B\,DE\,C\,GH\,IJ\,LT\,NO\,QR";
   print OUT "\,YZ\,WX\,B\,DE\,C\,GH\,IJ\,LT\,NO\,QR";
   print OUT "$ISOGGheaders";
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      my $loc = $linedata[0];
      my @datums;
      for (my $i = 16; $i <=25; $i++){
         my @datum = split(/_/, $linedata[$i]);
         $datums[$i - 16][0] = $datum[0];
         $datums[$i - 16][1] = $datum[1];
      }
      print OUT $line;
      for(my $i = 0; $i <= 9; $i++){ print OUT "\,$datums[$i][0]"; }
      for (my $i = 0; $i <= 9; $i++){ print OUT "\,$datums[$i][1]"; }
      if($ISOGG{$loc} ne ""){
         print OUT "\,$ISOGG{$loc}";
         delete $ISOGG{$loc};
      }
      print OUT "\n";
   }
   print OUT "\nUnused ISOGG data\n";
   foreach my $ISOGG (sort {$a<=>$b} keys %ISOGG){
      print OUT "$ISOGG{$ISOGG}\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Nearest(){

   my $file = "d:/perl64/programs/1KG/Y NN Table CP2.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $headers = shift @data;
   my @inds = split (/\,/, $headers);
   my $numinds = @inds; my $counter = 0; my %NNs;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      for (my $i = $counter; $i <= $numinds - 1; $i++){
         $NNs{$i}{$counter} = $linedata[$i];
         $NNs{$counter}{$i} = $linedata[$i];
      }
      $counter++;
   }
   my $file = "d:/perl64/programs/1KG/Nearest.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   foreach my $ind (sort {$a cmp $b} keys %NNs){
      print OUT $inds[$ind];
      my $minval = 1_000_000; my $minind;
      for (my $i = 0; $i <= $numinds - 1; $i++){
         if($ind != $i){
            if($NNs{$ind}{$i} < $minval){
               $minval = $NNs{$ind}{$i}; $minind = $inds[$i];
            }
         }
      }
      print OUT "\,$minind\,$minval\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Combine(){

   my %CombinedData;
   my $file = "d:/perl64/programs/1KG/PoznikData.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $headers = shift @data;
   my $numinds = @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $CombinedData{$linedata[0]}{"pop"} = $linedata[1];
      $CombinedData{$linedata[0]}{"basal"} = $linedata[2];
      $CombinedData{$linedata[0]}{"macro"} = $linedata[3];
      $CombinedData{$linedata[0]}{"haplo"} = $linedata[4];
      $CombinedData{$linedata[0]}{"coverage"} = $linedata[5];
      $CombinedData{$linedata[0]}{"missing"} = $linedata[6];
      $CombinedData{$linedata[0]}{"PoznikPrivates"} = $linedata[7];
      $CombinedData{$linedata[0]}{"subtree"} = $linedata[8];
      $CombinedData{$linedata[0]}{"branch"} = $linedata[9];
   }
   my $file = "d:/perl64/programs/1KG/YvariantsCP privates.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $numinds = @data;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $CombinedData{$linedata[0]}{"MyPrivates"} = $linedata[1];
   }
   my $file = "d:/perl64/programs/1KG/CombinedPoznikMyData.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "Ind\,coverage\,missing\,PoznikPrivates\,MyPrivates\,Population\,BasalSNP\,MacroHaplogroup\,Haplogroup\,subtree\,branch\n";
   foreach my $ind (sort {$a cmp $b} keys %CombinedData){
      print OUT "$ind";
      print OUT "\,", $CombinedData{$ind}{"coverage"};
      print OUT "\,", $CombinedData{$ind}{"missing"};
      print OUT "\,", $CombinedData{$ind}{"PoznikPrivates"};
      print OUT "\,", $CombinedData{$ind}{"MyPrivates"};
      print OUT "\,", $CombinedData{$ind}{"pop"};
      print OUT "\,", $CombinedData{$ind}{"basal"};
      print OUT "\,", $CombinedData{$ind}{"macro"};
      print OUT "\,", $CombinedData{$ind}{"haplo"};
      print OUT "\,", $CombinedData{$ind}{"subtree"};
      print OUT "\,", $CombinedData{$ind}{"branch"};
      print OUT "\n";
   }
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub TestCGVCF(){

   my $file = "d:/perl64/programs/1KG/CG/vcfBeta-HG00731-200-37-ASM.vcf";
   open (IN, "<$file") || die("Couldn't open <$file");
   $file = "d:/perl64/programs/1KG/CG/HG00731 Y.vcf";
   open (OUT, ">$file") || die("Couldn't open >$file");

   my $counter;
   while(<IN>){
      my $line = $_;
      if(substr($line, 0, 1) eq "Y"){
         print OUT $line;
      }
      $counter++;
   }
   close IN;
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub TestCGTSV(){

   my $file = "d:/perl64/programs/1KG/CG/Complete_Public_Genomes_54genomes_VQHIGH_testvars.tsv";
   open (IN, "<$file") || die("Couldn't open <$file");
   $file = "d:/perl64/programs/1KG/CG/Y.tsv";
   open (OUT, ">$file") || die("Couldn't open >$file");

   my $headers = <IN>;
   print OUT $headers;
   while(<IN>){
      my $line = $_;
      my @linedata = split(/\t/, $line);
      if($linedata[1] eq "chrY"){
         print OUT $line;
      }
   }
   close IN;
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub CGPrivates(){

   print "   Loading data\n";
   my $file = "d:/perl64/programs/1KG/CG Y variants.fas";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN;
   my $data = join(//, @data);
   @data = split(/\>/, $data);
   my $numinds = @data;
   my @sequences; my @inds; my $counter = 0; my $numnucs;
   foreach my $line (@data){
      my @linedata = split(/\n/, $line);
      push @inds, $linedata[0];
      my @seq = split(//, $linedata[1]);
      $numnucs = @seq;
      for (my $loc = 0; $loc <= $numnucs - 1; $loc++){
         $sequences[$counter][$loc] = $seq[$loc];
      }
      $counter++;
   }
   
   print "   Calculating Privates...\n";
   my %privates; my $private;
   for (my $loc = 0; $loc <= $numnucs - 1; $loc++){
      my $count = 0;
      for (my $peep = 0; $peep <= $numinds - 1; $peep++){
         if($sequences[$peep][$loc] == 0){ $private = $peep; $count++;}
         if ($count > 1){last;}
      }
      if($count == 1){$privates{$private}++;}
      $count = 0;
      for (my $peep = 0; $peep <= $counter - 1; $peep++){
         if($sequences[$peep][$loc] == 1){ $private = $peep; $count++;}
         if ($count > 1){last;}
      }
      if($count == 1){$privates{$private}++;}
   }
   my $file = "d:/perl64/programs/1KG/CG Y variants privates.csv";
   open (OUT, ">$file") || die("Couldn't open >$file");
   foreach my $ind (sort {$a cmp $b} keys %privates){
      print OUT "$ind\,$inds[$ind]\,$privates{$ind}\n";
   }
   close OUT;
   print "DONE\n";
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub Expand1KG(){

   my %CGlocs;
   my $file = "d:/perl64/programs/1KG/CG/CG Y.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $temp = shift @data; print OUT "$temp\n";
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $CGlocs{$linedata[0] + 1}{1} = $line;
   }
   my $file = "d:/perl64/programs/1KG/CG/1KG Y.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;

   my $file = "d:/perl64/programs/1KG/CG/CG Y expanded.csv";
   open (OUT1, ">$file") || die("Couldn't open >$file");
   my $file = "d:/perl64/programs/1KG/CG/1KG Y expanded.csv";
   open (OUT2, ">$file") || die("Couldn't open >$file");
   my $temp = shift @data; print OUT1 "$temp\n"; print OUT2 "$temp\n";
   my $temp = shift @data;  my $temp = shift @data;  my $temp = shift @data;  my $temp = shift @data;

   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $CGlocs{$linedata[0]}{2} = $line;
   }
   foreach my $loc (sort {$a<=>$b} keys %CGlocs){
      print OUT1 "$loc\,$CGlocs{$loc}{1}\n";
      print OUT2 "$loc\,$CGlocs{$loc}{2}\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub TestYLocs(){

   print "   Loading Reference Y...\n";
   my $file = "d:/perl64/programs/1KG/Homo_sapiens.GRCh38.Y.fa";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data; shift @data;
   my $Y = join('',@data);

   print "   Loading 1KG loc data...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my $inds  = <IN>; chomp $inds;
   my $pops  = <IN>; chomp $pops;
   my $spops = <IN>; chomp $spops;
   my $shaps = <IN>; chomp $shaps;
   my $haps  = <IN>; chomp $haps;

   my @refs; my @alts; my @locs; my $counter = 0;
   while(<IN>){
      my $line = $_;
      chomp $line;
      my @linedata = split(/\,/, $line);
      my $loc = shift @linedata;
      my $ref = shift @linedata;
      my $alt = shift @linedata;
      $locs[$counter] = $loc;
      $refs[$counter] = $ref;
      $alts[$counter] = $alt;
      $counter++;
   }
   close IN;
   
   print "   Testing\n";
   for (my $i = 0; $i <= $counter - 1; $i++){
      my $nuc = substr($Y, $locs[$i], 1);
      print "      $i $locs[$i] $refs[$i] $alts[$i] $nuc\n";
      <STDIN>;
   }
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub CreateFlier(){

   print "   Creating flier...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP Flier 500.meg";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN;
   my $header = shift @data; my $title = shift @data; shift @data;
   my $sequences = join ('', @data); @data = "";
   my @sequences = split("#", $sequences); $sequences = "";

   my $file = "d:/perl64/programs/1KG/YvariantsCP Flier 1000.meg";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "#mega\n";
   print OUT "Title:YvariantsCP Flier.meg\n";
   print OUT "\n";

   my $sequence = shift @sequences;
   my @sequence = split(/\n/, $sequence);
   my $seqname = shift @sequence; print OUT "#$seqname\n";
   $sequence = join('', @sequence);
   print OUT "$sequence"; print OUT "A" x 500; print OUT "\n";
   
   foreach my $sequence (@sequences){
      my @sequence = split(/\n/, $sequence);
      my $seqname = shift @sequence; print OUT "#$seqname\n";
      $sequence = join('', @sequence);
      print OUT "$sequence"; print OUT "G" x 500; print OUT "\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub RecodeMega(){

   print "   Loading Recode Data...\n";
   my $file = "d:/perl64/programs/1KG/1KG Individual and Population Codes.csv";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; shift @data; chomp @data;
   my %recode;
   foreach my $line (@data){
      my @linedata = split(/\,/, $line);
      $recode{$linedata[0]} = $linedata[5];
   }

   print "   Loading Seqeunce Data...\n";
   my $file = "d:/perl64/programs/1KG/YvariantsCP export.fas";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN;
   my $sequences = join ('', @data); @data = "";
   my @sequences = split(">", $sequences); $sequences = "";

   my $file = "d:/perl64/programs/1KG/YvariantsCP Recoded.meg";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT "#mega\nTitle:YvariantsCP Recoded.meg\n\n";
   
   print "   Recoding...\n";
   foreach my $sequence (@sequences){
      my @sequence = split(/\n/, $sequence);
      my $seqname = shift @sequence;
      my $dash = index($seqname, "-");
      $seqname = substr($seqname, $dash + 1, 7); print OUT "#$recode{$seqname}\n";
      print "$seqname $recode{$seqname}\n";
      $sequence = join("\n", @sequence);
      print OUT "$sequence\n";
   }
   close OUT;
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sub ConsensusAllelesToList(){

   my $file = "d:/perl64/programs/1KG/NO Alleles.txt";
   open (IN, "<$file") || die("Couldn't open <$file");
   my @data = <IN>; close IN; chomp @data;
   my $data = join('', @data);
   my $file = "d:/perl64/programs/1KG/NO Alleles list.txt";
   open (OUT, ">$file") || die("Couldn't open >$file");
   print OUT $data;
   close OUT;

}
