package DNAMelting::Meltingprofiles;

require Exporter;
@ISA = (Exporter);

@EXPORT = qw(calc_melting_temp_interval calc_seq_prob make_fixed_temp_list make_temp_list init_require_pack calc_tm_prof);

#use LibSN::DNAMelting::Gnuplot qw(make_xy_batchfile plot_a_graph plot_scaled_tempprofile random_file_name);

#use FindBin;
#use lib "$FindBin::Bin";
#use FileHandle;
# autoflush;
# version: Eivind Tostesen. January 29, 2004.
# Use with: perl meltingprofiles.pl bacteriophagelambda.txt

# package variables:
my $_results = "results"; # the catalog where the results go when done

# class member variables (cross global package variables.. fyfy):
$m_sequence = "";
$N = 0;

sub calc_Tm { 
  my ($pm,$tol)=@_; #for example (0.5,0.001)
  #*print "($pm,$tol)";
  my $Fm=log($pm/(1-$pm)); 

  #Important: Tm must be BETWEEN the brackets T1 and T2 !! 

  #left bracket: 
  my $T1=60+273.15; 
  my $V1=1/$T1; 
  $T=$T1; 
  &do_partition_function_recursions(); 
  my $p1=0; 
  foreach $i (1..$N) { 
    $p1+=&p_closed($i); #defined in the partition function module 
  } 
  $p1/=$N;
  my $F1=log($p1/(1-$p1)); 

  #right bracket: 
  my $T2=100+273.15; 
  my $V2=1/$T2; 
  $T=$T2; 
  &do_partition_function_recursions(); 
  my $p2=0; 
  foreach $i (1..$N) { 
    $p2+=&p_closed($i); #defined in the partition function module 
  } 
  $p2/=$N; 
  my $F2=log($p2/(1-$p2)); 

  my $pnew; 
  my $Tm=0; 
  my $Vm; 
  #my $dT=1000; 
  my $dp=1000; 
  while (abs($dp)>$tol) { 
    $Vm=$V1*(($Fm-$F2)/($F1-$F2))+$V2*(($Fm-$F1)/($F2-$F1)); 
    #$dT=1/$Vm-$Tm; 
    $Tm=1/$Vm; 
    $T=$Tm; 
    &do_partition_function_recursions(); 
    $pnew=0; 
    foreach $i (1..$N) { 
      $pnew+=&p_closed($i); #defined in the partition function module 
    } 
    $pnew/=$N; 
    $dp=$pnew-$pm; 
    if ($p1>$pnew && $pnew>$pm) { 
      $F1=log($pnew/(1-$pnew)); 
      $V1=1/($Tm); 
    } 
    if ($pm>$pnew && $pnew>$p2) { 
      $F2=log($pnew/(1-$pnew)); 
      $V2=1/($Tm); 
    } 
  } 
  return $Tm; 
} #end of calc_Tm 



sub calc_average_p_closed {
  foreach $i (1..$N) {
    $p_closed_average+=&p_closed($i); #defined in the partition function module
  }
  $p_closed_average/=$N;
}

sub make_pclosed_file {
  my $pclosedfilename=$_[0];
  # probability profile in gnuplot format:
  open POUT, ">$pclosedfilename"
    or die "cannot create file: $!";
  foreach $i (1..$N) {
    my $p=&p_closed($i); #defined in the partition function module
    print POUT "$i $p\n";
  }
  close POUT;
  #/ every new file needs proper acccess
  chmod(0755, $pclosedfilename);
}


# Melting curve - create a graph showing total possibility for each temperature
sub calc_melting_temp_interval{
  #--- arguments:
  my ($temp_min, $temp_max, $step_type, $step_value, $algorithm, $thermodyn, $seq, $nat_con) = @_;

  # dynamic initiation
  init_require_pack($algorithm, $thermodyn, $seq, $nat_con);

  # containing 2D-plots for the graph
  my %xyList = ();

  # abort if users input is invalid
  my  $abort = "false";

  # taking the time
  my $before = time();
  #*print "<p> temp: $temp_min $temp_max</p>\n";
  # find the maximum number of steps
  my $max_time = 120; #sec
  my $max_thoughput = 4000; #bp/sec
  my $max_num_of_steps = 40;
  my $maxsteps = ($max_time*$max_thoughput)/(length($seq));
  
  # minimum number of points/steps/positions in the plot
  my $minsteps = 2;
  
  # *** Relative stepsize
  if ($step_type eq "relative" && $maxsteps >= $minsteps){ # zero if the user has specified relative resolution && check if the sequence is too long
    # if the number of steps available ex. $maxsteps = 120, then the plot is going to be 'overplotted'.
    if ($maxsteps > $max_num_of_steps){
      $maxsteps = $max_num_of_steps;
    }
    # Find temperatures in the interval based on the resolution of the given stepsize
    if($step_value eq "Small" && ($maxsteps*0.5 >= $minsteps)){
      $maxsteps = $maxsteps*0.5; # 50% of max
    }
    elsif($step_value eq "Medium" && ($maxsteps*0.75 >= $minsteps)){
      $maxsteps = $maxsteps*0.75; # 75% of max
    }
    elsif($step_value eq "High"){
      # Max
    }
    else{
      print "<p> ERROR: The step resolution was too low. Probably because the DNA sequence was too long. Max one step is: 480000bp and select \"Relative Step: High\" to get the desired value. </p>\n";
      print "<p>Number of steps was:  $maxsteps</p>\n";
      $abort = "true";
    }
    if ($abort eq "false"){
      # contains the interval of temperatures
      @temp_interval = ();
      my $stepsize = (($temp_max - $temp_min)/$maxsteps);
      #*print "SS: $stepsize\n";
      #*print "MS: $maxsteps\n";

      # create an interval that overstep its border ex. [1-100] with stepsize 1.2 goes to 100.8
      $maxsteps = $maxsteps + 1;
      # make a list of temperatures between the interval [min_temp, max_temp] with a stepsize that determine the resolution of the plot
      for ($i = 0; $i < $maxsteps; $i ++) {
	my $new_temp = $temp_min + ($i*$stepsize);
	push(@temp_interval, $new_temp );
      }
      # Calculate the probability plot for each temperature:
      foreach $temp(@temp_interval){
	$T = $temp;
	$TC=$T-273.15;
	&do_partition_function_recursions(); #defined in the partition function module
	&calc_average_p_closed();
	$xyList{$temp} = $p_closed_average;
	#*print "<p> $p_closed_average is closed at $TC C</p>\n";
      } 
    }
  }
  # *** Manual stepsize
  elsif ($step_type eq "manually" && $maxsteps >= 1){
    $stepsize = $step_value; # to be used in the error message
    # number of steps in the given interval. Counts from 0.
    my $num_steps = (($temp_max - $temp_min)/$step_value);
    # $temp_max - $temp_min => interval [1..$maxsteps]. Check if the interval was too short or if the interval was too high
    if ($num_steps >= $minsteps && $num_steps <= $maxsteps){
      # create an interval that overstep its border ex. [1-100] with stepsize 1.2 goes to 100.8 only if stepsize does not fill the interval.
      if(($num_steps*$step_value)+$temp_min > $temp_max){ # ex. (intv[80..90]C; step_value=1C;) => num_steps=10; (10x*1C)+ 80C = 90C
	$temp_max = $temp_max + $step_value;	
      }
      my $i = 0;
      # iterate through the interval with the given step-value:
      for ($i = $temp_min; $i <= $temp_max; $i = $i + $step_value) {
	$T = $i;
	$TC=$T-273.15;
	&do_partition_function_recursions(); #defined in the partition function module
	&calc_average_p_closed();
	$xyList{$i} = $p_closed_average;
	#*print "<p> $p_closed_average is closed at $TC C</p>\n";
      }
    }
    else{
      $min_stepsize = (($temp_max - $temp_min)/$maxsteps);
      $user_interval = ($temp_max - $temp_min);
      $user_steps = (($temp_max - $temp_min)/$stepsize);
      print "<p>ERROR: No plots to make! Interval: ($temp_max - $temp_min)=>$user_interval) and [min/max] stepsize was: [$min_stepsize/$user_interval], user had: $stepsize</p>\n";
      print "<p> Maximum number of steps : $maxsteps, and the user had: $user_steps steps <br/> User has either specified:
<br/>1. Too short interval with high resolution. <br/>2.Too long interval with high interval. <br/>3. Too short interval with low resolution. <br/>4. Too long interval with a low resolution.</p>\n";
      $abort = "true";
    }
  }
  # running time
  my $run_time =(time()-$before);
  #*print "<p>Runtime - Calc_melting_temp_interval: $run_time </p>\n";

  if ($abort eq "false"){
    # make a graph
    my @xyfile = (Gnuplot::make_xy_batchfile(\%xyList));
    # the function returns one value
    my $webplotfile = $xyfile[0];
    # converting kelvin to celsius and sending the coordinates list to the Gnuplot
    my $graphfilename =  Gnuplot::plot_a_graph($TC, ($temp_min - 273.15), ($temp_max - 273.15), \@xyfile, "prob_prof");
    return ($graphfilename, $webplotfile);
  }
  else {
    return ("abort");
  }
}

sub calc_seq_prob{
  #*print "<p>PP: calc_seq_prob </p>\n";
  #--- arguments:
  ($tempList, $algorithm, $thermodyn, $seq, $nat_con) = @_;
  
  # ET and VN 2010-12-08: modified check for "Exact" and "Approx"
  
  my $max_time = 120; #sec (This is the Apache server timeout)
  my $estcomptime; # estimated comp time per temperature
  
  if ($algorithm eq "Exact") {
    $estcomptime=9.52e-7*length($seq)**2;
  } elsif ($algorithm eq "Approx") {
    $estcomptime=6.66e-5*length($seq);
  }
  # check total comp time is below maximum.
  if ($max_time > @{$tempList}*$estcomptime) {
    # taking the time
    my $before = time();
    init_require_pack($algorithm, $thermodyn, $seq, $nat_con);
    
    # a list of all the files to be plotted
    @plotList = ();

    foreach $temp(@{$tempList}){
      $T = $temp;
      &do_partition_function_recursions(); #defined in the partition function module
      $TC=$T-273.15;
      # unique filenaming
      my $randfilename = Gnuplot::random_file_name(1,1200000000);
      my $closed_file = "$_results/pclosed_profileT$TC\_$randfilename";
      &make_pclosed_file("$closed_file");
      # put the filenames in a list
      push(@plotList, "$closed_file");
    }
    # running time
    my $run_time =(time()-$before);
    #*print "<p>Runtime for \"calc_seq_prob\": $run_time </p>\n";

    # plot the files into one graph
    my $graphfilename = Gnuplot::plot_a_graph("calc_seq_prob", 1, $N, \@plotList, "prob_prof");
    return $graphfilename;
  }
  else{
    # ET and VN 2010-12-08: changed error message
    print "<p>Error! Too long sequence or too many temperatures or both.</p>\n";
    return ("abort");
    
  }
}

# make a list of temperatures based on a fixed list of helicities
sub make_fixed_temp_list{
  # sending the arguments into the initiator
  init_require_pack(@_);
  my @temp_list = ();
  $tmp_temp = &calc_Tm(0.10,0.0001);
  push (@temp_list, $tmp_temp);
  $tmp_temp = &calc_Tm(0.25,0.0001);
  push (@temp_list, $tmp_temp);
  $tmp_temp = &calc_Tm(0.50,0.0001);
  push (@temp_list, $tmp_temp);
  $tmp_temp = &calc_Tm(0.75,0.0001);
  push (@temp_list, $tmp_temp);
  $tmp_temp = &calc_Tm(0.90,0.0001);
  push (@temp_list, $tmp_temp);
 return @temp_list;
}

# make a list of temperatures based on a list of helicities (almost the same as make_fixed_temp_list)
# this function is linear in time with a processing at 0.625kpb pr sec on the pubgeneserver
sub make_temp_list{
  ($hel_list, $algorithm, $thermodyn, $seq, $nat_con) = @_;
  #*print "Args: make_temp_list = ($hel_list, $algorithm, $thermodyn, $nat_con and a sequence input)\n";
  # find the maximum number of calculations
  my $max_time = 120; #sec
  my $max_thoughput = 600; #bp/sec
  my $max_calculations = ($max_time*$max_thoughput)/(length($seq));
  # check number of calulations is below maxium.
  if($max_calculations > @{$hel_list}){
    my $before = time();
    # sending the arguments into the initiator
    init_require_pack($algorithm, $thermodyn, $seq, $nat_con);
    my @temp_list = ();
    foreach $hel(@{$hel_list}){
      #*print "$hel\n";
      $tmp_temp = &calc_Tm($hel,0.001);
      push (@temp_list, $tmp_temp);
    }
    # running time
    my $run_time =(time()-$before);
    #*print "<p>Runtime for \"make_temp_list\": $run_time </p>\n";
    return @temp_list;
  }
  else{
    my $len_user_seq = length($seq);
    my $max_seq_len = ($max_time*$max_thoughput);
    if($max_seq_len <  $len_user_seq){
      print "<p>Error! Too long sequence. User had: $len_user_seq. Maximum sequence length: $max_seq_len</p>\n";
    }
    else{
      my $cuml_time = ((@{$hel_list}*$len_user_seq)/$max_thoughput);
      print "<p>Error! Too many sequences. Cumulated amount of calculation was: $cuml_time. Max time was: $max_time</p>\n";
    }
    return("abort");
  }
}

sub calc_tm_prof{
  my ($helicity, $algorithm, $thermodyn, $seq, $nat_con) = @_;
  my $before = time();
  # sending the arguments into the initiator
  init_require_pack($algorithm, $thermodyn, $seq, $nat_con);
  #*print "Tm_prof: ($helicity, $algorithm, $thermodyn, $seq, $nat_con) \n";
  my @TL = &calc_temp_profile($helicity);

  # unique filenaming
  ##my $randomname = Gnuplot::random_file_name(1,1200000000);
  ##my $tempprof_filename = "$_results/tempprof_profile_$randomname.txt";
  ##open POUT, ">$tempprof_filename"
    ##or die "cannot create file: $!";
  ##foreach $i (1..$N) {
    ##print POUT "$i $TL[$i]\n";
  ##}
  ##close POUT;
  ### make the right rights to the right user
  ##chmod(0755, "$tempprof_filename");

  ##my $after =(time()-$before);
  ###*print "<p>Runtime - calc_tm_prof: $after sec\n</p>";
  ##my $pngfile = "$_results/T-prof1_$randomname.png";
  ### create the plotfile
  ##Gnuplot::plot_scaled_tempprofile("$tempprof_filename","$pngfile",1,$N);
  ##return ($pngfile, $tempprof_filename);

  return @TL;
}


sub init_require_pack{
  ($algorithm, $thermodyn, $seq, $nat_con) = @_; # ET: my?
  # print "\ninit_argv: @_\n";

  #--- Read in the sequence:
  $m_sequence = $seq;
  #an input dependent constant:
  $N=length($m_sequence);

  #---Choose one thermodynamics module:
  #if ($thermodyn  eq  "Blake and Delcourt 1998") {
    #require BlakeDelcourt98;}
  #if ($thermodyn  eq  "Fixman and Freire 1977") {
    #require FixmanFreire77;}
  #if ($thermodyn  eq  "SantaLucia (polymer) 1998") {
    #require SantaLucia98;}
  #if ($thermodyn  eq  "SantaLucia (oligo) 1998") {
    #require SantaLucia_oligo;}
  if ($thermodyn  eq  "Blossey and Carlon 2003") {
    require DNAMelting::BlosseyCarlon03;}
  #if ($thermodyn  eq  "Gotoh and Tagashira 1981") {
    #require GotohTagashira81; }

  #---CHOOSE ONE PARTITION FUNCTION MODULE:
  if ($algorithm eq "Approx") {
    require DNAMelting::Yeramian_Tostesen_approx;
  }
  elsif ($algorithm eq "Exact") {
    require Yeramian_Tostesen_exact;
  }
  &initialize_loop_entropy_factor();
  $Na_conc=$nat_con;
  #$T=86+273.15; #the temperature = 359,15
  $T = $temp_min; #G1
}

# HOW TO CALCULATE A SINGLE T-PROFILE
sub calc_temp_profile {
  my ($pL)=@_; #for example, 0.75, 0.5, or 0.25
  # print "calc_temp_prof: $pL\n";
  #Initialize
  $T=59+273;
  @TL = (); @p1 = ();
  &do_partition_function_recursions(); #defined in the partition function module
  foreach $i (1..$N) {
    $p1[$i]=&p_closed($i); #defined in the partition function module
  }
  #print "<p>p1# @p1 N: $N\n</p>";
  #print "<p>pL: $pL\n</p>";
  #temperature scan
  for $T (40+273..120+273) {
		&do_partition_function_recursions(); #defined in the partition function module
		foreach $i(1..$N){
		#foreach $i (1..$N) {
			 $p2=&p_closed($i); #defined in the partition function module
			 if (($p1[$i]-$pL)*($p2-$pL)<=0) {   #if pL between p1 and p2 then find TL
				  #Now do sigmoidal interpolation
				  if (abs($p1[$i]-0.5)<0.5 && abs($p2-0.5)<0.5) { #if not 0 or 1 then
						# from sigmoidal to linear quantities
						$FL=log($pL/(1-$pL));
						$F1=log($p1[$i]/(1-$p1[$i]));
						$F2=log($p2/(1-$p2));
						$V1=1/($T-1);
						$V2=1/($T);
						#weighted interpolation
						$VL=$V1*(($FL-$F2)/($F1-$F2))+$V2*(($FL-$F1)/($F2-$F1));
						$TL_2 = 1/$VL-273;
				  } else { #else (0 or 1) do linear interpolation
						$TL_2 = -273+($T-1)*(($pL-$p2)/($p1[$i]-$p2))+$T*(($pL-$p1[$i])/($p2-$p1[$i]));
				  }
				  $TL[$i]=$TL_2;
				  # print "<p>TL_2# $TL_2  N: $N\n</p>";
			 } # if
			 $p1[$i]=$p2;
		}
  }
  return (@TL);
} #end of calc_temp_profile


1;
