#################################    PARTITION FUNCTION MODULE

# CONTENTS:
# This Perl module contains the exact but slow variant of the computational methods described in our article:
# E. Tøstesen, F. Liu, T.-K. Jenssen, and E. Hovig. Speed-Up of DNA Melting Algorithm with Complete Nearest Neighbor Properties. Biopolymers, vol. 70, 364-376 (2003).
# The do_partition_function_recursions subroutine runs in time O(N^2) using the exact loop entropy factor

############################################## CONSTANTS

#rescaling parameters:
$norm_factor=1e-30;
$bignumber=1/$norm_factor**2;    #should be some fraction of the machine limit

##############################################   LOOP ENTROPY FACTOR

sub initialize_loop_entropy_factor {
  #calculate omega array (the exact loop entropy factor)
  foreach $k (2..$N) {
    $omega[2*$k]=$sigma*(2*$k+$d)**(-$alfa);
  }
}


############################################## SLOW CALCULATION OF LR AND RL PARTITION FUNCTIONS

sub do_partition_function_recursions {


  ##############  CALCULATE $T AND SALT DEPENDENT QUANTITIES

  &calc_statistical_weights(); #defined in the thermodynamics module


  ##################################  Left-Right (FORWARD) RECURSION INITIALIZATIONS

  $V_LR[1]=1;
  $V_LR[2]=$beta*$s_solo_LR[1];
  $U1_LR[1]=1; #to be used only in &p_helix()
  $U1_LR[2]=$beta;
  $U2_LR[2]=$beta*$s_end_LR[1]*$s_NN_LR[2];
  $V_LR[3]=$U2_LR[2]*$s_end_LR[2]+$U1_LR[2]*$s_solo_LR[2];
  $Z_total_LR=$V_LR[1]+$V_LR[2]+$V_LR[3];
  $level=0;
  $level[1]=0;
  $level[2]=0;
  %normalization_times_RL=();
  %normalization_times_LR=();

  ######################################   Left-Right (FORWARD) ITERATION

  foreach $i (3..$N) {
    $W=0;
    foreach $j (2..$i-1) {                #calculate W_i-1
      $W+=$V_LR[$j]*$omega[2*($i+1-$j)];
      if ($normalization_times_LR{$j}) {
        $W*=$norm_factor;
      }
    }
    $U1_LR[$i]=$beta*$V_LR[1]+$W;
    $U2_LR[$i]=$s_NN_LR[$i]*($U1_LR[$i-1]*$s_end_LR[$i-1]+$U2_LR[$i-1]);
    $V_LR[$i+1]=$U2_LR[$i]*$s_end_LR[$i]+$U1_LR[$i]*$s_solo_LR[$i];
    $Z_total_LR+=$V_LR[$i+1];

    if ($Z_total_LR>$bignumber) {  #normalize when above threshold
      $Z_total_LR*=$norm_factor;
      $U1_LR[$i]*=$norm_factor;
      $U2_LR[$i]*=$norm_factor;
      $V_LR[$i+1]*=$norm_factor;
      $V_LR[1]*=$norm_factor;
      $normalization_times_LR{$i}++;      #use for W sum (forward)
      $normalization_times_RL{$N+2-$i}++; #remember for RL (backward)
      $level++;
    }
    $level[$i]=$level;  #to be used in loop probs calculation
  }


  ##################################  Right-Left (BACKWARD) RECURSION INITIALIZATIONS

  $V_RL[1]=1;
  $V_RL[2]=$beta*$s_solo_RL[1];
  $U1_RL[1]=1; #to be used only in &p_helix()
  $U1_RL[2]=$beta;
  $U2_RL[2]=$beta*$s_end_RL[1]*$s_NN_RL[2];
  $V_RL[3]=$U2_RL[2]*$s_end_RL[2]+$U1_RL[2]*$s_solo_RL[2];
  $Z_total_RL=$V_RL[1]+$V_RL[2]+$V_RL[3];

  if ($normalization_times_RL{2}) { #normalize if it's time
    $Z_total_RL*=$norm_factor;
    $U1_RL[2]*=$norm_factor;
    $U2_RL[2]*=$norm_factor;
    $V_RL[3]*=$norm_factor;
    $V_RL[1]*=$norm_factor;
  }

  ######################################   Right-Left (BACKWARD) ITERATION

  foreach $i (3..$N) {
    $W=0;
    foreach $j (2..$i-1) {                #calculate W_i-1
      $W+=$V_RL[$j]*$omega[2*($i+1-$j)];
      if ($normalization_times_RL{$j}) {
        $W*=$norm_factor;
      }
    }
    $U1_RL[$i]=$beta*$V_RL[1]+$W;
    $U2_RL[$i]=$s_NN_RL[$i]*($U1_RL[$i-1]*$s_end_RL[$i-1]+$U2_RL[$i-1]);
    $V_RL[$i+1]=$U2_RL[$i]*$s_end_RL[$i]+$U1_RL[$i]*$s_solo_RL[$i];
    $Z_total_RL+=$V_RL[$i+1];

    if ($normalization_times_RL{$i}) { #normalize when it's time
      $Z_total_RL*=$norm_factor;
      $U1_RL[$i]*=$norm_factor;
      $U2_RL[$i]*=$norm_factor;
      $V_RL[$i+1]*=$norm_factor;
      $V_RL[1]*=$norm_factor;
    }
  }


  ######################################### END OF SLOW LR AND RL

} #End of do_partition_function_recursions


############################################ HOW TO CALCULATE VARIOUS PROBABILTIES
sub p_closed {
  my $i=$_[0];
  if($i==0){
    $p_closed=0;
  }
  elsif ($i==1) {
    $p_closed=$V_RL[$N+1]/$Z_total_RL;
  } 
  elsif ($i==$N) {
    $p_closed=$V_LR[$N+1]/$Z_total_LR;
  } 
  else {
    $p_closed=( $U1_LR[$i]*$s_solo_LR[$i]*$U1_RL[$N+1-$i]
               +$U1_LR[$i]* $s_end_LR[$i]*$U2_RL[$N+1-$i]
	       +$U2_LR[$i]* $s_end_LR[$i]*$U1_RL[$N+1-$i]
	       +$U2_LR[$i]               *$U2_RL[$N+1-$i])/($beta*$Z_total_LR);
  }
}

sub p_right {
  my $x=$_[0];
  $p_right=($V_LR[$x+1]/$Z_total_LR)*$norm_factor**($level-$level[$x]);
}

sub p_left {
  my $y=$_[0];
  $p_left=($V_RL[$N+2-$y]/$Z_total_RL)*$norm_factor**$level[$y];
}

sub p_loop {
  my ($x,$y)=@_;
    $Z_xy=$V_LR[$x+1]*$omega[2*($y-$x)]*$V_RL[$N+2-$y]/$beta;
    $p_loop=$Z_xy/$Z_total_LR;
    $p_loop*=$norm_factor**($level[$y]-$level[$x]);
}

sub multiply_helix_chain {
  #calculates @xi_LR, $xi_total_LR, %xi_normalization_LR (also for RL) and @xi_level
  my $level=0;
  %xi_normalization_RL=();
  %xi_normalization_LR=();
  # FORWARD:
  $xi_LR[1]=1;
  $xi_total_LR=1;
  $xi_level[1]=0;
  foreach $i (2..$N) {
    $xi_total_LR*=$s_NN_LR[$i];
    if ($xi_total_LR>$bignumber) {
      $xi_total_LR*=$norm_factor;
      $xi_normalization_LR{$i}++;
      $xi_normalization_RL{$N+2-$i}++;
      $level++;
    }
    if ($xi_total_LR<(1/$bignumber)) {
      $xi_total_LR/=$norm_factor;
      $xi_normalization_LR{$i}--;
      $xi_normalization_RL{$N+2-$i}--;
      $level--;
    }
    $xi_LR[$i]=$xi_total_LR;
    $xi_level[$i]=$level;
  }
  # BACKWARD:
  $xi_RL[1]=1;
  $xi_total_RL=1;
  foreach $i (2..$N) {
    $xi_total_RL*=$s_NN_RL[$i];
    if ($xi_normalization_RL{$i}) {
      $xi_total_RL*=$norm_factor**$xi_normalization_RL{$i};
    }
    $xi_RL[$i]=$xi_total_RL;
  }
}
sub p_helix {
  my ($x,$y)=@_;
  &multiply_helix_chain() if (!@xi_LR); #if xi is empty (first time)
  if ($x==$y) { #an isolated bp
    $Z_xy=$U1_LR[$x]*$s_solo_LR[$x]*$U1_RL[$N+1-$y]/$beta;
    $p_helix=$Z_xy/$Z_total_LR;
  } elsif ($x<$y) { #helix with at least two bp
    $Z_xy=$U1_LR[$x]*$s_end_LR[$x]*$s_end_LR[$y]*$U1_RL[$N+1-$y]/$beta;
    $Z_xy*=$xi_total_LR/($xi_LR[$x]*$xi_RL[$N+1-$y]); #the product from x+1 to y
    $p_helix=$Z_xy/$Z_total_LR;
    $p_helix*=$norm_factor**($xi_level[$x]-$xi_level[$y]+$level[$y]-$level[$x]);
  }
  return $p_helix;
}
###########  ###########  THE END...! Eivind Tostesen, 2004 :-)  ###########  ###########


