##################  THERMODYNAMICS MODULE

#This is the set of R. Blossey and E. Carlon (2003): Reparametrizing loop entropy weights: Effect on DNA melting curves. Phys. Rev. E, 68, 061911.
#It's a reparametrization of the parameters of R. D. Blake and S. G. Delcourt (1998): Thermal stability of DNA. Nucl. Acids Res. 26, 3323-3332.
#such that alfa=2.15 and sigma=1.26e-4

#The usage is: open-to-closed: sij=exp((-dHij/RTij)*(1-Tij/T)),
#where Tij=T0ij+log10[Na+]*dTij/dlog10[Na+]. T0ij is for 1.0 M salt.
#Below are tables of dHij, T0ij and dTij/dlog10[Na+] (in kcal/mol and degrees C).

my $beta = 1;

%T0ij=(
  "AA=TT" => 89.08,
  "AT" => 89.38,
  "TA" => 79.47,
  "CA=TG" => 89.71,
  "AC=GT" => 121.17,
  "AG=CT" => 98.49,
  "GA=TC" => 105.09,
  "CG" => 105.28,
  "GC" => 143.73,
  "CC=GG" => 118.49,
);

%dHij=(
  "AC=GT" => 10.51,
  "TA" => 7.81,
  "AA=TT" => 8.45,
  "GC" => 11.91,
  "CC=GG" => 10.34,
  "GA=TC" => 9.47,
  "AT" => 8.50,
  "CG" => 9.53,
  "AG=CT" => 9.10,
  "CA=TG" => 8.51,
);

%dTij_dlog10Na=(
  "AC=GT" => 13.71,
  "TA" => 22.08,
  "AA=TT" => 19.78,
  "GC" => 10.21,
  "CC=GG" => 14.18,
  "GA=TC" => 16.94,
  "AT" => 19.19,
  "CG" => 16.01,
  "AG=CT" => 17.30,
  "CA=TG" => 19.35,
);

%nn2ten=(
  AA => "AA=TT",
  TT => "AA=TT",
  AT => "AT",
  TA => "TA",
  CA => "CA=TG",
  TG => "CA=TG",
  "GT" => "AC=GT",
  AC => "AC=GT",
  CT => "AG=CT",
  AG => "AG=CT",
  GA => "GA=TC",
  TC => "GA=TC",
  CG => "CG",
  GC => "GC",
  GG => "CC=GG",
  CC => "CC=GG",
);

sub NN {
# &NN(i,seq) returns the two nearest neighbor bases (i-1,i) in the sequence ending at position i (i>1 i.e. it starts at second base)
  substr($_[1],$_[0]-2,2);
}

$R=1.987;  #the gas constant
$d=1;
$sigma=1.26e-4;
$alfa=2.15; #to be used with exact loop entropy factor

#Yeramian_Tostesen_approx.pm uses the following variable:
$expfname="ET-2.15-14.muex"; #file should correspond to value of $alfa

sub calc_statistical_weights {
# Calculation of salt- and temperature dependent statistical weights

  $beta=1; #should depend on T and Na+ for short sequences

  # define solo- and end-statistical weights (trivial constants):
    foreach $i (1..$N) {
      $s_solo_LR[$i]=0;
      $s_end_LR[$i]=1;
      $s_solo_RL[$N+1-$i]=0;
      $s_end_RL[$N+1-$i]=1;
    }

  #calculate Tij at this salt concentration:
  while (($key,$value)=each %T0ij) {
    $Tij{$key}=$value+273.15+(log($Na_conc)/log(10))*$dTij_dlog10Na{$key};
  }

  # how to calculate LR and RL NN-statistical weights (BD98)
  foreach $i (2..$N) {
    $ij=$nn2ten{&NN($i,$m_sequence)};
    $s_NN_LR[$i]=exp((-1000*$dHij{$ij}/($R*$Tij{$ij}))*(1-$Tij{$ij}/$T));
    $s_NN_RL[$N+2-$i]=$s_NN_LR[$i];
  }
}
