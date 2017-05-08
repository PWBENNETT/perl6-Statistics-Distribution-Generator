use v6.c;

unit module Statistics::Distribution;

constant pi = 3.14159265358979323846264338327950288419716939937510;
constant two_pi = 2 * pi;
constant e = exp 1;

class Statistics::Distribution::Generator {

  has @!alts;

  has multi infix:<±> is export = gaussian.new;
  has multi infix:<+/-> is export = gaussian.new;
  has multi infix:<|> is export = _add_alternative;

  use List::AllUtils qw( reduce );

  method FALLBACK {
    if (@!alts) {
      my $accum = reduce { $a + $b } map { $_.weight // 1 } @!alts;
      my $n = rand($accum);
      my $answer;
      for @!alts -> $alt {
        $n -= ($alt.weight // 1);
        if ($n <= 0) {
          $answer = $alt.FALLBACK;
          last;
        }
      }
      return $answer;
    }
    else {
      die "Something horrible has happened";
    }
  }

  sub gaussian is export {
    my ($μ, $σ) = @_;
    $μ ||= 0;
    $σ ||= 1 / 3;
    return Statistics::Distribution::Generator::gaussian.new(:$μ, :$σ);
  }

  sub uniform is export {
    my ($min, $max) = @_;
    $min //= 0;
    $max //= 1;
    return Statistics::Distribution::Generator::uniform.new(:$min, :$max);
  }

  sub logistic is export {
    return Statistics::Distribution::Generator::logistic.new();
  }

  sub supplied is export {
    my ($iv) = @_;
    my $code;
    if (ref $iv eq 'CODE') {
      $code = $iv;
    }
    else {
      $code = sub { return $iv };
    }
    return Statistics::Distribution::Generator::supplied.new(:$code);
  }

  sub gamma is export { return γ(@_) }

  sub γ is export {
    my ($Rorder, $θ) = map { $_ // 1 } @_;
    my $Zorder = $Rorder.Int;
    return Statistics::Distribution::Generator::γ.new(:$Rorder, :$θ, :$Zorder);
  }
  
  sub exponential is export {
    my ($λ) = map { $_ // 1 } @_;
    return Statistics::Distribution::Generator::exponential(:$λ);
  }

  sub _add_alternative {
    my ($lhs, $rhs, $swapped) = @_;
    ($lhs, $rhs) = ($rhs, $lhs) if $swapped;
    $rhs = supplied($rhs) unless ref($rhs) ~~ /^Statistics::Distribution::Generator/;
    my @alts
      = ref($lhs) eq 'Statistics::Distribution::Generator'
      ?? $lhs
      !! ($lhs)
      ;
    $self = Statistics::Distribution::Generator.new(:@alts);
    push @($self.alts), $rhs;
    return $self;
  }

}

class Statistics::Distribution::Generator::gaussian is Statistics::Distribution::Generator {

  has Numeric $.μ;
  has Numeric $.σ;

  method FALLBACK {
    my $U = rand;
    my $V = rand;
    return $.μ + (sqrt(-2 * log $U) * cos($two_pi * $V) * $.σ * 3);
  }

}

class Statistics::Distribution::Generator::uniform is Statistics::Distribution::Generator {

  has Numeric $.min;
  has Numeric $.max;

  method FALLBACK {
    return ($!max - $!min) * rand() + $!min;
  }

}

class Statistics::Distribution::Generator::logistic is Statistics::Distribution::Generator {

  method FALLBACK {
    return -log((1 / _rand_nonzero()) - 1);
  }

}

class Statistics::Distribution::Generator::supplied is Statistics::Distribution::Generator {

  has $.code;

  method FALLBACK {
    return $!code.();
  }

}

class Statistics::Distribution::Generator::γ is Statistics::Distribution::Generator {

  has Numeric $.Rorder;
  has Int $.Zorder;
  has Numeric $.θ;

  sub _tan { sin($_[0]) / cos($_[0]); }

  sub _rand_nonzero {
    my $rv;
    repeat { 1 } while (!($rv = rand));
    return $rv;
  }

  sub _γ_int {
    my $Rorder = shift;
    if ($Rorder < 12){
      my $prod = 1;
      loop (my $i = 0; $i < $Rorder; $i++) {
        $prod *= _rand_nonzero();
      }
      return -log($prod);
    }
    else {
      return _γ_large_int($Rorder);
    }
  }

  sub _γ_large_int {
    my $Rorder = shift;
    my $sqrt = sqrt(2 * $Rorder - 1);
    my ($x,$y,$v);
    repeat {
      repeat {
        $y = _tan(pi * rand);
        $x = $sqrt * $y + $Rorder - 1;
      } while ($x <= 0);
      $v = rand;
    } while ($v > (1 + $y * $y) * exp(($Rorder - 1) * log($x / ($Rorder - 1)) - $sqrt * $y));
    return $x;
  }

  sub _γ_frac {
    my $Rorder = shift;
    my $p = e / ($Rorder + e);
    my ($q, $x, $u, $v);
    do {
      $u = rand;
      $v = _rand_nonzero();
      if ($u < $p){
        $x = exp((1 / $Rorder) * log($v));
        $q = exp(-$x);
      }
      else {
        $x = 1 - log($v);
        $q = exp(($Rorder - 1) * log($x));
      }
    } while (rand >= $q);
    return $x;
  }

  method FALLBACK {
    my $rv;
    if ($!Rorder == $!Zorder) {
      $rv = $!θ * _γ_int($!Zorder);
    }
    elsif ($!Zorder == 0) {
      $rv = $!θ * _γ_frac($!Rorder);
    }
    else {
      $rv = $!θ * (_γ_int($!Zorder) + _γ_frac($!Zorder - $!Rorder));
    }
    return $rv;
  }

}

class Statistics::Distribution::Generator::exponential is Statistics::Distribution::Generator {

  has Numeric $.λ;

  method FALLBACK {
    my $rv = -log(rand) / $!λ;
  }

}
