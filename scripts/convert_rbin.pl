#!/usr/bin/perl -w
use Cwd;
use File::Path;
use File::Basename;

my $dirname = dirname(__FILE__);
my $rbin = "$dirname/../examples/rbin_to_ascii_particles";
unless (-x $rbin) {
    print "Couldn't find rbin converter in expected location ($rbin).\n";
    print "Make sure that you ran \"make\" in the examples/ subdirectory";
}

unless (@ARGV>0) {
    die "Usage: $0 /path/to/rockstar.cfg\n";
}

$c = new ReadConfig($ARGV[0]);
my $ns = $c->{NUM_SNAPS};
my $nb = $c->{NUM_WRITERS};

my @snapnames = (0..($ns-1));
if ($c->{SNAPSHOT_NAMES}) {
    open IN, "<", $c->{SNAPSHOT_NAMES} or die "Could not open $c->{SNAPSHOT_NAMES}!\n";
    @snapnames = <IN>;
    chomp(@snapnames);
    close IN;
}

close STDOUT;
for my $x (0..$#snapnames) {
    my @fns = map { "$c->{OUTBASE}/halos_$snapnames[$x].$_.rbin_wl" } (0..($nb-1));
    my $out_fn = "$c->{OUTBASE}/wl_particles_$x.txt";
    open STDOUT, ">", "$out_fn" or die ("Couldn't open $out_fn for writing!\n");
    print STDERR "Creating $out_fn...";
    system($rbin, @fns);
    close STDOUT;
    print STDERR "done\n";
}


package ReadConfig;

sub new {
    my ($class, $file) = @_;
    return bless {}, $class unless defined $file;
    open FILE, "<", $file or 
	die "Couldn't open file $file for reading!\n";
    local $_;
    my %config;
    while (<FILE>) {
	s/\#.*//;
	my ($key, $value) = /([^=]+)=(.*)/;
	next unless defined $value;
	#print "$key = $value\n";
	($key, $value) = map { trim($_) } ($key, $value);
	next unless length($key) and length($value);
	$config{$key} = $value;
    }
    close FILE;
    bless \%config, $class;
}

sub trim {
    my $val = shift;
    $val =~ s/^\s*['"]?//;
    $val =~ s/['"]?\s*$//;
    return $val;
}

sub set_defaults {
    my $self = shift;
    my %opts = @_;
    %$self = (%opts, %$self);
}

sub print_config {
    my ($self, $fh) = @_;
    $fh = \*STDOUT unless $fh;
    for (sort keys %$self) {
	print $fh "\"$_\" = \"$self->{$_}\"\n";
    }
}
