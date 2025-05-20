#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Sort::Naturally;
use File::Path qw(remove_tree);
use File::Copy "cp";
$|++;

=head1  NAME

NGS_MauveReorder_on_ref.pl

=head1 DESCRIPTION

Use Mauve to align draft assembly(ies) (fasta) on ref genome (gbk or fasta),
and reorder contigs,
eventually using GenBank file to transfer annotations.

=head1  USAGE

=over

=item -fa_ref <s>

Fasta file of the reference genome.

OR

=item -gb_ref <s>

Genbank file for annotation (better, will provide orthologs).

=item -draft <s>n

Fasta file(s) of the draft genome(s) (treated independantly).

=item [-p <i>]

Default:14, number of threads to use.

=item [-res_dir <s>]

To provide results directory (must already exist).

=item [-test]

Run a test of the program.

=back

example:
./NGS_MauveReorder_on_ref.pl -gb_ref genome.gbk -draft assembly_mira.fasta

=cut

my $prog_tag               = '[NGS_MauveReorder_on_ref.pl]';
my $b_test                 = 0; # ok 2017 11 27
my @res                    = `which Mauve`;
chomp($res[0]);
my $mauve_jar              = $res[0];
chomp($mauve_jar);
# $mauve_jar .= '.jar';
$mauve_jar .= 'CM';
# -e $mauve_jar or die "$prog_tag [Error] Cannot locate Mauve.jar:$mauve_jar, line ".__LINE__."\n";
-e $mauve_jar or die "$prog_tag [Error] Cannot locate MauveCM:$mauve_jar, line ".__LINE__."\n";

my $fa_ref                 = undef;
my $gb_ref                 = undef;
my $nb_threads             = 14;
my @drafts                 = ();

my $sample_name            = 'sample';
my $kmer_size              = undef;
my $res_dir                = undef;

# **********************************************************************
# CHECK OPTIONS
# **********************************************************************
my $nbargmini = 1;
if(scalar(@ARGV) < $nbargmini)
{ 
  print "Bad number of arguments: ".scalar(@ARGV)." found, at least $nbargmini wanted\n";
  foreach(0..$#ARGV){
    print "$_: $ARGV[$_]\n";
  }
  die `perldoc $0`;
}
GetOptions(
            "fa_ref=s"      => \$fa_ref,
	    "gb_ref=s"      => \$gb_ref,
	    "draft=s{1,}"   => \@drafts,
	    "p=i"           => \$nb_threads,
	    "res_dir=s"     => \$res_dir,
	    "test"          => sub { $b_test = 1 }
	    );
# **********************************************************************

if($b_test)
{
	my $test_dir = './test_NGS_MauveReorder_on_ref/';
	my $fa_ref   = $test_dir.'ref_Campylobacter_jejuni_complete_genome.fasta';
	my $draft    = $test_dir.'AC1363_LargeContigs_out_AllStrains.unpadded.fasta';
	my $res_dir  = $test_dir.'res/';
	my $nb_threads = 8;

        -e $res_dir or mkdir $res_dir;

	my $prefix = 'test_NGS_MauveReorder_on_ref';
	my $cmd = join(' ', 'perl ./NGS_MauveReorder_on_ref.pl',
						'-fa_ref',   $fa_ref,
						'-draft',    $draft,
						'-p',        $nb_threads,
						'-res_dir',  $res_dir
	);
	my @res = ();
	print "$prog_tag [TEST] cmd:$cmd, line ".__LINE__."\n";
	@res = `$cmd`;
	print join("\n", "$prog_tag [TEST] ran, res:", @res, '');
	exit;
}

# **********************************************************************
# to get short name of a file, keeping extension
# **********************************************************************
sub short_name_with_ext($)
{
    my($sh) = @_;
    $sh =~ s/^.+\///;
    return $sh;
}
# **********************************************************************    

# defined $run_name or die "$prog_tag [Error] run_name must be provided, line ".__LINE__."\n";

if( (not defined $fa_ref)and
	(not defined $gb_ref))
{
	die "$prog_tag [Error] Neither fa_f nor gb_ref defined, you must provide at least one of thiese files, line ".__LINE__."\n";
}

foreach($fa_ref, $gb_ref, @drafts)
{
	defined $_ or next;
	-e $_ or die "$prog_tag [Error] $_ file does not exist, line ".__LINE__."\n";
}

defined $res_dir or die "$prog_tag [Error] $res_dir not defined, please check, line ".__LINE__."\n";
-e $res_dir or die "$prog_tag [Error] $res_dir does not exist, please check, line ".__LINE__."\n";

my $mauve_reorder_dir      = undef;
$mauve_reorder_dir = $res_dir;
-e $mauve_reorder_dir or mkdir $mauve_reorder_dir;

foreach my $draft(@drafts)
{
	$kmer_size = undef;
	
	# if kmer size forced, we find it in directory name
	if($draft =~ /_kmer(\d+)\//)
	{
		$kmer_size = $1;
		print "$prog_tag kmer_size:$kmer_size, line ".__LINE__."\n";
	}
	else
	{
		print "$prog_tag kmer_size deduced from kmergenie (no kmer keyword in draft directory name), line ".__LINE__."\n";
	}
	# ------------------------------------------------------------------
	
	# deduce strict sample name for mira
	# $sample_name =~ s/^([^_]+)_.*?$/$1/;
	defined $kmer_size and $sample_name .= "_kmer${kmer_size}";
	print "$prog_tag sample_name:$sample_name, line ".__LINE__."\n";
	# ------------------------------------------------------------------

	if(defined $res_dir)
	{
		if($res_dir !~ /\/$/)
		{
			$res_dir .= '/';
		}
	}
	else
	{
		$res_dir = $mauve_reorder_dir.$sample_name.'/';
		-e $res_dir or mkdir $res_dir;
		$res_dir .= 'reordered/';
	}
	-e $res_dir and remove_tree($res_dir); 
        mkdir $res_dir;
	
	defined $gb_ref or $gb_ref = $fa_ref;

	# reorder with Mauve
	# my $cmd = "java -Xmx500m -cp $mauve_jar org.gel.mauve.contigs.ContigOrderer -output $res_dir -ref $gb_ref -draft $draft";
        my $cmd = "$mauve_jar -output $res_dir -ref $gb_ref -draft $draft";
	print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
	my @res = `$cmd`;
	print(join("\n",'$prog_tag ran:', @res, ''));

	# get result directory names (alignment\d+), the last is the good one
	my @alignment_dirs = glob("${res_dir}alignment*");
	@alignment_dirs = nsort @alignment_dirs;
	if(scalar(@alignment_dirs) == 0)
	{
		die "$prog_tag [Error] No alignment directory found (res of Mauve reordering) using glob ${res_dir}alignment* in $res_dir, line ".__LINE__."\n";
	}
	my $res_align_f = $res_dir.short_name_with_ext($alignment_dirs[$#alignment_dirs]);
	my $res_fa_f    = $res_dir.short_name_with_ext($draft);
	
	# copy alignment data in result directory
	print "$prog_tag We copy $alignment_dirs[$#alignment_dirs]/* in $res_dir, line ".__LINE__."\n";
	# get files we want to copy
	my @cp_f = glob("$alignment_dirs[$#alignment_dirs]/*");
	foreach(@cp_f)
	{
		# real copy
		cp($_,$res_dir);
	}
	# remove useless alignment directories
	print join("\n",	"$prog_tag We remove:",
						@alignment_dirs,
						"directories line ".__LINE__."\n");
	remove_tree(@alignment_dirs);

	# ------------------------------------------------------------------
	# final progressiveMauve alignement of reordered contigs on the ref genome

	# then realign to ref genome with progressiveMauve to know which contigs align and which not (the last for those)
	# to know for which contig we need species identification
	my $name_of_alignment_out_f = $res_fa_f;
	$name_of_alignment_out_f =~ s/\.(?:fa|fasta)$/_vs_ref.xmfa/;

	if($name_of_alignment_out_f eq $res_fa_f)
	{
		die "$prog_tag [Error] $res_fa_f file does not end with fa or fasta extension, needed to create xmfa output file name for progressiveMauve line ".__LINE__."\n";
	}
		
	$cmd = "progressiveMauve --seed-family --output=$name_of_alignment_out_f $gb_ref $res_fa_f";
	print "$prog_tag cmd:$cmd, line ".__LINE__."\n";
	# run progressiveMauve alignment
	my $prefix = 'progressiveMauve_'.$sample_name.'_';
	@res = `$cmd`;
        print "$prog_tag launched to create $name_of_alignment_out_f alignment file\n";
        print(@res);
	# ------------------------------------------------------------------
}
