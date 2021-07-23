#!/usr/bin/env python3

import os.path
import tempfile
import subprocess
from Bio import SeqIO


"""
-------------------------------------------------------------------------------
PGA (Perl)
Project: https://github.com/quxiaojian/PGA
File: PGA.pl, lines ~1550:
system ("makeblastdb -in $IR_temp_random -hash_index -dbtype nucl -blastdb_version 4 -logfile $screenlog");
system ("blastn -task blastn -query $IR_temp_random -db $IR_temp_random -outfmt 6 -perc_identity 99 -out $blast_IR_temp_random");
...
open (my $input_IR,"<",$blast_IR_temp_random);
my %IR;
while (<$input_IR>) {
    $_=~ s/\r|\n//g;
    my ($c1,$c2,$c3,$ir_length,$c5,$c6,$qs,$qe,$ss,$se,$c11,$c12)=split /\t/,$_;
    if (($ir_length != 2 * $length_cp) and ($ir_length != $length_cp) and ($ir_length >= $inverted_repeat) and ($qe <= $length_cp) and ($se <= $length_cp) and ($qs < $ss)) {
        $IR{$ir_length}=$qs."\t".$qe."\t".$ss."\t".$se;
    }
}
close $input_IR;
unlink ("$blast_IR_temp_random");
my (@IR_length,@boundary);
foreach my $key (sort {$b <=> $a} keys %IR) {
    push @IR_length,$key;
    push @boundary,$IR{$key};
}
my ($AA,$BB,$CC,$DD)=split /\t/,$boundary[0];
my ($JLB,$JSB,$JLA,$JSA);
if ($CC <= $length_cp) {
    $JLB=$AA;
    $JSB=$BB;
    $JLA=$CC;
    $JSA=$DD;
}elsif ($CC > $length_cp) {
    $JLB=$AA;
    $JSB=$BB;
    $JLA=$CC-$length_cp;
    $JSA=$DD;
}


Clearer form:
blastn -query <input_fasta_file> -subject <input_fasta_file> -perc_identity 99 -outfmt <fmt>
# Treba li -task blastn
Output process ...


-------------------------------------------------------------------------------
plann (Perl)
Project: https://github.com/daisieh/plann
File plann.pl, lines ~100:

system("blastn -query $fastafile -subject $fastafile -outfmt 5 -out $refblast -evalue 1e-200");

my $self_array = Blast::parse_xml ("$refblast");
my @irs = ();
foreach my $hit (@$self_array) {
    my @hsps = sort Blast::sort_regions_by_start @{$hit->{'hsps'}};
    foreach my $hsp (@hsps) {
        my $querylen = $hsp->{'query-to'} - $hsp->{'query-from'};
        # IRs are between 10,000 and 50,000 bp and are inverted.
        if (($querylen < 50000) && ($querylen > 10000) && ($hsp->{'hit-frame'} == -1)) {
            push @irs, $hsp;
        }
    }
}

if (@irs > 2) {
    print "Warning! There seem to be more than two inverted repeats (".@irs." found). Are you sure this is a plastome sequence?\n";
}
@irs = sort Blast::sort_hsps_by_query_start @irs;

# all of the IR information is in one of the IRs, so shift.
my $ir = shift @irs;

Clearer form:
blastn -query <input_fasta_file> -subject <input_fasta_file> -evalue 1e-200
"""

VARIANT_2_BLASTN_ARGS = dict(
    pga=('-perc_identity', '99'),
    plann=('-evalue', '1e-200')
)


def self_blast(variant, seq_filename, min_ir_length=1000, print_blast_output=False, leave_tmp_file=False):
    # Note: seq_filename can be in fasta or genbank format
    if any(seq_filename.endswith(e) for e in ('.fa', '.fas', '.fasta', '.fs')):
        fa_filename = seq_filename
    else:
        # Convert file into fasta format
        assert any(seq_filename.endswith(e) for e in ('.gb', '.genbank', '.gbs')), seq_filename
        fa_filename = os.path.join(tempfile.gettempdir(), f'{next(tempfile._get_candidate_names())}.fa')
        SeqIO.convert(seq_filename, 'genbank', fa_filename, 'fasta')

    # 5 XML output format
    # 6 tabular output format
    outfmt = '6 qstart qend sstart send sstrand length nident mismatch'
    cmd = ['blastn', '-query', fa_filename, '-subject', fa_filename, '-outfmt', outfmt]
    cmd.extend(VARIANT_2_BLASTN_ARGS[variant])
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return

    res_irs = None
    output = result.stdout.decode('utf-8')
    if print_blast_output:
        print(output)
    lines = [(int(f[0]), int(f[1]), int(f[2]), int(f[3]), f[4], int(f[5]), int(f[6]), int(f[7]))
             for l in output.split('\n') if (f := l.split())]

    if variant == 'pga':
        # Hack to find sequence length. First match is for sure whole sequence.
        f_line = lines[0]
        if f_line[0] == f_line[2] == 1 and \
           f_line[1] == f_line[3] == f_line[5] == f_line[6] and \
           f_line[4] == 'plus' and f_line[7] == 0:
            seq_length = f_line[1]
        else:
            seq_length = len(SeqIO.read(seq_filename, 'genbank').seq)

        irs = []
        for line in lines:
            qs, qe, ss, se, _, ir_length, _, _ = line
            if ir_length != 2 * seq_length and \
               ir_length != seq_length and \
               ir_length >= min_ir_length and \
               qe <= seq_length and \
               se <= seq_length and \
               qs < ss:
                irs.append((ir_length, qs, qe, ss, se))
        if irs:
            _, AA, BB, CC, DD = max(irs)  # Max by ir_length
            AA -= 1  # Our indexing
            DD -= 1
            res_irs = ((AA, BB), (DD, CC)) if CC <= seq_length else ((AA, BB), (DD, CC - seq_length))

    elif variant == 'plann':
        irs = []
        for line in lines:
            qs, qe, ss, se, strand, ir_length, _, _ = line
            querylen = qe - qs
            if querylen < 50000 and \
               querylen > 10000 and \
               strand == 'minus':
                irs.append((qs, qe, ss, se))
        qs, qe, ss, se = min(irs)  # sort Blast::sort_hsps_by_query_start @irs;
        res_irs = ((qs - 1, qe), (se - 1, ss))

    else:
        assert False, variant

    #
    if fa_filename != seq_filename:
        if leave_tmp_file:
            print(f'Note: tmp file {fa_filename} is not removed!')
        else:
            os.remove(fa_filename)

    return res_irs


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="Run self-blast IR detection on sequence stored as fasta or genbank file.")
    parser.add_argument('variant', choices=tuple(VARIANT_2_BLASTN_ARGS.keys()), help='Annotation variant')
    parser.add_argument('seq_filename', help='Sequence filename')
    parser.add_argument('-m', '--min-ir-length', default=1000, help='Minimal ir length for PGA variant.')
    parser.add_argument('-B', '--print-blast-output', action='store_true', help='Print Blastn output.')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')

    params = parser.parse_args()
    print(self_blast(params.variant, params.seq_filename,
                     min_ir_length=params.min_ir_length,
                     print_blast_output=params.print_blast_output,
                     leave_tmp_file=params.leave_tmp_file))
