#from cyvcf2 import VCF
import vcf
import sys
import argparse
import pysam
from math import log
from math import log10
from math import exp
from collections import Counter

parser = argparse.ArgumentParser(description='''
POOHA: Parent Of Origin Haplotype Annotator ---
Assign parent of origin to heterozygous variants in a child using read-backed phasing
''')

parser.add_argument('vcf', type=str,
                    help='VCF-file of called variants.')
parser.add_argument('father', type=str,
                    help='Father of this family.')
parser.add_argument('mother', type=str,
                    help='Mother of this family')
parser.add_argument('child', type=str,
                    help='Child of this family')
parser.add_argument('child_bam', type=str,
                    help='Bam file with the reads from the child')
parser.add_argument('--verbose', action='store_true',
                    help="print progress output")
parser.add_argument('--min-parents-GQ',
                    default=20, action='store', type=int,
                    help='Minimum GQ for variants in parents when assigning obvious phase.')
parser.add_argument('--min-child-GQ',
                    default=10, action='store', type=int,
                    help='Minimum GQ for variants in child to be concidered.')
parser.add_argument('--region',
                    default="", action='store', type=str,
                    help='Analyse only the specified region. Either a chromosome: "chr1" or an interval: "chr1:1000-2000"')
parser.add_argument('--max-marker-distance',
                    default=2000, action='store', type=int,
                    help='Maximal distance between two het-markers that '
                    'will be considered.')
parser.add_argument('--conflicting_pairs', type=argparse.FileType('w'), default=None,
                    help='Name of file to write conflicting pairs of variants in')
parser.add_argument('--output_variants', choices=['all_het', 'all_mut', 'germline', 'somatic'],
                    help='Which variants should be output. all_het: all heterozygous variants in child. all_mut: all mutations (germline and somatic).', default='germline')

args = parser.parse_args()

if args.region == "":
    reader = vcf.Reader(filename=args.vcf)
elif ':' in args.region:
    chrom = args.region.split(":")[0]
    pos1, pos2 = (int(x) for x in args.region.split(":")[1].split('-'))
    reader = vcf.Reader(filename=args.vcf).fetch(chrom, pos1, pos2)
else:
    reader = vcf.Reader(filename=args.vcf).fetch(args.region)

bamfile = pysam.Samfile(args.child_bam, 'rb')

def get_pileup_column(site):
    '''Extract a specific column of the alignment and get
    all pileup reads there. Extract the query name and the called
    nucleotide'''
    pileup = bamfile.pileup(site[0], site[1]-1, site[1])
    len_diff = len(site[3]) - len(site[2])
    reflen = site[4]
    if len_diff == 0:
        assert reflen == 1
        for col in pileup:
            if col.pos == site[1]-1:
                return dict((read.alignment.qname, (read.alignment.seq[read.query_position],
                                                    read.alignment.mapping_quality,
                                                    read.alignment.is_read1))
                            for read in col.pileups if not read.is_del and
                            not read.query_position is None and
                            (read.alignment.seq[read.query_position] == site[2] or
                             read.alignment.seq[read.query_position] == site[3]) and
                            read.alignment.mapping_quality > 0)
    else:
        for col in pileup:
            # skal jeg se paa positionen efter naar jeg ser paa indels?
            if col.pos == site[1]-1:
                return dict((read.alignment.qname, (site[2] if read.indel == (len(site[2])-reflen)
                                                    else site[3],
                                                    read.alignment.mapping_quality,
                                                    read.alignment.is_read1))
                            for read in col.pileups if
                            (read.indel == len(site[2]) - reflen or
                             read.indel == len(site[3]) - reflen) and
                            read.alignment.mapping_quality > 0)
    return {}

def LGQ(gt):
    L = gt.data.PL
    L.sort()
    return L[1]

def VAF(gt):
    L = gt.data.AD
    if L[1] == 0:
        return 0.0
    else:
        return float(L[1])/(L[0]+L[1])

def phred_2_prob(Q):
    return 10**(-float(Q)/10)

def prob_2_phred(P):
    if P<1e-100:
        return 1000
    else:
        return -10*log10(P)

def trio_qual(Q1, Q2, Q3):
    return ((1 - phred_2_prob(Q1))*(1 - phred_2_prob(Q2))*(1 - phred_2_prob(Q3)))

def get_obvious_PoO():
    assigned_hets = []
    unassigned_hets = []
    PoO = {}
    call_qual = {}
    last_chrom = 'na'
    last_pos = 0
    n_som = 0
    n_germ = 0
    for record in reader:
        try:
            # or ((record.FILTER is not None) and len(record.FILTER)>0)):
            if (record is None):
                continue
            if (len(assigned_hets) + len(unassigned_hets)) > 0 and (last_chrom != record.CHROM or (record.POS - last_pos) > args.max_marker_distance):
                if args.verbose:
                    print("Read to", record.CHROM, record.POS, file=sys.stderr)
                yield assigned_hets, unassigned_hets, PoO, call_qual, n_som, n_germ
                assigned_hets = []
                unassigned_hets = []
                PoO = {}
                call_qual = {}
                n_som = 0
                n_germ = 0
            gt_father = record.genotype(args.father)
            gt_mother = record.genotype(args.mother)
            gt_child = record.genotype(args.child)
            if not gt_child.is_het:
                continue
            father_alleles = gt_father.gt_bases.split('/')
            mother_alleles = gt_mother.gt_bases.split('/')
            child_a1, child_a2 = gt_child.gt_bases.split('/')
            if (not record.is_snp) and len(child_a1) == len(child_a2) and len(child_a1) != 1:
                continue

            var = (record.CHROM, record.POS, child_a1,
                   child_a2, len(record.REF))

            if LGQ(gt_child) < args.min_child_GQ:
                continue

            qual = trio_qual(LGQ(gt_father), LGQ(gt_mother), LGQ(gt_child))
            qual = qual + (1.0-qual)*0.5
            if min(LGQ(gt_father), LGQ(gt_mother)) < args.min_parents_GQ or qual < 0.55:
                unassigned_hets.append(var)
                call_qual[var] = ('low_parent_GQ', qual)
                continue

            if not ((child_a1 in father_alleles and child_a2 in mother_alleles) or
                    (child_a2 in father_alleles and child_a1 in mother_alleles)):
                if gt_father.gt_type == 0 and gt_mother.gt_type == 0:
                    if VAF(gt_child) < 0.3:
                        call_qual[var] = ('ref_somatic', qual)
                        n_som += 1
                    else:
                        call_qual[var] = ('ref_denovo', qual)
                        n_germ += 1
                else:
                    if VAF(gt_child) > 0.7:
                        call_qual[var] = ('alt_somatic', qual)
                        n_som += 1
                    else:
                        call_qual[var] = ('alt_denovo', qual)
                        n_germ += 1
                unassigned_hets.append(var)
            else:
                call_qual[var] = ('no_mendelian_violation', qual)
                if (father_alleles.count(child_a2) < mother_alleles.count(child_a2) or
                        father_alleles.count(child_a1) > mother_alleles.count(child_a1)):
                    assigned_hets.append(var)
                    try:
                        LR = max(log(1.0-qual) - log(qual), -500.0)
                    except ValueError:
                        LR = -500.0
                    PoO[var] = ('M', 'obvious', LR, 1.0-qual, qual, 0, 0)
                elif (father_alleles.count(child_a2) > mother_alleles.count(child_a2) or
                      father_alleles.count(child_a1) < mother_alleles.count(child_a1)):
                    assigned_hets.append(var)
                    try:
                        LR = min(log(qual) - log(1.0-qual), 500)
                    except ValueError:
                        LR = 500.0
                    PoO[var] = ('P', 'obvious', LR, qual, 1.0-qual, 0, 0)
                else:
                    unassigned_hets.append(var)
            if var not in PoO:
                PoO[var] =  ('U', 'unassigned', 0.0, 0, 0, 0, 0)
            last_chrom = record.CHROM
            last_pos = record.POS
        except AttributeError:
            # if args.verbose:
            #    print "Unexpected error:", sys.exc_info()
            continue
        except TypeError:
            # if args.verbose:
            #    print "Unexpected error:", sys.exc_info()
            continue
    if (len(assigned_hets) + len(unassigned_hets)) > 0:
        if args.verbose:
            print("Read to", last_chrom, last_pos, file=sys.stderr)
        yield assigned_hets, unassigned_hets, PoO, call_qual, n_som, n_germ


def assign_PoO(assigned_hets, unassigned_hets, PoO, call_qual):
    start_j = 0
    unassigned_hets.sort()
    assigned_hets.sort()
    new_unassigned_hets = []
    new_assigned_hets = []
    for i in range(len(unassigned_hets)):
        if unassigned_hets[i] not in pileups:
            pileups[unassigned_hets[i]] = get_pileup_column(unassigned_hets[i])
        i_reads = pileups[unassigned_hets[i]]

        var_i = unassigned_hets[i]
        i_chrom, i_pos, i_ref, i_alt, i_reflen = var_i

        # Collect the j indices relevant for this i index.
        valid_j_indices = []

        for j in range(start_j, len(assigned_hets)):
            var_j = assigned_hets[j]
            if (assigned_hets[j][0] < i_chrom or
                    (assigned_hets[j][0] == i_chrom and i_pos - assigned_hets[j][1] > args.max_marker_distance)):
                start_j += 1
                continue

            var_j_type = call_qual[var_j][0]

            if assigned_hets[j][0] > i_chrom:
                break  # different chromosome, no need to continue here

            if assigned_hets[j][1] - i_pos > args.max_marker_distance:
                break

            valid_j_indices.append(j)

        var_i_type = call_qual[var_i][0]
        total_reads = 0
        n_j = 0
        log_LR_x = 0.0

        # now phase the i/j sites.
        for j in valid_j_indices:
            if assigned_hets[j] not in pileups:
                pileups[assigned_hets[j]] = get_pileup_column(assigned_hets[j])
            j_reads = pileups[assigned_hets[j]]

            names_overlap = set(i_reads.keys()).intersection(j_reads.keys())

            if len(names_overlap) == 0:
                continue

            j_chrom, j_pos, j_ref, j_alt, j_reflen = assigned_hets[j]
            n_j += 1
            P_j_pat = PoO[assigned_hets[j]][3]
            P_j_mat = PoO[assigned_hets[j]][4]

            total_reads += len(names_overlap)
            # probability of data given that o(i)==o(j)
            log_P_reads_given_oi_eq_oj = 0.0
            # probability of data given that o(i)!=o(j)
            log_P_reads_given_oi_neq_oj = 0.0
            # log_P_reads_given_i_pat_j_pat = 0.0
            # log_P_reads_given_i_mat_j_pat = 0.0
            # log_P_reads_given_i_pat_j_mat = 0.0
            # log_P_reads_given_i_mat_j_mat = 0.0
            n_eq = 0
            n_neq = 0
            n_aa = 0 # alternative allele at both i and j
            n_rr = 0 # refernce allele at both i and j
            n_ar = 0 # refernce allele at i and alternative allele at j
            n_ra = 0 # refernce allele at i and alternative allele at j

            for qname in names_overlap:
                i_allele, i_mapq, i_is_read1 = i_reads[qname]
                j_allele, j_mapq, j_is_read1 = j_reads[qname]

                if i_is_read1 == j_is_read1 and i_mapq == j_mapq:
                    #assert i_mapq == j_mapq
                    P_r_correct = (1.0 - phred_2_prob(i_mapq))
                else:
                    P_r_correct = (1.0 - phred_2_prob(i_mapq))*(1.0 - phred_2_prob(j_mapq))

                if i_allele == i_ref and not var_i_type == 'ref_somatic':
                    if j_allele == j_ref and not var_j_type == 'ref_somatic':
                        n_rr += 1
                        # i_ref phased med j_ref
                        assert P_r_correct + (1.0 - P_r_correct) * 0.5 >= 0.5
                        n_eq += 1
                        log_P_reads_given_oi_eq_oj += log(
                            P_r_correct + (1.0 - P_r_correct) * 0.5)
                        log_P_reads_given_oi_neq_oj += log(
                            (1.0 - P_r_correct) * 0.5)
                    elif j_allele == j_alt:
                        if call_qual[var_j][0] == 'ref_somatic':
                            continue
                        n_ra += 1
                        # i_ref phased med j_alt
                        n_neq += 1
                        log_P_reads_given_oi_neq_oj += log(
                            P_r_correct + (1.0 - P_r_correct) * 0.5)
                        log_P_reads_given_oi_eq_oj += log(
                            (1.0 - P_r_correct) * 0.5)
                elif i_allele == i_alt:
                    if j_allele == j_ref and not var_j_type == 'ref_somatic':
                        n_ar += 1
                        # i_alt phased med j_ref
                        n_neq += 1
                        log_P_reads_given_oi_neq_oj += log(
                            P_r_correct + (1.0 - P_r_correct) * 0.5)
                        log_P_reads_given_oi_eq_oj += log(
                            (1.0 - P_r_correct) * 0.5)
                    elif j_allele == j_alt:
                        # i_alt phased med j_alt
                        n_aa += 1
                        n_eq += 1
                        log_P_reads_given_oi_eq_oj += log(
                            P_r_correct + (1.0 - P_r_correct) * 0.5)
                        log_P_reads_given_oi_neq_oj += log(
                            (1.0 - P_r_correct) * 0.5)

            if log_P_reads_given_oi_neq_oj > 500:
                log_P_reads_given_oi_neq_oj = 500
            elif log_P_reads_given_oi_neq_oj < -500:
                log_P_reads_given_oi_neq_oj = -500
            P_reads_given_oi_neq_oj = exp(log_P_reads_given_oi_neq_oj)

            if log_P_reads_given_oi_eq_oj > 500:
                log_P_reads_given_oi_eq_oj = 500
            elif log_P_reads_given_oi_eq_oj < -500:
                log_P_reads_given_oi_eq_oj = -500
            P_reads_given_oi_eq_oj = exp(log_P_reads_given_oi_eq_oj)

            if (n_aa+n_rr) > (n_ar+n_ra):
                if n_ar == 0 and n_ra > 0:
                    if var_i_type == 'ref_denovo':
                        likely_somatic_reads[var_i] += 1
                elif n_ar > 0 and n_ra == 0:
                    if var_j_type == 'ref_denovo':
                        likely_somatic_reads[var_j] += 1
                elif n_ar > 0 and n_ra > 0:
                    conflicting_pairs.append(((i_chrom, i_pos, i_alt), (j_chrom, j_pos, j_alt), min(n_ar, n_ra)))
            elif (n_aa+n_rr) < (n_ar+n_ra):
                if n_aa == 0 and n_rr > 0:
                    # TODO: Tjek disse if sætnigner
                    if var_i_type == 'ref_denovo' and var_j_type == 'no_mendelian_violation':
                        likely_somatic_reads[var_i] += 1
                    elif var_j_type == 'ref_denovo' and var_i_type == 'no_mendelian_violation':
                        likely_somatic_reads[var_j] += 1
                elif n_aa > 0 and n_rr == 0:
                    # TODO: Tjek disse if sætnigner
                    if var_i_type == 'ref_denovo' and var_j_type == 'no_mendelian_violation':
                        likely_somatic_reads[var_i] += 1
                    elif var_j_type == 'ref_denovo' and var_i_type == 'no_mendelian_violation':
                        likely_somatic_reads[var_j] += 1
                elif n_ar > 0 and n_ra > 0:
                    conflicting_pairs.append(((i_chrom, i_pos, i_alt), (j_chrom, j_pos, j_alt), min(n_ar, n_ra)))

            try:
                log_P_reads_given_i_pat = log(P_j_pat * P_reads_given_oi_eq_oj +
                                              P_j_mat * P_reads_given_oi_neq_oj)
                log_P_reads_given_i_mat = log(P_j_pat * P_reads_given_oi_neq_oj +
                                              P_j_mat * P_reads_given_oi_eq_oj)
            except ValueError:
                if args.verbose:
                    print(unassigned_hets[i], file=sys.stderr)
                    print(n_eq, n_neq, file=sys.stderr)
                    print(log_P_reads_given_oi_neq_oj,
                          log_P_reads_given_oi_eq_oj, file=sys.stderr)
                    print(n_j, len(names_overlap),
                          P_r_correct, file=sys.stderr)
                    print(P_reads_given_oi_neq_oj, P_reads_given_oi_eq_oj,
                          P_j_pat, P_j_mat, file=sys.stderr)
                continue
            log_LR_x += log_P_reads_given_i_pat - log_P_reads_given_i_mat

        if log_LR_x > 500:
            LR_x = exp(500)
        else:
            LR_x = exp(log_LR_x)

        if total_reads == 0 or log_LR_x == 0:
            new_unassigned_hets.append(unassigned_hets[i])
            PoO[unassigned_hets[i]] = (
                'U', 'unassigned', 0.0, 0, 0, total_reads, n_j)
        elif log_LR_x > 0.0:
            new_assigned_hets.append(unassigned_hets[i])
            PoO[unassigned_hets[i]] = (
                'P', 'phased', log_LR_x, LR_x/(1+LR_x), 1.0 - (LR_x/(1+LR_x)), total_reads, n_j)
        else:
            new_assigned_hets.append(unassigned_hets[i])
            PoO[unassigned_hets[i]] = (
                'M', 'phased', log_LR_x, LR_x/(1+LR_x), 1.0 - (LR_x/(1+LR_x)), total_reads, n_j)

    return len(unassigned_hets) - len(new_unassigned_hets), assigned_hets + new_assigned_hets, new_unassigned_hets


pileups = {}
conflicting_pairs = []
likely_somatic_reads = Counter()

print('chrom', 'pos', 'ref', 'alt', 'child', 'ParentOfOrigin', 'method', 'PoO_qual', 'variant_type', 'n_informative_reads', 'n_informative_vars', 'n_likely_somatic_reads')

if args.verbose:
    print("Assigning PoO to variants with obvious PoO...", file=sys.stderr)
for assigned_hets, unassigned_hets, PoO, call_qual, n_som, n_germ in get_obvious_PoO():
    if args.verbose:
        print(len(assigned_hets),
              "heterozygous variants with obvious PoO", file=sys.stderr)
        print(len(unassigned_hets),
              "heterozygous variants without obvious PoO", file=sys.stderr)
    if (len(assigned_hets) > 0 and
        ((args.output_variants == 'all_het') or
         (args.output_variants == 'somatic' and n_som>0) or
         (args.output_variants == 'germline' and n_germ>0) or
         (args.output_variants == 'all_mut' and (n_germ>0 or n_som>0)))):
        n_new = 1
        i = 1
        while n_new > 0  and len(unassigned_hets) > 0:
            n_new, assigned_hets, unassigned_hets = assign_PoO(
                assigned_hets, unassigned_hets, PoO, call_qual)
            if args.verbose:
                print('Assigning PoO using read-backed phasing, round',
                      i, '...', file=sys.stderr)
                print('assigned PoO for ', n_new,
                      ' new variants.', file=sys.stderr)
                print("len(assigned_hets) =", len(assigned_hets), file=sys.stderr)
                i += 1
    if args.verbose:
        print('DONE', file=sys.stderr)
    if not args.conflicting_pairs is None:
        for v1, v2, nm in conflicting_pairs:
            print('_'.join(str(x) for x in v1), '_'.join(str(x) for x in v2), nm, file=args.conflicting_pairs)
    variants = sorted(PoO)
    for var in variants:
        parent, method, log_LR_x, P_prob, M_prob, nreads, n_j = PoO[var]

        if method == 'unassigned':
            Q = 0
        else:
            Q = prob_2_phred(min(P_prob, M_prob))

        mtype, cqual = call_qual[var]
        if(args.output_variants == 'all_het' or
            (('denovo' in mtype) and args.output_variants in ['all_mut', 'germline']) or
                (('somatic' in mtype) and args.output_variants in ['all_mut', 'somatic'])):
            print(var[0], var[1], var[2], var[3], args.child, parent,
                  method, int(Q), mtype, nreads, n_j, likely_somatic_reads[var])
    pileups = {}
    conflicting_pairs = []
    likely_somatic_reads = Counter()
