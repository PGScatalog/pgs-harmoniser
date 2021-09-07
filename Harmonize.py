import argparse, os, sys, gzip
from tqdm import tqdm
from pgs_harmonizer.harmonize import *
from datetime import datetime

# Inputs
parser = argparse.ArgumentParser(
    description='Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome build.')
subparsers = parser.add_subparsers(help='Harmonization Commands help', dest='HmAction')

# Sub-parser for HmPOS
parser_POS = subparsers.add_parser('HmPOS', help='HmPOS - Harmonizing position information (adding/updating chr/pos information)')
parser_POS.add_argument(dest="pgs_id", help="PGS Catalog Score ID", metavar="PGS######", type=str)
parser_POS.add_argument(dest="target_build", help="Target genome build choices: 'GRCh37'or GRCh38'", metavar="GRCh3#",
                        choices=['GRCh37', 'GRCh38'])
parser_POS.add_argument("-loc_files", dest="loc_scorefiles",
                        help="Root directory where the PGS files are located, otherwise assumed to be in: ../pgs_ScoringFiles/",
                        metavar="DIR",
                        default='../pgs_ScoringFiles/', required=False)
parser_POS.add_argument("-source_build", dest="source_build",
                        help="Source genome build [overwrites information in the scoring file header]",
                        metavar="GENOMEBUILD",
                        default=None, required=False)
parser_POS.add_argument("-loc_hmoutput", dest="loc_outputs",
                        help="Directory where the harmonization output will be saved (default: PGS_HmPOS/)",
                        metavar="DIR",
                        default='./PGS_HmPOS/', required=False)
parser_POS.add_argument('--var2location',
                        help='Root directory to the annotations from the var2location.pl script (ENSEMBL SQL connection)',
                        metavar="DIR", required=False)
parser_POS.add_argument('--silent_tqdm', help='Disables tqdm progress bar',
                        action='store_true', required=False)
parser_POS.add_argument('--ignore_rsid', help='Ignores rsID mappings and harmonizes variants using only liftover',
                        action='store_true', required=False)
parser_POS.add_argument('--gzip', help='Writes gzipped harmonized output',
                        action='store_true', required=False)

# Sub-parser for HmVCF
parser_VCF = subparsers.add_parser('HmVCF', help='HmVCF - Checking positional information and/or adding other_alleles')
parser_VCF.add_argument(dest="pgs_id", help="PGS Catalog Score ID", metavar="PGS######", type=str)
parser_VCF.add_argument("-loc_files", dest="loc_scorefiles",
                        help="Root directory where the PGS files are located, otherwise assumed to be in: PGS_HmPOS/",
                        metavar="DIR",
                        default='./PGS_HmPOS/', required=False)
parser_VCF.add_argument("-loc_hmoutput", dest="loc_outputs",
                        help="Directory where the harmonization output will be saved (default: PGS_HmPOS/)",
                        metavar="DIR",
                        default='./PGS_HmVCF/', required=False)
parser_VCF.add_argument(dest="target_build",
                        help="Target genome build choices: 'GRCh37'or GRCh38'",
                        metavar="GRCh3#",
                        choices=['GRCh37', 'GRCh38'])
parser_VCF.add_argument("-loc_vcfs", dest="loc_vcfref",
                        help="Directory where the VCF files are located, otherwise assumed to be in: map/vcf_ref/",
                        metavar="DIR",
                        default='map/vcf_ref/'),
parser_VCF.add_argument("-cohort_vcf", dest="cohort_name",
                        help="Cohort VCF: Used to check if a variant is present in the genotyped/imputed variants for "
                             "a cohort and add other allele when the information from ENSEMBL is ambiguous "
                             "(multiple potential alleles)",
                        metavar="COHORT",
                        default=None, required=False)
parser_VCF.add_argument('--addOtherAllele',
                        help='Adds a other_allele(s) column for PGS that only have a recorded effect_allele',
                        action='store_true', required=False)
parser_VCF.add_argument('--addVariantID',
                        help='Returns a column with the ID from the VCF corresponding to the match variant/allele(s)',
                        action='store_true', required=False)
parser_VCF.add_argument('--author_reported',
                        help='Replaces unmappable variants (hm_code = -5) with the author-reported code (hm_code = 0)',
                        action='store_true', required=False)
parser_VCF.add_argument('--skip_strandflips',
                        help='This flag will stop the harmonizing from trying to correct strand flips',
                        action='store_true', required=False)
parser_VCF.add_argument('--silent_tqdm', help='Disables tqdm progress bar',
                        action='store_true', required=False)
parser_VCF.add_argument('--gzip', help='Writes gzipped harmonized output',
                        action='store_true', required=False)

args = parser.parse_args()

class HarmonizationError(Exception):
    """Base class for exceptions in this module."""
    pass

def variant_HmPOS(v, rsIDmaps=None, liftchain=None, isSameBuild=False, inferOtherAllele=False):
    """Finds Harmonized Position (HmPOS) information for a variant using Ensembl variation/liftover"""
    hm_source = ''  # {'Author-reported', 'ENSEMBL Variation', 'liftover' }
    hm_rsID = ''
    hm_chr = ''
    hm_pos = ''
    hm_inferOtherAllele = None  # Field to capture the inferred other/reference allele

    if rsIDmaps and (v['rsID'] in rsIDmaps):
        v_map = rsIDmaps.get(v['rsID'])
    else:
        v_map = None

    if v_map is not None:
        hm_chr, hm_pos, hm_alleles = list(v_map.select_canonical_data(chromosomes))
        hm_source = 'ENSEMBL'
        hm_rsID = v_map.id
        if inferOtherAllele:
            if 'other_allele' in v:
                if pd.isnull(v.get('other_allele')) is True:
                    hm_inferOtherAllele = v_map.infer_OtherAllele(v['effect_allele'])
            else:
                hm_inferOtherAllele = v_map.infer_OtherAllele(v['effect_allele'])  # Based on the rsID
    elif 'chr_name' and 'chr_position' in v:
        if isSameBuild:
            hm_chr = v['chr_name']
            hm_pos = v['chr_position']
            hm_source = 'Author-reported'  # Author-reported
        elif (liftchain is not None) and (liftchain.chain is not None):
            if (pd.isnull(v['chr_name']) is False) and (pd.isnull(v['chr_position']) is False):
                hm_chr, hm_pos, hm_liftover_multimaps = list(liftchain.lift(v['chr_name'], v['chr_position']))  # liftover
                hm_source = 'liftover'
    if all([x == '' for x in [hm_chr, hm_pos]]):
        hm_source = 'Unknown'

    if hm_pos != '':
        hm_pos = str(hm_pos)

    if inferOtherAllele:
        if hm_inferOtherAllele is None:
            hm_inferOtherAllele = ''
        return pd.Series([hm_source, hm_rsID, hm_chr, hm_pos, hm_inferOtherAllele])
    else:
        return pd.Series([hm_source, hm_rsID, hm_chr, hm_pos])


def run_HmPOS(args, chunksize=100000):
    # Module-specifc imports
    from pgs_harmonizer.ensembl_tools import ensembl_post, clean_rsIDs, parse_var2location
    from pgs_harmonizer.liftover_tools import liftover, map_release

    ## Set I/O File focations
    # Scoring file location
    if 'loc_scorefiles' in args:
        if not args.loc_scorefiles.endswith('/'):
            args.loc_scorefiles += '/'
        loc_scorefile = args.loc_scorefiles + args.pgs_id + '.txt.gz'
    else:
        loc_scorefile = '../pgs_ScoringFiles/{}.txt.gz'.format(args.pgs_id)
    # Define output location
    ofolder = args.loc_outputs
    if ofolder.endswith('/'):
        ofolder = ofolder[:-1]
    if os.path.isdir(ofolder) is False:
        os.mkdir(ofolder)
    loc_hm_out = '{}/{}_hmPOS{}.txt'.format(ofolder, args.pgs_id, args.target_build)
    if args.gzip is True:
        loc_hm_out += '.gz'

    # Read Score File
    print('Reading Score File')
    header, df_scoring = read_scorefile(loc_scorefile)
    tqdm.pandas()

    # Get consistent source build (e.g. NCBI/GRC)
    source_build = header['genome_build']
    source_build_mapped = None
    if args.source_build is not None:
        source_build = args.source_build

    if source_build in map_release.values():
        for grc, hg in map_release.items():
            if hg == source_build:
                source_build_mapped = grc
        print('PGS ID: {} | Build: {}/{}'.format(header['pgs_id'], source_build, source_build_mapped))
    elif source_build in map_release.keys():
        source_build_mapped = source_build
        print('PGS ID: {} | Build: {}'.format(header['pgs_id'], source_build))
    else:
        print('PGS ID: {} | Build: {}'.format(header['pgs_id'], source_build))
    print('Number of variants (score file lines) = {}'.format(header['variants_number']))
    # ToDo - print columns available for mapping

    # Sorting out the genome build
    mappable = False
    tf_unmappable2authorreported = False
    if source_build_mapped is None:
        if 'rsID' in df_scoring.columns:
            mappable = True
        else:
            sys.exit(
                '{} CAN NOT BE HARMONIZED: Need to specify the source genome build (-source_build)'.format(
                    header['pgs_id']))
            # ToDo possibly implement a method to guess the genome build using GRCh37/38 VCFs
    else:
        mappable = True

    if mappable is False:
        sys.exit(
            '{} CAN NOT BE HARMONIZED: Insufficient variant or genome build data for mapping'.format(header['pgs_id']))

    if args.target_build == source_build_mapped:
        print('Harmonizing -> {}'.format(args.target_build))
        isSameBuild = True
    else:
        print('Re-Mapping/Lifting + Harmonizing -> {}'.format(args.target_build))
        isSameBuild = False

    # Load Liftover Chains
    if source_build is not None:
        build_map = liftover(source_build, args.target_build)  # Get the chain file
    else:
        build_map = None

    # Source ENSEMBL DB/API variant mappings if required
    mapping_ensembl = None
    if 'rsID' in df_scoring.columns and args.ignore_rsid is False:
        tomap_rsIDs = clean_rsIDs(list(df_scoring['rsID']))
        use_rest = True
        if args.var2location:
            loc_var2location = args.var2location
            if os.path.isdir(loc_var2location):
                use_rest = False
                mapping_ensembl = parse_var2location(loc_var2location)
        if use_rest == True:
            print('Retrieving rsID mappings from ENSEMBL API')
            mapping_ensembl = ensembl_post(tomap_rsIDs, args.target_build)  # Retrieve the SNP info from ENSEMBL
    # Start Output

    if args.gzip is True:
        hm_out = gzip.open(loc_hm_out, 'wt')
    else:
        hm_out = open(loc_hm_out, 'w')
    # Append information to header:
    header.update({'HmPOS_build': args.target_build,
                   'HmPOS_date': str(
                       datetime.date(datetime.now()))})  # ToDo Consider adding information about the ENSEMBL build?
    hm_out.write('\n'.join(create_scoringfileheader(header)))
    hm_out.write('\n')

    hm_Passed = True
    hm_counts = {}
    hm_chunks = int(np.ceil(df_scoring.shape[0] / chunksize))
    while hm_Passed is True:
        for ic in tqdm(range(0, hm_chunks)):
            start = ic*chunksize
            end = start + chunksize
            try:
                df_chunk = df_scoring.iloc[start:end,:].copy()
                #print(start, end, df_chunk.index[0], df_chunk.index[-1])
                df_chunk[['hm_source', 'hm_rsID', 'hm_chr', 'hm_pos', 'hm_inferOtherAllele']] = df_chunk.apply(variant_HmPOS,
                                                                                                               axis=1,
                                                                                                               rsIDmaps=mapping_ensembl,
                                                                                                               liftchain=build_map,
                                                                                                               isSameBuild=isSameBuild,
                                                                                                               inferOtherAllele=True)
                for hm_source, hm_count in dict(df_chunk['hm_source'].value_counts()).items():
                    if hm_source in hm_counts:
                        hm_counts[hm_source] += hm_count
                    else:
                        hm_counts[hm_source] = hm_count
                if ic == 0:
                    df_chunk.to_csv(hm_out, mode='a', index=False, sep='\t')  # Write output using pandas
                else:
                    df_chunk.to_csv(hm_out, mode='a', index=False, header=False, sep='\t')  # Write output using pandas
            except:
                hm_Passed = False
        hm_Passed = 'COMPLETED'

    if hm_Passed == 'COMPLETED':
        hm_out.close()
        print('Mapped {} -> {}'.format(hm_counts, loc_hm_out))
        return
    else:
        hm_out.close()
        os.remove(loc_hm_out)
        print('FAILED')
        raise HarmonizationError
        return


def variant_HmVCF(v, vcfs_targetbuild, CohortVCF=None, returnOtherAllele=True):
    """Determines whether the variant maps correctly to a reference VCF"""
    hm_source = v['hm_source']  # {'Author-reported', 'ENSEMBL Variation', 'liftover' }
    hm_matchesVCF = False  # T/F whether the variant is consistent with the VCF/Variant Lookup
    hm_isPalindromic = False  # T/F whether the alleles are consistent with being palindromic
    hm_isFlipped = False  # T/F whether the alleles are consistent with the negative strand (from VCF)
    hm_vid = None
    hm_code = None  # Derived from the above True/False information

    # Sort out non-effect/other allele
    other_allele = None
    if 'other_allele' in v:
        if pd.isnull(v['other_allele']) is False:
            other_allele = v['other_allele']
    hm_inferOtherAllele = None  # Field to capture the inferred other/reference allele
    if 'hm_inferOtherAllele' in v:
        if pd.isnull(v['hm_inferOtherAllele']) is False:
            hm_inferOtherAllele = v['hm_inferOtherAllele']

    # Check/select allales
    if pd.isnull(v['hm_source']) is False:
        v_records = vcfs_targetbuild.vcf_lookup(chromosome=v['hm_chr'], position=v['hm_pos'], rsid=v['hm_rsID'])
        if CohortVCF is not None:
            hm_source += '+{}'.format(CohortVCF)

        if other_allele is None:
            if returnOtherAllele is True:
                other_allele, hm_TF, hm_vid, hm_code = v_records.infer_OtherAllele(eff=v['effect_allele'],
                                                                                          oa_ensembl=hm_inferOtherAllele)
                hm_matchesVCF, hm_isPalindromic, hm_isFlipped = hm_TF
            else:
                hm_TF, hm_vid = v_records.check_alleles(eff=v['effect_allele'])
                hm_matchesVCF, hm_isPalindromic, hm_isFlipped = hm_TF
        else:
            hm_TF, hm_vid = v_records.check_alleles(eff=v['effect_allele'],
                                                    oa=other_allele)
            hm_matchesVCF, hm_isPalindromic, hm_isFlipped = hm_TF

    if hm_code is None:
        if other_allele is not None:
            hm_code = DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped,
                                                 alleles=[v['effect_allele'], other_allele])
        else:
            hm_code = DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped,
                                                 alleles=[v['effect_allele']])

    # ToDo handle INDEL lookups in VCFs (e.g. ENSEMBL) better
    # ToDo (use allele frequency to resolve ambiguous variants hm_code=3)

    #Output
    if returnOtherAllele is True:
        return pd.Series([hm_source, hm_vid, hm_code, other_allele])
    else:
        return pd.Series([hm_source, hm_vid, hm_code])


def run_HmVCF(args):
    from pgs_harmonizer.variantlookup_tools import VCFs

    ## Set I/O File focations
    # Scoring file location
    if 'loc_scorefiles' in args:
        if not args.loc_scorefiles.endswith('/'):
            args.loc_scorefiles += '/'
        loc_scorefile = args.loc_scorefiles + '{}_hmPOS{}.txt.gz'.format(args.pgs_id, args.target_build)
    else:
        loc_scorefile = 'PGS_HmPOS/{}_hmPOS{}.txt.gz'.format(args.pgs_id, args.target_build)
    try:
        print('Reading Score File')
        header, df_scoring = read_scorefile(loc_scorefile)
    except:
        print('There was an error opening the file!')
        raise IOError
    tqdm.pandas()

    # Define output location
    ofolder = args.loc_outputs
    if ofolder.endswith('/'):
        ofolder = ofolder[:-1]
    if os.path.isdir(ofolder) is False:
        os.mkdir(ofolder)

    # Load Variant References (VCF & Cohort)
    print('Load Variant References (VCF & Cohort)')
    usingCohortVCF = None
    if args.cohort_name is not None:
        vcfs_targetbuild = VCFs(build=args.target_build, cohort_name=args.cohort_name, loc_vcfref=args.loc_vcfref)
        usingCohortVCF = args.cohort_name
        loc_hm_out = '{}/{}_hmVCF{}_{}.txt'.format(ofolder, args.pgs_id, args.target_build, usingCohortVCF)
        args.addOtherAllele = True
    else:
        vcfs_targetbuild = VCFs(build=args.target_build, loc_vcfref=args.loc_vcfref)  # ENSEMBL VCF
        loc_hm_out = '{}/{}_hmVCF{}.txt'.format(ofolder, args.pgs_id, args.target_build)
    if (vcfs_targetbuild.VCF is None) and (len(vcfs_targetbuild.by_chr) == 0):
        print('ERROR: Could not find the VCF')
        raise IOError

    # Start Output
    if args.gzip is True:
        loc_hm_out += '.gz'
        hm_out = gzip.open(loc_hm_out, 'wt')
    else:
        hm_out = open(loc_hm_out, 'w')

    header['HmVCF_date'] = str(datetime.date(datetime.now()))
    if usingCohortVCF is not None:
        header['HmVCF_ref'] = usingCohortVCF
    else:
        header['HmVCF_ref'] = 'Ensembl Variation / dbSNP'  # ToDo Consider adding information about the ENSEMBL build?
    hm_out.write('\n'.join(create_scoringfileheader(header)) + '\n')

    allcols = list(df_scoring.columns)
    if args.addOtherAllele is True:
        allcols += ['hm_source', 'hm_vid', 'hm_code', 'other_allele']
    else:
        allcols += ['hm_source', 'hm_vid', 'hm_code']
    hm_formatter = Harmonizer(df_scoring.columns, returnVariantID=args.addVariantID)
    chrcount = 0
    hm_counts = {}
    hm_Passed = True
    while hm_Passed is True:
        try:
            df_scoring['hm_chr'].fillna('', inplace=True)
            for hm_chr, df_chrom in df_scoring.groupby('hm_chr'):
                if hm_chr == '':
                    print('Harmonizing Chromosome: No HM_CHR')
                else:
                    print('Harmonizing Chromosome: {}'.format(hm_chr))
                df_chrom = df_chrom.copy()
                if args.addOtherAllele is True:
                    df_chrom[['hm_source', 'hm_vid', 'hm_code', 'other_allele']] = df_chrom.progress_apply(variant_HmVCF,
                                                                                                           axis=1,
                                                                                                           vcfs_targetbuild=vcfs_targetbuild,
                                                                                                           CohortVCF=usingCohortVCF,
                                                                                                           returnOtherAllele=True)
                else:
                    df_chrom[['hm_source', 'hm_vid', 'hm_code']] = df_chrom.progress_apply(variant_HmVCF, axis=1,
                                                                                           vcfs_targetbuild=vcfs_targetbuild,
                                                                                           CohortVCF=usingCohortVCF,
                                                                                           returnOtherAllele=False)
                # Post-Harmonization Fixes
                df_chrom.loc[df_chrom['hm_code'].isnull() == False, 'hm_code'] = [conv2int(x) for x in df_chrom.loc[df_chrom['hm_code'].isnull() == False, 'hm_code']]
                if args.skip_strandflips is False:
                    df_chrom = FixStrandFlips(df_chrom)  # also returns new column 'hm_fixedStrandFlip'

                # ToDo unmappable2authorreported
                # if args.author_reported is True:
                #     df_chrom = unmappable2authorreported(df_chrom)

                df_chrom = df_chrom.apply(hm_formatter.format_line, axis=1,
                                          original_build=header['genome_build'])
                df_chrom.columns = hm_formatter.cols_order

                for hm_code, hm_count in dict(df_chrom['hm_code'].value_counts()).items():
                    if hm_code in hm_counts:
                        hm_counts[hm_code] += hm_count
                    else:
                        hm_counts[hm_code] = hm_count

                if chrcount == 0:
                    df_chrom.to_csv(hm_out, mode='a', index=False, sep='\t')  # Write output using pandas
                else:
                    df_chrom.to_csv(hm_out, mode='a', index=False, header=False, sep='\t')  # Write output using pandas
                chrcount += 1
            hm_Passed = 'COMPLETED'
        except:
            hm_Passed = False

    if hm_Passed == 'COMPLETED':
        hm_out.close()
        print('Harmonized {} -> {}'.format(hm_counts, loc_hm_out))
        return
    else:
        hm_out.close()
        os.remove(loc_hm_out)
        print('FAILED')
        raise HarmonizationError
        return


if __name__ == "__main__":
    if args.HmAction == 'HmPOS':
        run_HmPOS(args)
    elif args.HmAction == 'HmVCF':
        run_HmVCF(args)
    else:
        print('Not a valid method, try running: `python Harmonize.py -h` for valid options and more details')