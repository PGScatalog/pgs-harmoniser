import pandas as pd
import numpy as np
import re
import json
import gzip

remap_header = {
    'PGS ID': 'pgs_id',
    'PGS Name': 'pgs_name',
    'Reported Trait': 'trait_reported',
    'Original Genome Build': 'genome_build',
    'Number of Variants': 'variants_number',
    'PGP ID': 'pgp_id',
    'Citation': 'citation',
    'LICENSE': 'pgs_license',
    # Harmonization related
    'HmPOS Build': 'HmPOS_build',
    'HmPOS Date':'HmPOS_date',
    'HmVCF Reference': 'HmVCF_ref',
    'HmVCF Date': 'HmVCF_date'
} # ToDo remove once Scoring File headers are fixed


chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']
acceptable_alleles = re.compile('[ACGT]*$')  # alleles that can be reverse complemented


def reversecomplement(x):
    if acceptable_alleles.match(x):
        return ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
    else:
        return None


def conv2int(n):
    try:
        return int(n)
    except:
        return n


def read_scorefile(loc_scorefile):
    """Loads PGS Catalog Scoring file and parses the header into a dictionary"""
    if loc_scorefile.endswith('.gz'):
        f = gzip.open(loc_scorefile,'rt')
    else:
        f = open(loc_scorefile, 'rt')

    header = {}
    lastline = '#'
    while lastline.startswith('#'):
        lastline = f.readline()
        line = lastline.strip()
        if line.startswith('#'):
            if '=' in line:
                line = line[1:].split('=')
                field, val = [x.strip() for x in line]
                header[remap_header[field]] = val # ToDo change once Scoring File headers are fixed
    f.close()

    if ('genome_build' in header) and (header['genome_build'] == 'NR'):
        header['genome_build'] = None

    df_scoring = pd.read_table(loc_scorefile, float_precision='round_trip', comment='#',
                               dtype = {'chr_position' : 'object',
                                        'chr_name' : 'object',
                                        'hm_chr': 'object',
                                        'hm_pos': 'object'})

    # Make sure certain columns maintain specific datatypes
    if 'reference_allele' in df_scoring.columns:
        df_scoring = df_scoring.rename(columns={"reference_allele": "other_allele"})

    for poscol in ['chr_position', 'hm_pos']:
        if poscol in df_scoring.columns:
            df_scoring.loc[df_scoring[poscol].isnull() == False, poscol] = [conv2int(x) for x in df_scoring.loc[df_scoring[poscol].isnull() == False, poscol]]

    return(header, df_scoring)


def create_scoringfileheader(h):
    """Function to extract score & publication information for the PGS Catalog Scoring File commented header"""
    # Recreate original header
    lines = [
        '### PGS CATALOG SCORING FILE - see www.pgscatalog.org/downloads/#dl_ftp for additional information',
        '## POLYGENIC SCORE (PGS) INFORMATION',
        '# PGS ID = {}'.format(h['pgs_id'])
    ]
    if 'pgs_name' in h:
        lines.append('# PGS Name = {}'.format(h['pgs_name']))
    lines += ['# Reported Trait = {}'.format(h['trait_reported']),
              '# Original Genome Build = {}'.format(str(h['genome_build']).replace('None', 'NR')),
              '# Number of Variants = {}'.format(h['variants_number']),
              '## SOURCE INFORMATION',
              '# PGP ID = {}'.format(h['pgp_id']),
              '# Citation = {}'.format(h['citation'])
    ]
    if 'pgs_license' in h:
        lines.append('# LICENSE = {}'.format(h['pgs_license']))  # Append to header

    # Add Harmonization Details
    ## HmPOS
    if 'HmPOS_build' in h:
        lines += ['## HARMONIZATION DETAILS',
                  '# HmPOS Build = {}'.format(h['HmPOS_build']),
                  '# HmPOS Date = {}'.format(h['HmPOS_date'])
                  ]
        if 'HmVCF_ref' in h:
            lines += ['# HmVCF Reference = {}'.format(h['HmVCF_ref']),
                      '# HmVCF Date = {}'.format(h['HmVCF_date'])
                      ]
    return lines


def DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped,alleles = [], ):
    hm_coding = {
        (True, False, False): 5,
        (True, True, False): 4,
        (True, False, True): -4,
        (False, False, False): -5
    }
    hm_code = hm_coding.get((hm_matchesVCF, hm_isPalindromic, hm_isFlipped), -1)
    if hm_code == 4:
        if len(alleles) > 1:
            rc_alleles = [reversecomplement(x) for x in alleles]
            for A in alleles:
                if A in rc_alleles:
                    hm_code = 3
    return hm_code


def FixStrandFlips(df):
    df['hm_fixedStrandFlip'] = np.nan

    # Correct the effect_allele
    df['hm_reported_effect_allele'] = np.nan
    df.loc[df['hm_code'] == -4, 'hm_reported_effect_allele'] = df.loc[df['hm_code'] == -4, 'effect_allele']
    df.loc[df['hm_code'] == -4, 'effect_allele'] = df.loc[df['hm_code'] == -4, 'effect_allele'].apply(reversecomplement)

    # Correct the other_allele
    if 'other_allele' in df.columns:
        df['hm_reported_other_allele'] = np.nan
        df.loc[df['hm_code'] == -4, 'hm_reported_other_allele'] = df.loc[df['hm_code'] == -4, 'other_allele']
        df.loc[df['hm_code'] == -4, 'other_allele'] = df.loc[df['hm_code'] == -4, 'other_allele'].apply(reversecomplement)

    df.loc[df['hm_code'] == -4, 'hm_fixedStrandFlip'] = True
    return df


# def unmappable2authorreported(df):
#     # ToDo unmappable2authorreported
#     return df


class Harmonizer:
    """Class to select and harmonize variant locations in a PGS Scoring file."""
    def __init__(self, cols, returnVariantID=False):
        """Used to select the columns and ordering of the output PGS Scoring file"""
        self.cols_previous = cols

        # Decide what to return in variant_id column:
        self.has_rsID = False
        if 'rsID' in cols:
            self.has_rsID = True

        self.variant_id_vcf = returnVariantID  # return the variant_id from the VCF file or rsID as default

        # Standard set of columns (vaguely like VCF)
        self.cols_order = ['chr_name', 'chr_position', 'variant_id',
                           'effect_allele', 'other_allele', 'effect_weight',
                           'hm_code', 'hm_info']
        self.cols_extra = []

        # Check which other columns need to be added
        for c in self.cols_previous:
            if (c not in self.cols_order) and (c.startswith('hm_') is False) and (c != 'rsID'):
                self.cols_order.append(c)
                self.cols_extra.append(c)

    def format_line(self, v, original_build):
        """Method that takes harmonized variant location and compares it with the old information.
        Outputs any changes to an hm_info dictionary"""
        # Initialize Output
        l_output = ['']*len(self.cols_order)
        hm_info = {'hm_source': v['hm_source']}
        missing_val = '.'

        # Decide how to output the variant (pass/fail harmonization)
        PassHM = True
        hm_code = v['hm_code']
        if hm_code is None:
            v['hm_code'] = -1

        if hm_code < 0:  # Proceed as if it's not mapped
            PassHM = False
            if hm_code == -4:
                if 'hm_fixedStrandFlip' in v:
                    if v['hm_fixedStrandFlip'] is True:
                        hm_info['fixedStrandFlip'] = True
                        hm_info['reported_effect_allele'] = v['hm_reported_effect_allele']
                        if 'reported_other_allele' in v:
                            hm_info['reported_other_allele'] = v['hm_reported_other_allele']
                        PassHM = True
                    else:
                        hm_info['fixedStrandFlip'] = False

        if original_build is None:
            original_build = 'NR'

        # Define output
        if PassHM is True:
            # MANDATORY:chr_name|chr_position|effect_allele|effect_weight|hm_code
            l_output[self.cols_order.index('chr_name')] = v['hm_chr']
            l_output[self.cols_order.index('chr_position')] = v['hm_pos']
            l_output[self.cols_order.index('effect_allele')] = v['effect_allele']
            l_output[self.cols_order.index('effect_weight')] = v['effect_weight']
            l_output[self.cols_order.index('hm_code')] = str(hm_code)

            # MANDATORY: other_allele
            if ('other_allele' in v) and (pd.isnull(v['other_allele']) is False):
                l_output[self.cols_order.index('other_allele')] = v['other_allele']
            else:
                l_output[self.cols_order.index('other_allele')] = missing_val

            # MANDATORY: variant_id
            if self.variant_id_vcf is True:
                l_output[self.cols_order.index('variant_id')] = v['hm_vid']  # Return the ID from VCF

            if self.has_rsID:
                rsID = v['hm_rsID']
                old_rsID = v['rsID']
                if pd.isnull(rsID) is False:
                    if rsID != old_rsID:
                        hm_info['previous_rsID'] = old_rsID  # Compare to old rsID and write old-value to hm_info
                    # Decide what to output
                    if self.variant_id_vcf is False:
                        l_output[self.cols_order.index('variant_id')] = rsID
                    else:
                        hm_info['rsID'] = rsID  # Write to hm_info for provenance
            if l_output[self.cols_order.index('variant_id')] == '':
                l_output[self.cols_order.index('variant_id')] = missing_val  # Return missing val

            # MANDATORY: hm_info
            if len(hm_info) > 0:
                l_output[self.cols_order.index('hm_info')] = json.dumps(hm_info)

            # OPTIONAL: Add information from extra columns
            for colname in self.cols_extra:
                l_output[self.cols_order.index(colname)] = v[colname]
        else:
            hm_info['reported_build'] = original_build

            if self.has_rsID is True:
                rsID = v['hm_rsID']
                old_rsID = v['rsID']
                if pd.isnull(rsID) is False:  # mapped by rsID
                    hm_info['rsID'] = rsID
                    if rsID != old_rsID:
                        hm_info['previous_rsID'] = old_rsID
                else:
                    hm_info['rsID'] = old_rsID

            for i, colname in enumerate(self.cols_order):
                if colname == 'chr_name':
                    if pd.isnull(v['hm_chr']) is False:
                        hm_info['hm_chr'] = v['hm_chr']
                elif colname == 'chr_position':
                    if pd.isnull(v['hm_pos']) is False:
                        hm_info['hm_pos'] = v['hm_pos']
                elif colname == 'variant_id':
                    if pd.isnull(v['hm_vid']) is False:
                        hm_info['variant_id'] = v['hm_vid']
                    l_output[i] = missing_val
                elif colname in ['effect_allele', 'other_allele']:
                    val = v[colname]
                    if pd.isnull(val) is False:
                        hm_info[colname] = val
                elif colname == 'hm_info':
                    l_output[i] = json.dumps(hm_info)
                elif colname == 'hm_code':
                    l_output[i] = str(hm_code)
                else:
                    l_output[i] = v[colname]  # Output everything that isn't a mandatory column as a value

        return pd.Series(l_output)




