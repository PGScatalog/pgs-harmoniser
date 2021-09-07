import os
import pandas as pd
from cyvcf2 import VCF
from pgs_harmonizer.harmonize import reversecomplement, chromosomes, DetermineHarmonizationCode


class VCFResult:
    """Class to parse the results of a VCF Lookup and compare alleles"""

    def __init__(self, chr, pos, build, vcf_result=None):
        self.vcf_query = [chr, pos, build]
        self.vcf_result = vcf_result # List of variants
        if self.vcf_result is None:
            if chr in chromosomes:
                self.vcf_result = vcf_lookup(*self.vcf_query)

    def check_alleles(self, eff, oa=None):
        """Check if this variant exists in the ENSEMBL VCF files"""
        if (self.vcf_result is None) or (len(self.vcf_result) == 0):
            return (False, False, False), None
        else:
            # Collect variants that match or mismatch the VCF
            v_hm = []
            v_hm_code = []

            # Loop through all overlapping variants
            for v in self.vcf_result:
                hm_matchesVCF = False
                hm_isPalindromic = False
                hm_isFlipped = False

                alleles = [v.REF] + v.ALT
                alleles_rc = [reversecomplement(x) for x in alleles]
                # print('Alleles: {} | RC: {}'.format(alleles, alleles_rc))

                # Check allele(s) against VCF
                if oa is not None:
                    if (eff in alleles) and (oa in alleles):
                        hm_matchesVCF = True
                        if (eff in alleles_rc) and (oa in alleles_rc):
                            hm_isPalindromic = True
                    elif (eff in alleles_rc) and (oa in alleles_rc):
                        hm_matchesVCF = True
                        hm_isFlipped = True
                else:
                    if eff in alleles:
                        hm_matchesVCF = True
                        if eff in alleles_rc:
                            hm_isPalindromic = True
                    elif eff in alleles_rc:
                        hm_matchesVCF = True
                        hm_isFlipped = True

                hm_tuple = (hm_matchesVCF, hm_isPalindromic, hm_isFlipped)
                v_hm.append([hm_tuple, v.ID])
                v_hm_code.append(DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped))
        return v_hm[v_hm_code.index(max(v_hm_code))]

    def infer_OtherAllele(self, eff, oa_ensembl=None, allowINDELs = False):
        """Try to infer the reference_allele from a vcf lookup (assumes the vcf is sorted by variant quality)"""
        if oa_ensembl is not None:
            oa_ensembl = oa_ensembl.split('/')  # [list of potential alleles, e.g. from ENSEMBL]
        oa_mapped = []  # list of potential reference alleles (if there are more variants)
        oa_mapped_hmcodes = []  # list of potential reference alleles (if there are more variants)

        # Loop through cohort variants to look for acceptable other_alleles
        for v in self.vcf_result:
            alleles = [v.REF] + v.ALT
            if len(alleles) < 2:
                continue
            alleles_rc = [reversecomplement(x) for x in alleles]

            minorallele = v.INFO.get('MA')
            if minorallele is not None:
                minorallele_rc = reversecomplement(v.INFO.get('MA'))
            else:
                minorallele_rc = None

            hm_matchesVCF = False
            hm_isPalindromic = False
            hm_isFlipped = False
            oa = None

            # Assume effect_allele is on positive strand
            if eff in v.ALT:
                oa = v.REF
            elif eff == v.REF:
                if len(v.ALT) == 1:
                    oa = v.ALT[0]
                else:
                    if (minorallele is not None) and (minorallele != eff):
                        oa = minorallele
                    else:
                        oa = v.ALT[0] # This is just choosing the first one

            if oa is not None:
                hm_matchesVCF = True
                if eff in alleles_rc:
                    hm_isPalindromic = True
            else:
                if eff in alleles_rc:
                    hm_matchesVCF = True
                    hm_isFlipped = True
                    alleles_rc.remove(eff)
                    if (minorallele_rc is not None) and (minorallele_rc in alleles_rc):
                        oa = minorallele_rc
                    else:
                        oa = alleles_rc[0] # ToDo Make more informed choice about the other_allele (should match ENSEMBL, MAF)

            hm = (hm_matchesVCF, hm_isPalindromic, hm_isFlipped)

            # If we know the potential alleles (e.g. from ENSEMBL API/DB) check them
            # before adding this variant/allele combination as an acceptable match (should
            # get rid of weird INDELs that overlap without using rsID that isn't present in each VCF)
            if oa_ensembl is not None:
                oa_append = False
                if oa in oa_ensembl:
                    oa_append = True
                elif (hm_isFlipped is True) and (oa in [reversecomplement(x) for x in oa_ensembl]):
                    oa_append = True
            else:
                oa_append = True

            # allowINDELs (e.g. variant alleles with length > 1)
            if (oa_append is True) and (allowINDELs is False) and (oa is not None):
                if len(oa) > 1:
                    oa_append = False

            if oa_append is True:
                oa_mapped.append([oa, hm, v.ID])

                if hm == (True, False, False):
                    oa_mapped_hmcodes.append(5)
                elif hm == (True, True, False):
                    oa_mapped_hmcodes.append(4)
                elif hm == (True, False, True):
                    oa_mapped_hmcodes.append(-4)
                else:
                    oa_mapped_hmcodes.append(-5)

        # Select the best OA
        #print(list(zip(oa_mapped, oa_mapped_hmcodes)))
        if len(oa_mapped_hmcodes) >= 1:
            i_bestmap = oa_mapped_hmcodes.index(max(oa_mapped_hmcodes))
            return oa_mapped[i_bestmap][0], oa_mapped[i_bestmap][1], oa_mapped[i_bestmap][2], oa_mapped_hmcodes[i_bestmap]
        else:
            return None, (False, False, False), None, -5


def vcf_lookup(chromosome, position, build, loc_vcfref='map/vcf_ref/'):
    """Retrieve variants overlapping with a position in a VCF File"""
    if build not in ['GRCh37', 'GRCh38']:
        raise ValueError("Invalid Genome Build. Expected one of: GRCh37, GRCh38")

    if chromosome not in chromosomes:
        raise ValueError("Invalid Chromosome. Expected one of: {}".format(chromosomes))

    loc_vcf = f'{loc_vcfref}{build}/homo_sapiens-chr{chromosome}.vcf.gz'
    # If the given path hasn't a subdirectory for the assembly
    if not os.path.isfile(loc_vcf):
        loc_vcf = f'{loc_vcfref}/homo_sapiens-chr{chromosome}.vcf.gz'
    vcf = VCF(loc_vcf)
    if (type(position) is str) and ('-' in position):
        return list(vcf('{}:{}'.format(chromosome, position)))
    else:
        return list(vcf('{}:{}-{}'.format(chromosome, position, position)))


class VCFs:
    """Class to open and hold all VCF files for a genome build"""
    def __init__(self, build, loc_vcfref='map/vcf_ref/', cohort_name=None):
        self.VCF = None
        self.by_chr = {}
        self.build = build
        if cohort_name is None:
            for chr in chromosomes:
                loc_vcf = f'{loc_vcfref}{self.build}/homo_sapiens-chr{chr}.vcf.gz'
                # If the given path hasn't a subdirectory for the assembly
                if not os.path.isfile(loc_vcf):
                    loc_vcf = f'{loc_vcfref}/homo_sapiens-chr{chr}.vcf.gz'
                self.by_chr[chr] = VCF(loc_vcf)
        else:
            loc_vcf = '{}/{}/cohort_ref/{}.vcf.gz'.format(loc_vcfref, build, cohort_name)
            if os.path.isfile(loc_vcf):
                self.VCF = VCF(loc_vcf)

    def vcf_lookup(self,chromosome, position, rsid = None):
        """Lookup a variant in a specific genome build"""
        if chromosome not in chromosomes:
            return VCFResult(chromosome, position, self.build, [])

        if (type(position) is str) and ('-' in position):
            query = '{}:{}'.format(chromosome, position)
        else:
            query = '{}:{}-{}'.format(chromosome, position, position)

        if self.VCF is None:
            r_lookup = list(self.by_chr[chromosome](query))
        else:
            r_lookup = list(self.VCF(query))

        # Filter results by rsID (if supplied)
        if rsid is not None:
            v_rsid = [v for v in r_lookup if v.ID == rsid] # List of variants that match rsID
            if len(v_rsid) > 0:
                r_lookup = v_rsid # Replace with list of variant(s) that match rsID, else supplies all variants at locus

        return VCFResult(chromosome, position, self.build, r_lookup)


# def guess_build(loc_file, vcf_root):
