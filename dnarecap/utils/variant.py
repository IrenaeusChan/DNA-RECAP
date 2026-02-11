class Variant(object):
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = str(chrom)
        self.pos = int(pos)
        self.ref = str(ref)
        self.alt = str(alt)
        if self.is_snp(): self.mut_type = "snp"
        if self.is_insertion(): self.mut_type = "insertion"
        if self.is_deletion(): self.mut_type = "deletion"
        if self.is_complex(): self.mut_type = "complex"

    @classmethod
    def createVariant(self, variant): return Variant(variant.split(":")[0], variant.split(":")[1], variant.split(":")[2], variant.split(":")[3])

    @classmethod
    def createVariantManual(self, chrom, pos, ref, alt): return Variant(chrom, pos, ref, alt)

    # Single base mutation e.g. A>T
    def is_snp(self):
        if len(self.ref) == 1 and len(self.alt) == 1: return True
        else: return False

    # Multiple base mutation e.g. ACGT>ACGTT
    def is_indel(self):
        if len(self.ref) | len(self.alt) > 1: return True
        else: return False

    # There are two types of INDELs: Insertions and Deletions
    # Insertions are when the reference has a single base and the ALT has multiple bases
    def is_insertion(self):
        if (len(self.ref) == 1 and len(self.alt) > 1 and self.ref[0] == self.alt[0]): return True
        else: return False

    # Deletions are when the reference has multiple bases and the ALT has a single base
    def is_deletion(self):
        if (len(self.ref) > 1 and len(self.alt) == 1 and self.ref[0] == self.alt[0]): return True
        else: return False

    # Complex is when the reference has multiple bases and the ALT has multiple bases
    def is_complex(self):
        if (len(self.ref) > 1 and len(self.alt) > 1): return True
        elif (self.is_indel() and self.ref[0] != self.alt[0]): return True
        else: return False

    def toString(self): return(self.chrom + ":" + str(self.pos) + ":" + self.ref + ":" + self.alt)