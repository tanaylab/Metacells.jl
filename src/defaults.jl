"""
Default values for parameters.
"""
module Defaults

export GENE_FRACTION_REGULARIZATION
export MIN_SIGNIFICANT_GENE_UMIS

"""
When computing log (base 2) of the fraction of the gene expression out of the total, we use
this regularization factor by default top avoid zero values.
"""
GENE_FRACTION_REGULARIZATION = 1e-5

"""
When comparing gene expression levels (e.g., when computing fold factors), we do not consider the result to be
significant unless at least this number of UMIs were used to compute the compared value(s).
"""
MIN_SIGNIFICANT_GENE_UMIS = 40

end
