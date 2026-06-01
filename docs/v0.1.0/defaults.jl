"""
Default values for parameters.
"""
module Defaults

export GENE_FRACTION_REGULARIZATION_FOR_CELLS
export GENE_FRACTION_REGULARIZATION_FOR_METACELLS
export MIN_SIGNIFICANT_GENE_UMIS

"""
When computing log (base 2) of the fraction of the gene expression out of the total of "many" cells (such as all the
cells in a metacell), we use this regularization factor by default to avoid zero values. This is lower than
[`GENE_FRACTION_REGULARIZATION_FOR_CELLS`](@ref).
"""
GENE_FRACTION_REGULARIZATION_FOR_METACELLS = 1e-5

"""
When computing log (base 2) of the fraction of the gene expression out of the total of a single cell, we use this
regularization factor by default to avoid zero values. This is higher than
[`GENE_FRACTION_REGULARIZATION_FOR_METACELLS`](@ref).
"""
GENE_FRACTION_REGULARIZATION_FOR_CELLS = 1e-4

"""
When comparing gene expression levels (e.g., when computing fold factors), we do not consider the result to be
significant unless at least this number of UMIs were used to compute the compared value(s).
"""
MIN_SIGNIFICANT_GENE_UMIS = 40

end
