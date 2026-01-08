"""
Do simple per-cell analysis.
"""
module AnalyzeCells

export compute_cells_total_UMIs!

using DataAxesFormats
using TanayLabUtilities

using ..Contracts

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_total_UMIs_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_excluded_vector

"""
    function compute_cells_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of UMIs of all the non-excluded genes in each cell.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(GuaranteedOutput),
    ],
) function compute_cells_total_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    total_UMIs_per_cell = daf["/ gene &! is_excluded / cell : UMIs %> Sum"].array
    set_vector!(daf, "cell", "total_UMIs", total_UMIs_per_cell; overwrite)
    return nothing
end

end  # module
