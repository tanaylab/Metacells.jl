"""
Do simple per-cell analysis.
"""
module AnalyzeCells

export compute_cells_total_UMIs!
export compute_cells_covered_UMIs!

using DataAxesFormats
using TanayLabUtilities

using ..Contracts

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_covered_UMIs_vector
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_total_UMIs_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_excluded_vector

"""
    function compute_cells_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute the total UMIs of the genes in each cell.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [cell_axis(OptionalInput), gene_axis(RequiredInput)],
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

"""
    function compute_cells_covered_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute the total UMIs of the covered genes in each cell.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [cell_axis(OptionalInput), gene_axis(RequiredInput)],
    data = [
        cell_gene_UMIs_matrix(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        cell_covered_UMIs_vector(GuaranteedOutput),
    ],
) function compute_cells_covered_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    covered_UMIs_per_cell = daf["/ gene & is_covered / cell : UMIs %> Sum"].array
    set_vector!(daf, "cell", "covered_UMIs", covered_UMIs_per_cell; overwrite)
    return nothing
end

end  # module
