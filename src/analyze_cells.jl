"""
Do simple per-cell analysis.
"""
module AnalyzeCells

export compute_vector_of_total_UMIs_per_cell!

using DataAxesFormats
using StatsBase
using TanayLabUtilities

using ..Contracts

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_total_UMIs_per_cell

"""
    function compute_vector_of_total_UMIs_per_cell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_total_UMIs_per_cell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(CreatedOutput),
    ],
) function compute_vector_of_total_UMIs_per_cell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    total_UMIs_per_cell = daf["@ cell @ gene [ ! is_excluded ] :: UMIs >| Sum"].array
    set_vector!(daf, "cell", "total_UMIs", total_UMIs_per_cell; overwrite)
    @debug "Mean (included) UMIs per cell: $(mean(total_UMIs_per_cell))" _group = :mcs_details  # NOLINT
    return nothing
end

end  # module
