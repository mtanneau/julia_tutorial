using DelimitedFiles
using SparseArrays
using LDLFactorizations

# warmup
A = [ 1.7     0     0     0     0     0     0     0   .13     0
        0    1.     0     0   .02     0     0     0     0   .01
        0     0   1.5     0     0     0     0     0     0     0
        0     0     0   1.1     0     0     0     0     0     0
        0   .02     0     0   2.6     0   .16   .09   .52   .53
        0     0     0     0     0   1.2     0     0     0     0
        0     0     0     0   .16     0   1.3     0     0   .56
        0     0     0     0   .09     0     0   1.6   .11     0
      .13     0     0     0   .52     0     0   .11   1.4     0
        0   .01     0     0   .53     0   .56     0     0   3.1 ];

Ti = Int
A = convert(SparseMatrixCSC{T,Ti}, A)  # only upper triangle will be accessed
b = T[.287, .22, .45, .44, 2.486, .72, 1.55, 1.424, 1.621, 3.759]

LDLT = ldl(A)
x = LDLT \ b

using MatrixMarket
A = MatrixMarket.mmread("K_0.mtx")
rhs = readdlm("rhs_0.rhs")[:]

ldl(A)
@time ldl(A)  # time in double precision

A = convert(SparseMatrixCSC{T,Ti}, A)
rhs = convert(Vector{T}, rhs)

@time ldl(A)  # time in extended precision
