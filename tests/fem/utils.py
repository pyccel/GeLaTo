# coding: utf-8

from numpy import allclose

def assert_identical_coo(A, B):

    # ...
    assert(A.shape == B.shape)
    assert(A.nnz == B.nnz)

    assert(allclose(A.row,  B.row))
    assert(allclose(A.col,  B.col))
    assert(allclose(A.data, B.data))

#    assert(allclose(A.todense(), B.todense()))
    # ...
