nscancor
========

An R package for non-negative and sparse canonical correlation
analysis (CCA).

CCA is a method for finding associations between paired data sets.
For example, a health study might record the gene expression levels and a
number of physiological parameters for a patient cohort. If one
conjectures that the cause for the physiological symptoms has a
genetic component, one could expect to find a correlation between the
expression of certain genes and the strength of certain symptoms. CCA
finds a pair of linear projections (called _canonical vectors_), one
for each data modality, such that the projected values (called
_canonical variables_) have maximum correlation. The next pair of
canonical variables is found by again maximizing their correlation,
under the additional constraint that the they have to be uncorrelated
to all previous ones, and so on.

CCA was first introduced by Hotelling in 1936, and has many
similarities to principal component analysis (PCA). Where the PCA
solution is computed from the eigenvalue decomposition (EVD) of the
covariance matrix of a single data set, the CCA solution is computed
from the EVD of the cross-covariance matrix of the two data sets. This
approach is very efficient, but one sometimes encounters the following
problems during an analysis. First, if at least one of the data sets
contains more features than samples (a common case for gene expression
data), there exist an infinite number of trivial projections that
achieve perfect correlation. Regularization of the canonical vectors
is necessary to again solve a well-posed problem.  Second, the
projections are typically linear combinations with non-zero weights
for all features, which makes an interpretation of the weights
difficult. A sparse solution which only includes a small number of
important features is often desirable.

This package implements a CCA algorithm called `nscancor` which can
enforce appropriate constraints on the canonical vectors to address
both aforementioned problems. Enforcing a bound on the Euclidean norm
(also called the L2 norm) of the projections avoids trivial
correlations. Enforcing a bound on the L1 norm leads to sparse
solutions, where many of the weights are exactly zero. And enforcing
non-negativity of the projection weights is useful for analysing data
where only positive influence of features is deemed appropriate. The
algorithm executes iterated regression steps, and the constraints
enter via the regression functions. `nscancor` is therefore modular,
and builds on the many regression methods that are
available, e.g. ridge regression or the elastic net. By using two
different regression functions, the proper constraints can be enforced
for each domain.

The package also provides a generalization of constrained CCA for
analyzing more than two data sets. The `mcancor` algorithm is
structurally analogous to `nscancor`, but it maximizes the sum of all
pairwise correlations of canonical variables. As with `nscancor`,
specifying the regression function for each domain makes it possible
to enforce appropriate constraints on each canonical vector.

[This blog
post](http://sigg-iten.ch/learningbits/2014/01/20/canonical-correlation-analysis-under-constraints/)
explains how to use the package and demonstrates its benefits.
