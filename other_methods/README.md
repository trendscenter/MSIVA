# Alternative methods
These are alternative methods used for comparison throughout the examples.

## `icatb_runica.m`
This is the standard Infomax ICA algorithm provided with the [GIFT] toolbox.

[GIFT]: http://trendscenter.org/software/gift/

## `icatb_iva_second_order.m`
This is the standard IVA-G algorithm provided with the [GIFT] toolbox.

## `icatb_iva_laplace.m`
This is the standard IVA-L algorithm provided with the [GIFT] toolbox.

## `icatb_iva_laplace_bkt.m`
This is a modification of the standard IVA-L algorithm provided with the [GIFT] toolbox. This modification implements a sufficient decrease check for step length as well as backtracking and rules for setting the step length.

## `jbd.m`
This is an implementation of [JBD-SOS] for ISA based on second-order statistics via joint block diagonalization of covariance matrices.

[JBD-SOS]: https://www.irit.fr/~Dana.Lahat/jbd.zip

## `isa_est.m`
This is an implementation of [ISA] based on uncorrelated multivariate Laplace higher-order statistics (uses no second-order satatistics).

[ISA]: http://ai.stanford.edu/~quocle/video_release.tar.gz
