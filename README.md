# is2

is2 provides a very efficient inference approach of nonlinear partially-observed Markov processes. It is inherited from R package pomp, so it can virtually be applied to every model therein. In is2, three second order iterated smoothing algorithm are provided, expoiting efficiently estimation of the observed information matrix to increase convergence rate.

At the moment, support is provided for

fixed lag smoothing
second order iterated smoothing algorithm of Nguyen and Ionides (2017),
second order iterated smoothing algorithm of Doucet et al. (2013),
reduced second order iterated smoothing algorithm of Nguyen et al. (2015)(submitted).
accelerate iterated filtering
average iterated filtering
momentum iterated filtering
particle iterated filtering

Simple worked examples are provided in the test directory of the installed package.

Future support for a variety of other algorithms is envisioned.

Please let the developers know if you find is2 useful, if you publish results obtained using it, if you come up with improvements, find bugs, or have suggestions or feature requests!

The package is provided under the GPL. Contributions are welcome, as are comments, suggestions, feature requests, examples, and bug reports. Please send these to dxnguyen at olemiss dot edu.

To keep abreast of new developments, subscribe to the is2-announce mailing list.
