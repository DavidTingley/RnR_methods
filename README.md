# RnR_methods
> A toolbox for detecting/quantifying/comparing reactivation and replay in neural data.
Some of these methods and tests were described in <a href="" target="_blank">Tingley & Peyrache (2019)</a>

- [Installation](#installation)
- [Reactivation](#Reactivation)
- [Replay](#replay)
- [Test and compare methods](#compare)
- [License](#license)

## Installation
- Clone this repo to your local machine using `https://github.com/davidtingley/RnR_metthods`
- Run `addpath(genpath(pwd))` from the \RnR_methods directory

## Reactivation
- 'ReactStrength', to compute reactivation (i.e. 0-lag neuronal correlation) with PCA or ICA.
Type 'help ReactStrength' for mor info

## Repaly

## Compare
- 'compareReplayMethods', to compare different reactivation methods (Fig. 3 in <a href="" target="_blank">Tingley&Peyrache</a>)
- 'compareNoiseDegradation', to compare how noise in data affects replay (and reactivation) (Fig. 4 in <a href="" target="_blank">Tingley&Peyrache</a>)
- 'rankOrder_falsePositive_test', to see how limited length sequences and within even shuffling can lead to high FP rates/
- 'compareBinSize', to compare how bin size affects replay detection

## License
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/DavidTingley/RnR_methods/master)
