from pathlib import Path

import acts
import acts.examples

from utils import run_fitting

bha = acts.examples.AtlasBetheHeitlerApprox.loadFromFiles(
    lowParametersPath="../config/GeantSim_LT01_cdf_nC6_O5.par",
    highParametersPath="../config/GeantSim_GT01_cdf_nC6_O5.par",
    lowLimit=0.1,
    highLimit=0.3,
)

gsfOptions = {
    "betheHeitlerApprox": bha,
    "maxComponents": 12,
    "componentMergeMethod": acts.examples.ComponentMergeMethod.maxWeight,
    "mixtureReductionAlgorithm": acts.examples.MixtureReductionAlgorithm.KLDistance,
    "weightCutoff": 1.0e-4,
    "level": acts.logging.INFO,
}

outputDir=Path(snakemake.output[0]).parent

run_fitting(
    gsfOptions=gsfOptions,
    outputDir=outputDir,
    inputParticles=snakemake.input[0],
    inputHits=snakemake.input[1]
)
