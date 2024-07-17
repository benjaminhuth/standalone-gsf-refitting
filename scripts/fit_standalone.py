from pathlib import Path
import sys

import acts
import acts.examples

from utils import run_fitting


try:
    inputDir = Path(sys.argv[1])
except:
    print(f"Usage: {sys.argv[0]} <path-to-root-files>")

particles = inputDir / "particles.root"
assert particles.exists()

hits = inputDir / "hits.root"
assert hits.exists()

gsfOptions = {
    "betheHeitlerApprox": acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
    "maxComponents": 12,
    "componentMergeMethod": acts.examples.ComponentMergeMethod.maxWeight,
    "mixtureReductionAlgorithm": acts.examples.MixtureReductionAlgorithm.KLDistance,
    "weightCutoff": 1.0e-4,
    "level": acts.logging.INFO,
}

outputDir=Path.cwd()

run_fitting(
    gsfOptions=gsfOptions,
    outputDir=outputDir,
    inputParticles=particles,
    inputHits=hits,
)
