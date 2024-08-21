from pathlib import Path
import os

import numpy as np

import acts
import acts.examples

from acts.examples.simulation import *
from acts.examples.reconstruction import *
from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants

acts_root = Path(os.environ["ACTS_ROOT"])
assert acts_root.exists()

def setup():
    mdecorator = acts.examples.RootMaterialDecorator(
        fileName=str(acts_root / "thirdparty/OpenDataDetector/data/odd-material-maps.root"),
        level=acts.logging.WARNING,
    )
    detector, trackingGeometry, decorators = getOpenDataDetector(
        mdecorator, logLevel=acts.logging.WARNING
    )

    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    return detector, trackingGeometry, field, rnd


def run_fitting(
    gsfOptions,
    outputDir,
    inputParticles,
    inputHits,
    n_events=None,
    n_jobs=-1,
):
    detector, trackingGeometry, field, rnd = setup()
    defaultLogLevel=acts.logging.INFO

    s = acts.examples.Sequencer(
        events=n_events,
        numThreads=n_jobs,
        outputDir=outputDir,
        trackFpes=False,
        logLevel=defaultLogLevel,
    )

    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.FATAL, # Investigate errors
            outputParticles="particles_initial",
            filePath=inputParticles,
        )
    )

    s.addReader(
        acts.examples.RootSimHitReader(
            level=defaultLogLevel,
            filePath=inputHits,
            outputSimHits="simhits",
        )
    )

    #digiConfigFile = acts_root / "Examples/Algorithms/Digitization/share/odd-digi-geometric-config.json"
    #digiConfigFile = acts_root / "Examples/Algorithms/Digitization/share/odd-digi-smearing-config-notime.json"
    digiConfigFile = acts_root / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config-notime.json"
    assert os.path.exists(digiConfigFile)
    addDigitization(
        s=s,
        trackingGeometry=trackingGeometry,
        field=field,
        rnd=rnd,
        logLevel=defaultLogLevel,
        digiConfigFile=digiConfigFile,
    )

    addSeedingTruthSelection(
        s,
        "particles_initial",
        "particles_initial_selected",
        truthSeedRanges=TruthSeedRanges(rho=(0, 2 * u.mm), nHits=(3, None)),
        logLevel=defaultLogLevel,
    )

    seedingSel = acts_root / "Examples/Algorithms/TrackFinding/share/odd-seeding-config.json"
    assert os.path.exists(seedingSel)
    spacePoints = addSpacePointsMaking(
        s,
        trackingGeometry,
        seedingSel,
        logLevel=defaultLogLevel,
    )

    s.addAlgorithm(
        acts.examples.TruthSeedingAlgorithm(
            level=defaultLogLevel,
            inputParticles="particles_initial_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputSpacePoints=[spacePoints],
            outputParticles="truth_seeded_particles",
            outputProtoTracks="truth_seeded_prototracks",
            outputSeeds="seeds",
        )
    )

    s.addAlgorithm(
        acts.examples.TruthTrackFinder(
            level=defaultLogLevel,
            inputParticles="truth_seeded_particles",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputProtoTracks="prototracks",
        )
    )

    s.addAlgorithm(
        acts.examples.TrackParamsEstimationAlgorithm(
            level=defaultLogLevel,
            inputSeeds="seeds",
            inputProtoTracks="prototracks",
            outputTrackParameters="track_parameters",
            outputProtoTracks="prototracks_with_params",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            initialVarInflation=[100.0] * 6,
            particleHypothesis=acts.ParticleHypothesis.electron,
        )
    )

    s.addAlgorithm(
        acts.examples.TrackFittingAlgorithm(
            level=defaultLogLevel,
            inputMeasurements="measurements",
            inputSourceLinks="sourcelinks",
            inputProtoTracks="prototracks_with_params",
            inputInitialTrackParameters="track_parameters",
            outputTracks="kalman_tracks",
            calibrator=acts.examples.makePassThroughCalibrator(),
            fit=acts.examples.makeKalmanFitterFunction(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                magneticField=field,
                energyLoss=True,
                multipleScattering=True,
                reverseFilteringMomThreshold=0.0,
                freeToBoundCorrection=acts.examples.FreeToBoundCorrection(False),
            )
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=defaultLogLevel,
            inputTracks="kalman_tracks",
            inputParticles="truth_seeded_particles",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="kalman_matching",
            outputParticleTrackMatching="kalman_matching2",
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=defaultLogLevel,
            inputTracks="kalman_tracks",
            inputParticles="truth_seeded_particles",
            inputTrackParticleMatching="kalman_matching",
            filePath=str(outputDir / "tracksummary_kalman.root"),
            writeGsfSpecific=False
        )
    )

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            level=defaultLogLevel,
            inputTracks="kalman_tracks",
            outputTracks="gsf_tracks",
            fit=acts.examples.makeGsfFitterFunction(trackingGeometry, field, **gsfOptions)
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=defaultLogLevel,
            inputTracks="gsf_tracks",
            inputParticles="truth_seeded_particles",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="gsf_matching",
            outputParticleTrackMatching="gsf_matching2",
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=defaultLogLevel,
            inputTracks="gsf_tracks",
            inputParticles="truth_seeded_particles",
            inputTrackParticleMatching="gsf_matching",
            filePath=str(outputDir / "tracksummary_gsf_refit.root"),
            writeGsfSpecific=True
        )
    )

    s.run()


