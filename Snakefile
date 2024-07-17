envvars:
    "ACTS_ROOT",


rule simulate_electrons:
    output:
        "{pt}GeV/particles.root",
        "{pt}GeV/hits.root",
    params:
        n_events=10,
        p_min=lambda wildcards: float(wildcards.pt),
        p_max=lambda wildcards: float(wildcards.pt),
        abs_eta=3,
        p_transverse=True,
        uniform_eta=True,
        pdg=11,
        randomize_charge=True,
    script:
        "scripts/simulation.py"

rule refit:
    input:
        "{pt}GeV/particles.root",
        "{pt}GeV/hits.root",
    output:
        "{pt}GeV/tracksummary_kalman.root",
        "{pt}GeV/tracksummary_gsf_refit.root",
    script:
        "scripts/fit.py"


rule all:
    default_target: True
    input:
        "10GeV/tracksummary_kalman.root",
        "10GeV/tracksummary_gsf_refit.root",
