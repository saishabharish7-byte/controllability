# Exact controllability of an *E. coli* metabolic network

## Project question
How does excluding high-connectivity currency metabolites change the estimated minimum number of driver metabolites in a central-metabolism network?

## Why this is a different application
The reference paper develops exact controllability for complex networks. This project applies the same mathematical tool to a **metabolic reaction network**, rather than a protein–protein interaction network. Metabolites are nodes; directed substrate-to-product relationships within reactions form edges.

## Dataset included
`data/ecoli_central_metabolism_subset.csv` is a small, hand-curated teaching subset of *E. coli* central metabolism (glycolysis, pentose phosphate pathway, TCA-cycle and anaplerotic reactions). It is intended for a transparent, reproducible first implementation; it is **not** a complete genome-scale reconstruction and should not be described as one.

For a larger, externally sourced dataset, the repository includes `scripts/download_bigg_model.py`, which downloads the public BiGG Models *E. coli* core model SBML file. Run the optional full-model workflow described below. The exact model URL and retrieval date are recorded by the downloader.

## Network construction
For each reaction, each substrate is connected to each product. Reversible reactions also receive reverse edges. The binary adjacency matrix uses `A[i, j] = 1` when metabolite i has an edge to metabolite j`. This is a directed projection of reaction stoichiometry, not a kinetic model.

Two networks are compared:
1. Full participant network.
2. Currency-metabolite-reduced network, excluding common carriers/cofactors (ATP, ADP, AMP, NAD(H), NADP(H), FAD(H2), CoA, water, phosphate, protons and CO2). The list is explicit in `scripts/analyze.py` and can be modified with justification.

## Exact controllability calculation
For adjacency matrix A, the minimum number of driver nodes is estimated as:

`N_D = max_lambda [ N - rank(lambda I - A) ]`

where the expression in brackets is the geometric multiplicity of eigenvalue lambda. The driver-node fraction is `n_D = N_D / N`.

Numerical eigenvalue clustering and SVD rank tolerances are used. These are numerical approximations; report the tolerance and perform the suggested sensitivity check before making strong claims. The network is unweighted in the included implementation.

## Run it (Python 3.10+)
```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate
pip install -r requirements.txt
python scripts/analyze.py
```

Outputs are written to `results/`:
- `controllability_comparison.csv`
- `network_comparison.png`
- `metabolite_edges_full.csv`
- `metabolite_edges_without_currency.csv`

## Optional: download and analyze the full BiGG model
```bash
python scripts/download_bigg_model.py
python scripts/analyze.py --sbml data/e_coli_core.xml
```
The SBML parser uses reaction reactants/products to build the same directed metabolite projection. Confirm the download and model metadata before submission; the curated subset remains available if the external endpoint changes.

## Required interpretation (complete after running)
In `results/interpretation_template.md`, fill in:
- Which network has more/fewer nodes and edges?
- How did N_D and n_D change after removing currency metabolites?
- Why might currency metabolites inflate apparent connectivity?
- What are two limitations of projecting reactions into a simple unweighted graph?
- Does the conclusion remain stable when numerical tolerances are changed?

Do not claim that a metabolite identified as a driver node is automatically a practical drug target or intervention point. The model is a linear network abstraction.

## Reproducibility and limitations
- This is a directed, unweighted metabolite projection. It ignores stoichiometric coefficients, reaction fluxes, kinetics, compartments and regulatory constraints.
- A reaction with multiple substrates/products is represented by all substrate-to-product pairs; this loses reaction-level hypergraph structure.
- Currency-metabolite exclusion is an analytical choice, not a universally correct biological rule.
- Exact controllability depends on edge weights and model structure. The included graph uses binary edges.
- SVD-based rank and eigenvalue multiplicity are sensitive to numerical tolerances. The script exports a tolerance sensitivity table.

## AI assistance disclosure
If you use AI assistance to generate or debug code or draft explanations, disclose it according to your course rules. Verify every result, understand the implementation, and write your own biological interpretation and synopsis.
