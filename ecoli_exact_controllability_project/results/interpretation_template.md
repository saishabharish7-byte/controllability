# Results interpretation — complete this after running the analysis

Do not invent or copy expected results. Use the values in `controllability_comparison.csv`.

1. **Network size:** Compare node count (N) and edge count (E) in the full and reduced networks.
2. **Driver nodes:** State the calculated ND and nD for each network. Explain the direction and magnitude of the change.
3. **Biological interpretation:** Explain why ubiquitous cofactors/currency metabolites can connect many otherwise separate reactions in a projected network.
4. **Limitations:** Discuss at least two: unweighted edges, reaction-to-pairwise-edge projection, no stoichiometry/flux/kinetics, no compartmentation, and numerical rank tolerance.
5. **Robustness:** Inspect `tolerance_sensitivity.csv`. Say whether ND changes across tolerances.
6. **Caution:** Driver metabolites are mathematical control inputs in this model, not automatically drug targets or experimentally validated intervention points.
