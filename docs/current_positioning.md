# Current Positioning

This file summarizes the maintained project scope after cleanup.

## Current Mainline
- Pipeline: `Stage1 -> Stage3 -> Stage4`
- Stage4 outputs: CytoSPACE baseline and SVTuner + CytoSPACE route2
- Core method: type-aware / mismatch-aware preprocessing before CytoSPACE mapping

## Core Contributions
1. Stage3 identifies unsupported or profile-masked cell types before spatial assignment.
2. Route2 modifies the input cell pool and labels without changing the CytoSPACE solver.
3. The maintained workflow stops after the two Stage4 mapping outputs used by the current analyses.

## Out Of Current Scope
- Stage2 SVG-aware / SVG+HVG exploration is historical.
- Stage5/Stage6/Stage7 post-processing scripts have been removed from the maintained workflow.
- Deprecated simulation-only reporting flows should not be reintroduced unless explicitly needed.
