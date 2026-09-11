# Region 6 manual trace provenance

`Region_6_trace_control_points.tsv` is the version-controlled manual spatial annotation used by the Region 6 Mayassi-style feasibility analysis.

The points were entered during visual review of the Region 6 pass-QC cell map. The guide population was defined exactly as:

```r
region_spatial$RefAll_subtype_predicted.id == "Smooth muscle"
```

The line starts at the outer free edge, follows the outer smooth-muscle/serosal band through successive turns, crosses two deliberately labelled inter-turn bridges, and ends at the inner hook. The coordinates were rounded plotting positions; they were not produced by an optimization algorithm and are not histological ground truth.

Run `notebooks/B3a_Region6_trace_definition.ipynb` before projection. It prints every point, calculates segment diagnostics, and creates a numbered overlay. If a point needs revision, copy the TSV to a candidate file, edit the candidate, regenerate the numbered overlay, and obtain biological review before replacing the canonical file. Never tune the trace only to improve projection metrics.
