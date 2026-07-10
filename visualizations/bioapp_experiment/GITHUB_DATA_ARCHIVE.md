# GitHub data archive

The three Phase 4 `assigned_expression/matrix.mtx` files exceed GitHub's
regular 100 MiB file limit. The repository therefore stores lossless gzip
copies named `matrix.mtx.gz` in the same directories.

To restore the original Matrix Market file:

```bash
gzip -dk matrix.mtx.gz
```

The uncompressed files remain in the local working dataset but are ignored by
Git. No matrix rows, columns, or values are removed by this archive step.
