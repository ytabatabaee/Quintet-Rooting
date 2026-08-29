# Changelog

## QR-STAR v1.3.0

### Added
- Direct support for multi-copy gene-family trees through integrated DISCO decomposition.
- Gene-copy to species mapping using delimiters or explicit mapping files.
- Optional saving of decomposed DISCO trees.
- Multi-copy preprocessing diagnostics.

### Behavior
- DISCO subtrees with fewer than five species are discarded because they cannot contribute quintets.
- Missing-taxon normalization is enabled automatically in multi-copy mode.
- Existing single-copy QR-STAR behavior is unchanged.
