"""Custom OntoWeaver transformers used by BaseAdapter's OntoWeaver-based node/edge generation.

The graph has a layered hub-and-spoke shape (binding -> {receptor_complex, pmhc, database, PMID},
receptor_complex -> {chain_1, chain_2}, chain_N -> {v_gene, j_gene}, pmhc -> {epitope, mhc}). OntoWeaver
maps one row-subject with radiating edges per mapping pass, so each hub in that hierarchy needs its own
mapping pass (see the `ontoweaver_mapping_*.yaml` files), reconciled together afterwards. Composite
entities (receptor_complex, pmhc, binding) don't have a column holding their id directly -- their id is
built from the same row data whether they're being used as a pass's row subject or referenced as another
pass's target, so the component builders below are shared by both roles to guarantee the two passes agree
on the same id string.
"""

from ontoweaver import base, transformer
from ontoweaver.make_value import ValueMaker


def _is_missing(value) -> bool:
    return value is None or isinstance(value, float) or str(value).strip().lower() in ("", "nan", "none")


def chain_component(row, type_col: str, cdr3_col: str, v_col: str, j_col: str) -> str | None:
    """Build `locus:cdr3[:v:v_call][:j:j_call]` (lowercase locus) for one chain, or None if its CDR3 is missing.

    v_call and j_call are each tagged with their own `v:`/`j:` marker (rather than appended bare)
    so that a chain missing one of them can never collide with a chain missing the other -- e.g.
    "locus:cdr3:v:X" (has v_call, no j_call) can't be confused with "locus:cdr3:j:X" (the reverse).
    This keeps a chain that's missing a gene call as a distinct node from an otherwise-identical
    chain that has it, instead of merging their (possibly conflicting) `complete` flags together.
    """
    if cdr3_col not in row or _is_missing(row[cdr3_col]):
        return None
    parts = [str(row[type_col]).lower(), str(row[cdr3_col])]
    if v_col in row and not _is_missing(row[v_col]):
        parts.append(f"v:{row[v_col]}")
    if j_col in row and not _is_missing(row[j_col]):
        parts.append(f"j:{row[j_col]}")
    return ":".join(parts)


def epitope_component(row) -> str:
    """Build `epitope:<iedb_iri>`. Always present: BaseAdapter.table requires a non-null epitope."""
    if _is_missing(row["iedb_iri"]):
        raise ValueError("epitope_component: 'iedb_iri' is missing, but BaseAdapter.table requires it")
    return f"epitope:{row['iedb_iri']}"


def _mhc_fields(row) -> tuple[str, str, str, bool]:
    """Return (class, gene_1, gene_2) as blank-safe strings, plus whether any of them is present."""
    raw = [row.get("MHC_class"), row.get("MHC_gene_1"), row.get("MHC_gene_2")]
    present = any(not _is_missing(f) for f in raw)
    safe = tuple("" if _is_missing(f) else str(f) for f in raw)
    return (*safe, present)


def _chains_raw(row) -> str:
    """Build `<chain_1-or-none>|<chain_2-or-none>`, the shared, unprefixed core of receptor_complex's id."""
    c1 = chain_component(row, "chain_1_type", "chain_1_cdr3", "chain_1_v_call", "chain_1_j_call") or "none"
    c2 = chain_component(row, "chain_2_type", "chain_2_cdr3", "chain_2_v_call", "chain_2_j_call") or "none"
    return f"{c1}|{c2}"


def mhc_component(row) -> str | None:
    """Build `mhc:<class>|<gene_1>|<gene_2>`, or None if all three fields are missing."""
    cls, g1, g2, present = _mhc_fields(row)
    if not present:
        return None
    return f"mhc:{cls}|{g1}|{g2}"


def pmhc_component(row) -> str:
    """Build `pmhc:<iedb_iri>|<class>|<gene_1>|<gene_2>`.

    Uses the raw epitope/mhc values directly rather than nesting their own prefixed
    `epitope:`/`mhc:` ids, so the id doesn't repeat each value's type name. Always present
    (epitope is mandatory); the three mhc fields are blank, not omitted, when absent.
    """
    cls, g1, g2, _ = _mhc_fields(row)
    return f"pmhc:{row['iedb_iri']}|{cls}|{g1}|{g2}"


def receptor_complex_component(row) -> str:
    """Build `receptor_complex:<chain_1-or-none>|<chain_2-or-none>`. Always present: at least one chain is mandatory."""
    return f"receptor_complex:{_chains_raw(row)}"


def binding_component(row) -> str:
    """Build `binding:<chain_1>|<chain_2>|<iedb_iri>|<class>|<gene_1>|<gene_2>`.

    Built from the same raw chain/epitope/mhc values as receptor_complex/pmhc, rather than nesting
    their prefixed ids, so the id doesn't repeat "receptor_complex:"/"pmhc:" internally. Deliberately
    excludes source/PMID/database so that the same receptor-complex-recognizes-pmhc claim reported
    by different sources collapses onto one node.
    """
    cls, g1, g2, _ = _mhc_fields(row)
    return f"binding:{_chains_raw(row)}|{row['iedb_iri']}|{cls}|{g1}|{g2}"


class _RowValueTransformer(base.Transformer):
    """Base for transformers whose id is computed from the whole row rather than from `columns`.

    Subclasses set `component_fn` to a `row -> str | None` function; a None result yields an empty
    string, which OntoWeaver's base Transformer filters out -- skipping the node/edge for that row.
    """

    component_fn: staticmethod

    class ValueMaker(ValueMaker):
        """Computes the id value for a single row by delegating to the outer class's `component_fn`."""

        def __init__(self, component_fn, raise_errors: bool = True):
            """Store the row -> id component function to use."""
            self.component_fn = component_fn
            super().__init__(raise_errors)

        def __call__(self, columns, row, i):
            yield self.component_fn(row) or ""

    def __init__(  # noqa: PLR0913 (signature mirrors ontoweaver.base.Transformer's, e.g. `cat`)
        self,
        properties_of,
        label_maker=None,
        branching_properties=None,
        columns=None,
        output_validator=None,
        multi_type_dict=None,
        raise_errors=True,
        **kwargs,
    ):
        """Initialize the transformer (see the class docstring)."""
        self.value_maker = self.ValueMaker(self.component_fn, raise_errors=raise_errors)
        super().__init__(
            properties_of,
            self.value_maker,
            label_maker,
            branching_properties,
            columns,
            output_validator,
            multi_type_dict,
            raise_errors=raise_errors,
            **kwargs,
        )


class chain_id(base.Transformer):  # noqa: N801 (lowercase name required: referenced by name from the mapping YAML, matching OntoWeaver's own `cat`/`cat_format`/`map` convention)
    """Builds `locus:cdr3[:v:v_call][:j:j_call]` from a `[type_col, cdr3_col, v_call_col, j_call_col]` column quadruple.

    Rows missing the CDR3 column yield an empty string, which OntoWeaver's base Transformer
    filters out -- mirroring BaseAdapter's old `.dropna(subset=unique_cols)` behavior.
    """

    class ValueMaker(ValueMaker):
        """Computes the id value for a single row; see the enclosing class's docstring."""

        def __call__(self, columns, row, i):
            type_col, cdr3_col, v_col, j_col = columns
            yield chain_component(row, type_col, cdr3_col, v_col, j_col) or ""

    def __init__(  # noqa: PLR0913 (signature mirrors ontoweaver.base.Transformer's, e.g. `cat`)
        self,
        properties_of,
        label_maker=None,
        branching_properties=None,
        columns=None,
        output_validator=None,
        multi_type_dict=None,
        raise_errors=True,
        **kwargs,
    ):
        """Initialize the transformer (see the class docstring)."""
        self.value_maker = self.ValueMaker(raise_errors=raise_errors)
        super().__init__(
            properties_of,
            self.value_maker,
            label_maker,
            branching_properties,
            columns,
            output_validator,
            multi_type_dict,
            raise_errors=raise_errors,
            **kwargs,
        )


class prefixed_id(base.Transformer):  # noqa: N801
    """Builds `<prefix>:<value>` from a single column; empty (skipped) when that column's value is missing.

    `prefix` is a mapping-level parameter (e.g. `prefixed_id: {columns: [chain_1_v_call], prefix: v_gene, ...}`),
    used for reference/leaf entities whose id is the whole of their informational content (v_gene, j_gene, PMID).
    """

    class ValueMaker(ValueMaker):
        """Computes the id value for a single row; see the enclosing class's docstring."""

        def __init__(self, prefix: str, raise_errors: bool = True):
            """Store the prefix to prepend to the column value."""
            self.prefix = prefix
            super().__init__(raise_errors)

        def __call__(self, columns, row, i):
            (col,) = columns
            if col not in row or _is_missing(row[col]):
                yield ""
                return
            yield f"{self.prefix}:{row[col]}"

    def __init__(  # noqa: PLR0913 (signature mirrors ontoweaver.base.Transformer's, e.g. `cat`)
        self,
        properties_of,
        label_maker=None,
        branching_properties=None,
        columns=None,
        output_validator=None,
        multi_type_dict=None,
        raise_errors=True,
        prefix=None,
        **kwargs,
    ):
        """Initialize the transformer (see the class docstring)."""
        self.value_maker = self.ValueMaker(prefix, raise_errors=raise_errors)
        super().__init__(
            properties_of,
            self.value_maker,
            label_maker,
            branching_properties,
            columns,
            output_validator,
            multi_type_dict,
            raise_errors=raise_errors,
            **kwargs,
        )


class mhc_id(_RowValueTransformer):  # noqa: N801
    """Builds `mhc:<class>|<gene_1>|<gene_2>`; empty (skipped) when no MHC field is present."""

    component_fn = staticmethod(mhc_component)


class pmhc_id(_RowValueTransformer):  # noqa: N801
    """Builds `pmhc:<epitope>|<mhc-or-none>`; always present since epitope is mandatory."""

    component_fn = staticmethod(pmhc_component)


class receptor_complex_id(_RowValueTransformer):  # noqa: N801
    """Builds `receptor_complex:<chain_1-or-none>|<chain_2-or-none>`; always present."""

    component_fn = staticmethod(receptor_complex_component)


class binding_id(_RowValueTransformer):  # noqa: N801
    """Builds `binding:<receptor_complex>|<pmhc>`; always present."""

    component_fn = staticmethod(binding_component)


for _cls in (chain_id, prefixed_id, mhc_id, pmhc_id, receptor_complex_id, binding_id):
    transformer.register(_cls)
