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

from ontoweaver import base, congregate, transformer
from ontoweaver.make_value import ValueMaker


def _patch_congregate_performance() -> None:
    """Replace `ontoweaver.congregate.Congregate.call` with an O(k) equivalent.

    The upstream implementation does `self._duplicates[elem] = self._duplicates.get(elem, []) + [elem]`,
    which rebuilds the whole list from scratch on every duplicate occurrence of the same node/edge id --
    O(k^2) for k duplicates of one id. Our `database`/`PMID` nodes have a near-constant id across an
    entire adapter's table (e.g. the same database name on every row), so on large adapters (VDJDB) this
    makes reconciliation impractically slow (multiple minutes and climbing). `setdefault(...).append(...)`
    produces the exact same final dict/order in O(k). Safe to remove once fixed upstream.
    """

    def _call(self, biocypher_tuples):
        for t in biocypher_tuples:
            elem = self._elem_cls.from_tuple(t, serializer=self.serializer)
            self._duplicates.setdefault(elem, []).append(elem)
            yield elem

    congregate.Congregate.call = _call


_patch_congregate_performance()


def _is_missing(value) -> bool:
    return value is None or isinstance(value, float) or str(value).strip().lower() in ("", "nan", "none")


def chain_component(row, type_col: str, cdr3_col: str, v_col: str) -> str | None:
    """Build `locus:cdr3[:v_call]` (lowercase locus) for one chain, or None if its CDR3 is missing."""
    if cdr3_col not in row or _is_missing(row[cdr3_col]):
        return None
    parts = [str(row[type_col]).lower(), str(row[cdr3_col])]
    if v_col in row and not _is_missing(row[v_col]):
        parts.append(str(row[v_col]))
    return ":".join(parts)


def epitope_component(row) -> str:
    """Build `epitope:<iedb_iri>`. Always present: BaseAdapter.table requires a non-null epitope."""
    return f"epitope:{row['iedb_iri']}"


def _mhc_fields(row) -> tuple[str, str, str, bool]:
    """Return (class, gene_1, gene_2) as blank-safe strings, plus whether any of them is present."""
    raw = [row.get("MHC_class"), row.get("MHC_gene_1"), row.get("MHC_gene_2")]
    present = any(not _is_missing(f) for f in raw)
    safe = tuple("" if _is_missing(f) else str(f) for f in raw)
    return (*safe, present)


def _chains_raw(row) -> str:
    """Build `<chain_1-or-none>|<chain_2-or-none>`, the shared, unprefixed core of receptor_complex's id."""
    c1 = chain_component(row, "chain_1_type", "chain_1_cdr3", "chain_1_v_call") or "none"
    c2 = chain_component(row, "chain_2_type", "chain_2_cdr3", "chain_2_v_call") or "none"
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
    """Builds `locus:cdr3[:v_call]` from a `[type_col, cdr3_col, v_call_col]` column triple.

    Rows missing the CDR3 column yield an empty string, which OntoWeaver's base Transformer
    filters out -- mirroring BaseAdapter's old `.dropna(subset=unique_cols)` behavior.
    """

    class ValueMaker(ValueMaker):
        """Computes the id value for a single row; see the enclosing class's docstring."""

        def __call__(self, columns, row, i):
            type_col, cdr3_col, v_col = columns
            yield chain_component(row, type_col, cdr3_col, v_col) or ""

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
