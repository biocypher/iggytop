"""Custom OntoWeaver transformer used by BaseAdapter's OntoWeaver-based node/edge generation.

Reproduces the `type:cdr3[:v_call]` id scheme (chain type lower-cased, V gene appended
only when present) previously produced by `BaseAdapter._generate_nodes_from_table`, so
the graph produced through OntoWeaver is a drop-in replacement for the hand-written
node/edge generation it replaces.
"""

from ontoweaver import base, transformer
from ontoweaver.make_value import ValueMaker


def _is_missing(value) -> bool:
    return value is None or isinstance(value, float) or str(value).strip().lower() in ("", "nan", "none")


class chain_id(base.Transformer):  # noqa: N801 (lowercase name required: referenced by name from the mapping YAML, matching OntoWeaver's own `cat`/`cat_format`/`map` convention)
    """Builds `type:cdr3[:v_call]` (lowercase type) from a `[type_col, cdr3_col, v_call_col]` triple.

    Rows missing the CDR3 column yield an empty string, which OntoWeaver's base Transformer
    filters out -- mirroring `BaseAdapter._generate_nodes_from_table`'s `.dropna(subset=unique_cols)`.
    """

    class ValueMaker(ValueMaker):
        """Computes the id value for a single row; see the enclosing class's docstring."""

        def __call__(self, columns, row, i):
            """Build the composite id string for one row (see class docstring)."""
            # By convention columns[0] is the "_type" column and columns[1] is the required
            # identifying column (e.g. cdr3): if it's missing, the entity doesn't exist for this row.
            required_col = columns[1]
            if required_col not in row or _is_missing(row[required_col]):
                yield ""
                return

            parts = []
            for c in columns:
                if c not in row or _is_missing(row[c]):
                    continue
                value = str(row[c])
                if c.endswith("_type"):
                    value = value.lower()
                parts.append(value)
            yield ":".join(parts)

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


transformer.register(chain_id)
