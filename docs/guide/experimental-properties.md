# Experimental properties from PubChem

PubChem's PUG-View service holds the experimental data deposited for a
compound — melting and boiling points, vapour pressure, solubility, logP,
density, viscosity and much else — each value with the reference it came from.
[`PubChemView`](../api/pubchemview.md) fetches it and turns the prose it is
written in back into numbers. This is an online service; see
[Network behaviour](network.md) and [Caching](caching.md) for how it is asked.

```python
from provesid import PubChemView

view = PubChemView()

table = view.get_property_table(2244, "Melting Point")      # aspirin
table[["StringWithMarkup", "ExperimentalValue", "Unit", "ValueSI", "UnitSI"]]
#      StringWithMarkup  ExperimentalValue Unit  ValueSI UnitSI
# 0  275 °F (NTP, 1992)              275.0   °F   408.15      K
# 1             138-140                NaN   °C      NaN      K
# 2              135 °C              135.0   °C   408.15      K
```

`view.get_available_properties(cid)` lists the headings PubChem holds for a
compound. Any PUG-View heading works, not only the experimental ones:
`"GHS Classification"` and `"Drug Indication"` return their values too.

## One value, as deposited and as a number

`extract_property_data()` returns one `PropertyData` per value PubChem holds:

| Field | Meaning |
|---|---|
| `value` | the string exactly as PubChem wrote it |
| `unit`, `conditions` | copies of `parsed.unit` and `parsed.conditions` |
| `reference`, `reference_number` | where the value came from |
| `description`, `name`, `heading` | PubChem's labels, e.g. heading `"Melting Point"` |
| `parsed` | a `ParsedValue`: the numbers recovered from `value` |

`value` is always the original text, so nothing the parser cannot read is
lost.

## Parsing the prose

PubChem reports experimental properties as prose written by whoever deposited
them: `"138-140 °C"`, `"8.5X10-5 mm Hg at 25 °C"`, `"greater than or equal to
100 mg/mL"`, `"Vapor pressure, kPa at 20°C: 24"`. `parse_value()` turns that
into numbers, and is the only parser for these strings in the package:

```python
from provesid import parse_value

v = parse_value("2.47 cP at 20 °C", "Viscosity")
v.value          # 2.47
v.unit           # 'cP'
v.value_si       # 0.00247
v.unit_si        # 'Pa·s'
v.temperature_c  # 20.0  -- the condition, not the value
```

| Field | Meaning |
| --- | --- |
| `text` | the original string; always populated |
| `value` | the single number, or None for a range |
| `value_min`, `value_max` | the bounds; equal to `value` for a single value, so a numeric filter needs no special case for ranges |
| `unit` | the unit as written, spelling normalised: `torr` and `mm Hg` both report as `mmHg` |
| `value_si`, `value_min_si`, `value_max_si`, `unit_si` | the same quantity in SI units: `K`, `Pa`, `kg/m³`, `mol/m³`, `Pa·s`, `m²/s` or `N/m` |
| `temperature_c` | the temperature the measurement was made at, in °C — a *condition*; for a melting point of 138 °C the value is 138 and this is None |
| `operator` | `>`, `<`, `>=`, `<=` or `~` when the entry bounds the value rather than stating it |
| `qualitative` | the word that replaced the number: `insoluble`, `miscible`, `negligible` |
| `conditions` | remaining qualifying text, such as a pressure or an `/Estimated/` note |

The second argument is the heading, and it settles what a string alone cannot:
a bare `138` under `"Melting Point"` is 138 °C, while a bare `1.19` under
`"LogP"` is dimensionless and complete. Passing None declines both hints.

Nothing is invented. An unrecognised unit is reported as written with no SI
conversion, and `%`, `ppm` and `ppb` are never converted — turning a
percentage into a concentration needs a density the string does not carry.

On a sample of 365 real value strings across ten compounds and ten properties,
90% yield a number and 8% a qualitative term; the rest are entries where
PubChem states no value at all, which are left out of the results. About three
quarters of the numbers also carry an SI conversion; the rest are
dimensionless (logP, pKa, refractive index) or compositions such as `%`. So
filter on `unit_si` rather than assuming `value_si` is populated.

## Property tables

`get_property_table()` returns the parsed values as a DataFrame, one row per
entry, with PubChem's reference numbers resolved to full citations in
`FullReference`. In the table above, row 1 shows the range convention:
`ExperimentalValue` is NaN because the entry states no single number, while
`ValueMin`/`ValueMax` hold 138 and 140 and `ValueMinSI`/`ValueMaxSI` hold
411.15 and 413.15.

Because the numbers are numbers, the frame is directly usable:

```python
# The mean melting point in kelvin, across entries reported in °C and °F alike
table["ValueSI"].mean()

# Only the entries that state a bound rather than a measurement
table[table["Operator"].notna()]

# Solubilities above 1 kg/m³, whatever unit they were deposited in
sol = view.get_property_table(2244, "Solubility")
sol[(sol["UnitSI"] == "kg/m³") & (sol["ValueSI"] > 1.0)]
```

`provesid.get_property_table(cid, heading)` does the same without building a
client first.

## Summaries and plain dicts

`get_property_summary()` reports the raw strings alongside the numbers:

```python
summary = view.get_property_summary(2244, "Melting Point")
summary["values"]             # ['275 °F (NTP, 1992)', '138-140', '135 °C', ...]
summary["numeric_values"]     # [275.0, 135.0, 135.0, 135.0, 275.0, 275.0]
summary["numeric_values_si"]  # [408.15, 408.15, 408.15, 408.15, 408.15, 408.15]
summary["units_si"]           # ['K']
```

`numeric_values` mixes the two scales, which is why `numeric_values_si`
exists: every entry agrees on 408.15 K.

`export_properties_to_dict()` flattens a list of `PropertyData` into plain
dicts, ready for `json.dumps` or `pd.DataFrame`:

```python
rows = view.export_properties_to_dict(view.get_melting_point(2244))
rows[1]
# {'value': '138-140', 'unit': '°C', 'heading': 'Melting Point',
#  'numeric_value': None, 'value_min': 138.0, 'value_max': 140.0,
#  'value_min_si': 411.15, 'value_max_si': 413.15, 'unit_si': 'K', ...}
```

The [PubChem View tutorial](../examples/pubchemview/pubchem_view_tutorial.md)
walks through a full example.
