# PubChem View API

The PubChem View module provides access to experimental property data from PubChem's PUG View service. This includes extraction of experimental values, units, and comprehensive reference information.

::: provesid.pubchemview

::: provesid.pubchemview_parse

## Quick Start

```python
from provesid import PubChemView
from provesid.pubchemview import get_property_table

# Initialize the view client
view = PubChemView()

# Get experimental melting points for aspirin (CID 2244)
for prop in view.extract_property_data(2244, 'Melting Point'):
    p = prop.parsed
    print(f"As deposited: {prop.value!r}")
    # An entry may state a single value or a range, so read the bounds --
    # they are filled for both.
    print(f"Parsed:       {p.value_min}-{p.value_max} {p.unit}"
          f" = {p.value_min_si}-{p.value_max_si} {p.unit_si}")
    print(f"Reference:    {prop.reference}")

# Or as a DataFrame, with the values parsed into numbers and the references
# resolved to full citations
df = view.get_property_table(2244, 'Melting Point')
print(df[['StringWithMarkup', 'ExperimentalValue', 'Unit', 'ValueSI', 'UnitSI']])

# The same table without building a client first
table = get_property_table(2244, 'Boiling Point')
print(table)
```

## Available Property Types

The PubChem View service provides access to various experimental properties:

### Physical Properties
- **Melting Point** - Melting point temperatures
- **Boiling Point** - Boiling point temperatures  
- **Density** - Density measurements
- **Vapor Pressure** - Vapor pressure data
- **Solubility** - Solubility in various solvents
- **LogP** - Partition coefficient data
- **Viscosity** - Viscosity measurements
- **Refractive Index** - Refractive index values

### Spectroscopic Properties
- **UV/Vis Spectrum** - UV-Visible spectroscopy data
- **IR Spectrum** - Infrared spectroscopy data
- **NMR Spectrum** - Nuclear magnetic resonance data
- **Mass Spectrum** - Mass spectrometry data

### Safety and Toxicity
- **Flash Point** - Flash point temperatures
- **Auto-Ignition Temperature** - Auto-ignition data
- **LD50** - Lethal dose data
- **LC50** - Lethal concentration data

## Data Structures

### PropertyData Class

`PropertyData` is one value PubChem holds for one property, with its provenance:

```python
@dataclass
class PropertyData:
    value: str                          # the string exactly as PubChem wrote it
    unit: Optional[str] = None          # normalised; a copy of parsed.unit
    conditions: Optional[str] = None    # a copy of parsed.conditions
    reference: Optional[str] = None
    reference_number: Optional[int] = None
    description: Optional[str] = None
    name: Optional[str] = None
    heading: Optional[str] = None       # e.g. "Melting Point"
    parsed: Optional[ParsedValue] = None
```

`value` is always the original text, so nothing the parser cannot read is lost.
`parsed` holds the numbers recovered from it.

### ParsedValue: the numbers behind the prose

PubChem reports experimental properties as prose written by whoever deposited
them — `"138-140 °C"`, `"8.5X10-5 mm Hg at 25 °C"`, `"greater than or equal to
100 mg/mL"`, `"Vapor pressure, kPa at 20°C: 24"`. `parse_value()` turns that
into numbers. It is the only parser in the package for these strings, so a fix
reaches every method that reports a value.

```python
from provesid import parse_value

v = parse_value("2.47 cP at 20 °C", "Viscosity")
v.value        # 2.47
v.unit         # 'cP'
v.value_si     # 0.00247
v.unit_si      # 'Pa·s'
v.temperature_c  # 20.0  -- the condition, not the value
```

| Field | Meaning |
| --- | --- |
| `text` | The original string. Always populated. |
| `value` | The single number, or None for a range. |
| `value_min`, `value_max` | The bounds. Equal to `value` for a single value, so a numeric filter needs no special case for ranges. |
| `unit` | The unit as written, spelling normalised: `torr` and `mm Hg` both report as `mmHg`. |
| `value_si`, `value_min_si`, `value_max_si`, `unit_si` | The same quantity in SI units: `K`, `Pa`, `kg/m³`, `mol/m³`, `Pa·s`, `m²/s` or `N/m`. |
| `temperature_c` | The temperature the measurement was made at, in °C — a *condition*. For a melting point of 138 °C the value is 138 and this is None. |
| `operator` | `>`, `<`, `>=`, `<=` or `~` when the entry bounds the value rather than stating it. |
| `qualitative` | The word that replaced the number: `insoluble`, `miscible`, `negligible`. |
| `conditions` | Remaining qualifying text, such as a pressure or a `/Estimated/` note. |

The second argument is the heading, and it settles two things a string alone
cannot: a bare `138` under `"Melting Point"` is 138 °C, while a bare `1.19`
under `"LogP"` is dimensionless and complete. Passing None declines both hints.

Nothing is invented. An unrecognised unit is reported as written with no SI
conversion, and `%`, `ppm` and `ppb` are never converted at all — turning a
percentage into a concentration needs a density the string does not carry, so
`unit_si` is None there rather than a guess.

On a sample of 365 real value strings across ten compounds and ten properties,
90% yield a number and 8% a qualitative term; the rest were entries where
PubChem states no value at all (a pointer to one of its own external tables),
and those are now left out of the results rather than returned blank.

About three quarters of the numbers also carry an SI conversion. The rest are
either dimensionless — logP, pKa, refractive index, where `value` is the whole
answer — or compositions such as `%` and `ppm`, which cannot become
concentrations without a density. So test `unit_si`, or filter on it, rather
than assuming `value_si` is populated:

```python
# Comparing SI values only makes sense within one quantity
concentrations = table[table["UnitSI"] == "kg/m³"]
```

### Property tables

`get_property_table()` returns the parsed values as a DataFrame, one row per
entry, with the reference resolved:

```python
import pandas as pd

table = view.get_property_table(2244, "Melting Point")
table[["StringWithMarkup", "ExperimentalValue", "Unit", "ValueSI", "UnitSI"]]
#      StringWithMarkup  ExperimentalValue Unit  ValueSI UnitSI
# 0  275 °F (NTP, 1992)              275.0   °F   408.15      K
# 1             138-140                NaN   °C      NaN      K
# 2              135 °C              135.0   °C   408.15      K
```

Row 1 shows the range convention: `ExperimentalValue` is NaN because the entry
states no single number, while `ValueMin`/`ValueMax` hold 138 and 140 (and
`ValueMinSI`/`ValueMaxSI` hold 411.15 and 413.15 — a temperature conversion is
affine, so a range cannot simply be scaled).

Because the numbers are numbers, the frame is directly usable:

```python
# The mean melting point in kelvin, across entries reported in °C and °F alike
table["ValueSI"].mean()

# Only the entries that state a bound rather than a measurement
table[table["Operator"].notna()]

# Solubilities above 1 g/L, whatever unit they were deposited in
sol = view.get_property_table(2244, "Solubility")
sol[sol["ValueSI"] > 1.0]

sol.to_csv("solubility.csv", index=False)
```

Any PUG-View heading works, not only the experimental ones. The parser walks the
record rather than following one hard-coded path, so headings that live
elsewhere in PubChem's table of contents — `"GHS Classification"` under *Safety
and Hazards*, `"Drug Indication"` under *Drug and Medication Information* — now
return their values instead of an empty frame.

### Summaries and serialisation

`get_property_summary()` reports the raw strings alongside the numbers:

```python
summary = view.get_property_summary(2244, "Melting Point")
summary["values"]             # ['275 °F (NTP, 1992)', '138-140', '135 °C', ...]
summary["numeric_values"]     # [275.0, 135.0, 135.0, 135.0, 275.0, 275.0]
summary["numeric_values_si"]  # [408.15, 408.15, 408.15, 408.15, 408.15, 408.15]
summary["units"]              # ['°F', '°C']
summary["units_si"]           # ['K']
```

`numeric_values` mixes the two scales, which is exactly why
`numeric_values_si` exists: every entry agrees on 408.15 K.

`export_properties_to_dict()` flattens the same information into plain dicts,
ready for `json.dumps` or `pd.DataFrame`:

```python
rows = view.export_properties_to_dict(view.get_melting_point(2244))
rows[1]
# {'value': '138-140', 'unit': '°C', 'heading': 'Melting Point',
#  'numeric_value': None, 'value_min': 138.0, 'value_max': 140.0,
#  'value_min_si': 411.15, 'value_max_si': 413.15, 'unit_si': 'K', ...}
```

### Reference information

Each value carries the reference it came from. `get_property_table()` resolves
PubChem's reference numbers to full citations in `FullReference`:

```python
for _, row in view.get_property_table(2244, "Melting Point").iterrows():
    print(row["ExperimentalValue"], row["Unit"], "--", row["FullReference"][:60])
```

## Error Handling

The module provides specific exception classes:

```python
from provesid.pubchemview import PubChemViewError, PubChemViewNotFoundError

try:
    properties = view.get_experimental_properties(999999, 'Melting Point')
except PubChemViewNotFoundError:
    print("Compound or property not found")
except PubChemViewError as e:
    print(f"API error: {e}")
```

## Batch Processing

Process multiple compounds efficiently:

```python
def batch_extract_properties(cids, property_type):
    """Extract properties for multiple compounds"""
    results = {}
    view = PubChemView()
    
    for cid in cids:
        try:
            df = view.get_property_table(cid, property_type)
            if not df.empty:
                results[cid] = {
                    'count': len(df),
                    # ValueSI, not ExperimentalValue: entries deposited in
                    # different units only average meaningfully in SI.
                    'mean_value_si': df['ValueSI'].mean(),
                    'unit_si': df['UnitSI'].dropna().unique().tolist(),
                    'units_as_deposited': df['Unit'].dropna().unique().tolist(),
                }
        except Exception as e:
            results[cid] = {'error': str(e)}
    
    return results

# Process multiple compounds
cids = [2244, 2519, 3672]  # Aspirin, Caffeine, Ibuprofen
melting_data = batch_extract_properties(cids, 'Melting Point')
```

## Integration with Other APIs

Combine with PubChem compound data:

```python
from provesid import PubChemAPI, PubChemView

def comprehensive_compound_analysis(cid):
    """Get both computed and experimental data"""
    api = PubChemAPI()
    view = PubChemView()
    
    # Get computed properties
    computed = api.get_compound_properties(
        [cid], 
        ['MolecularWeight', 'MolecularFormula', 'ConnectivitySMILES']
    )
    
    # Get experimental properties
    experimental = {}
    for prop_type in ['Melting Point', 'Boiling Point', 'Density']:
        try:
            df = view.get_property_table(cid, prop_type)
            if not df.empty:
                experimental[prop_type] = df
        except:
            pass
    
    return {
        'computed': computed,
        'experimental': experimental
    }

# Analyze aspirin
analysis = comprehensive_compound_analysis(2244)
```

## Performance Considerations

### Rate Limiting

The PubChemView client includes automatic rate limiting:

```python
# Adjust request frequency for large batch jobs
view = PubChemView(pause_time=1.0)  # 1 second between requests

# For development/testing with faster requests
view_fast = PubChemView(pause_time=0.1)  # 100ms between requests
```

### Caching

Consider implementing caching for frequently accessed data:

```python
import pickle
from pathlib import Path

def cached_property_extraction(cid, property_type, cache_dir='cache'):
    """Extract properties with file-based caching"""
    cache_path = Path(cache_dir) / f"{cid}_{property_type}.pkl"
    
    if cache_path.exists():
        with open(cache_path, 'rb') as f:
            return pickle.load(f)
    
    # Extract fresh data
    view = PubChemView()
    df = view.get_property_table(cid, property_type)
    
    # Cache the results
    cache_path.parent.mkdir(exist_ok=True)
    with open(cache_path, 'wb') as f:
        pickle.dump(df, f)
    
    return df
```

## See Also

- [PubChem API](pubchem.md) - For computed compound properties
- [PubChem Tutorial](../examples/pubchem/pubchem_tutorial.md) - Complete tutorial with detailed examples
