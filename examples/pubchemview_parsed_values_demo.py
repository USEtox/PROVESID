"""
Demo script for the typed values PubChemView now returns.

PubChem reports experimental properties as prose: "138-140", "275 °F",
"8.5X10-5 mm Hg at 25 °C", "greater than or equal to 100 mg/mL". This script
shows the single parser behind ``provesid.parse_value`` turning those into
numbers, and what that makes possible — comparing two compounds, averaging
entries deposited in different units, filtering on a threshold.

Run it against the live API; results are cached, so a second run is instant.
"""

from provesid import PubChemView, parse_value

print("=" * 80)
print("PubChemView Parsed Values Demo")
print("=" * 80)
print()

# ---------------------------------------------------------------------------
# 1. The parser on its own
# ---------------------------------------------------------------------------
print("1. Strings PubChem actually returns, parsed")
print("-" * 80)

samples = [
    ("138-140", "Melting Point"),
    ("275 °F (NTP, 1992)", "Melting Point"),
    ("8.5X10-5 mm Hg at 25 °C", "Vapor Pressure"),
    ("In water, 4,600 mg/L at 25 °C", "Solubility"),
    ("greater than or equal to 100 mg/mL", "Solubility"),
    ("Vapor pressure, kPa at 20°C: 24", "Vapor Pressure"),
    ("2.47 cP at 20 °C", "Viscosity"),
    ("Insoluble in water", "Solubility"),
    ("1.19", "LogP"),
]

for text, heading in samples:
    v = parse_value(text, heading)
    if v.is_range:
        number = f"{v.value_min}–{v.value_max}"
    elif v.value is not None:
        number = f"{v.operator or ''}{v.value}"
    else:
        number = v.qualitative or "—"
    si = "" if v.value_si is None else f"  = {v.value_si:.6g} {v.unit_si}"
    at = "" if v.temperature_c is None else f"  (at {v.temperature_c:g} °C)"
    print(f"  {text!r:38} -> {number} {v.unit or ''}{si}{at}")
print()
print("  Note the last three: a bound is not a measurement, 'insoluble' is a")
print("  result rather than a missing value, and a logP is dimensionless.")
print()

# ---------------------------------------------------------------------------
# 2. Why the SI columns exist
# ---------------------------------------------------------------------------
print("2. Entries deposited in different units, compared")
print("-" * 80)

view = PubChemView()
table = view.get_property_table(2244, "Melting Point")
print(table[["StringWithMarkup", "ExperimentalValue", "Unit",
             "ValueSI", "UnitSI"]].to_string(index=False))
print()
print(f"  Mean melting point: {table['ValueSI'].mean():.2f} K")
print("  Every entry agrees, although they were deposited in °C and °F --")
print("  which is the cross-check that the conversion is right both ways.")
print()

# ---------------------------------------------------------------------------
# 3. A range keeps both ends
# ---------------------------------------------------------------------------
print("3. Ranges")
print("-" * 80)

ranges = table[table["ExperimentalValue"].isna() & table["ValueMin"].notna()]
for _, row in ranges.iterrows():
    print(f"  {row['StringWithMarkup']!r} -> "
          f"{row['ValueMin']}–{row['ValueMax']} {row['Unit']} "
          f"= {row['ValueMinSI']}–{row['ValueMaxSI']} {row['UnitSI']}")
print("  ExperimentalValue is NaN for these: the entry states no single number.")
print("  ValueMin/ValueMax are filled for single values too, so a numeric")
print("  filter needs no special case.")
print()

# ---------------------------------------------------------------------------
# 4. Filtering and comparing compounds
# ---------------------------------------------------------------------------
print("4. Comparing compounds on a numeric property")
print("-" * 80)

compounds = {"Aspirin": 2244, "Ethanol": 702, "Caffeine": 2519, "Ibuprofen": 3672}
for name, cid in compounds.items():
    entries = view.get_property_table(cid, "Melting Point")
    numeric = entries["ValueSI"].dropna()
    if numeric.empty:
        print(f"  {name:<10} no numeric melting point")
        continue
    print(f"  {name:<10} {numeric.mean():7.2f} K  "
          f"({numeric.mean() - 273.15:6.1f} °C, from {len(numeric)} entries)")
print()

print("  Solubility entries above 1 g/L, whatever unit they were deposited in:")
sol = view.get_property_table(2244, "Solubility")
# Comparing SI values only makes sense within one quantity, so the mass
# concentrations are selected first -- a percentage has no SI equivalent and a
# molarity converts to mol/m³, not kg/m³.
concentrations = sol[sol["UnitSI"] == "kg/m³"]
above = concentrations[concentrations["ValueSI"] > 1.0]
for _, row in above.head(5).iterrows():
    print(f"    {row['StringWithMarkup'][:56]!r:58} = {row['ValueSI']:.4g} kg/m³")
skipped = len(sol) - len(concentrations)
if skipped:
    print(f"    ({skipped} further entries are percentages or qualitative, so they")
    print("     carry no comparable SI value -- reported as such, not guessed at)")
print()

# ---------------------------------------------------------------------------
# 5. Any heading, not only the experimental ones
# ---------------------------------------------------------------------------
print("5. Headings outside Experimental Properties")
print("-" * 80)

for heading in ["GHS Classification", "Drug Indication"]:
    try:
        rows = view.get_property_table(2244, heading)
    except Exception as exc:
        print(f"  {heading}: {type(exc).__name__}: {exc}")
        continue
    print(f"  {heading}: {len(rows)} rows")
    for _, row in rows.head(2).iterrows():
        text = row["StringWithMarkup"]
        if text:
            print(f"    {text[:70]!r}")
print()
print("  These live elsewhere in PubChem's table of contents. The parser walks")
print("  the record instead of following one fixed path, so they work too.")
print()

print("=" * 80)
print("Demo completed!")
print("=" * 80)
