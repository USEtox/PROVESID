"""
Tests for the single PUG-View value parser.

Every case here is a string shape PubChem actually returns, so the file doubles
as a record of what the free-text values look like. Nothing touches the network:
the parser takes a string and a heading and returns numbers.
"""

import pytest

from provesid.pubchemview_parse import (
    UNIT_ALIASES,
    UNIT_TO_SI,
    ParsedValue,
    parse_value,
)


# --------------------------------------------------------------------------
# Numbers
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_single_value_with_unit():
    """The base case, and the one the old parser returned as a string."""
    parsed = parse_value("135 °C", "Melting Point")

    assert parsed.value == 135.0
    assert parsed.unit == "°C"
    assert isinstance(parsed.value, float)


@pytest.mark.unit
def test_single_value_fills_both_bounds():
    """
    A single value is a degenerate range.

    Filling both bounds is what lets a caller filter on ``value_min`` and
    ``value_max`` without first asking which shape it got.
    """
    parsed = parse_value("135 °C", "Melting Point")

    assert parsed.value_min == parsed.value_max == 135.0
    assert parsed.is_numeric
    assert not parsed.is_range


@pytest.mark.unit
def test_range():
    """``"138-140"`` is two numbers, and neither of them is "the" value."""
    parsed = parse_value("138-140", "Melting Point")

    assert parsed.value is None
    assert (parsed.value_min, parsed.value_max) == (138.0, 140.0)
    assert parsed.is_range


@pytest.mark.unit
@pytest.mark.parametrize("text", ["20 to 25 °C", "20-25 °C", "20–25 °C", "25-20 °C"])
def test_range_spellings(text):
    """PubChem writes a range with a hyphen, an en dash or the word "to"."""
    parsed = parse_value(text, "Melting Point")

    assert (parsed.value_min, parsed.value_max) == (20.0, 25.0), text


@pytest.mark.unit
def test_x_notation():
    """
    PubChem writes scientific notation as ``8.5X10-5``.

    Rewriting it to ``8.5e-5`` before matching is what keeps the pattern list
    short enough to read.
    """
    parsed = parse_value("8.5X10-5 mm Hg at 25 °C", "Vapor Pressure")

    assert parsed.value == pytest.approx(8.5e-5)
    assert parsed.unit == "mmHg"


@pytest.mark.unit
def test_thousands_separator():
    """
    ``"In water, 4,600 mg/L"`` is one number, not a 4 followed by a unit.

    The comma after "water" must survive, so only a comma between a digit and
    exactly three digits is removed.
    """
    parsed = parse_value("In water, 4,600 mg/L at 25 °C", "Solubility")

    assert parsed.value == 4600.0
    assert parsed.unit == "mg/L"
    assert parsed.temperature_c == 25.0


@pytest.mark.unit
def test_negative_value():
    """LogP is routinely negative, and a leading minus is not a range."""
    parsed = parse_value("-2.6", "LogP")

    assert parsed.value == -2.6
    assert parsed.unit is None


@pytest.mark.unit
def test_no_number_at_all():
    """A string with no number gives None, never a zero."""
    parsed = parse_value("no numeric value here", "Melting Point")

    assert parsed.value is None
    assert parsed.value_min is None
    assert not parsed.is_numeric


@pytest.mark.unit
@pytest.mark.parametrize("text", [None, "", "   "])
def test_empty_input(text):
    """None and empty strings are answered, not raised on."""
    parsed = parse_value(text, "Melting Point")

    assert isinstance(parsed, ParsedValue)
    assert parsed.value is None


@pytest.mark.unit
def test_original_string_is_always_kept():
    """
    Nothing the parser fails to understand is lost.

    ``text`` is the caller's escape hatch: whatever the patterns miss is still
    there to be read by eye or by a better parser later.
    """
    weird = "1.95 mmÂ²/s at 20Â °C /mojibake/"
    assert parse_value(weird, "Viscosity").text == weird


# --------------------------------------------------------------------------
# Units
# --------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("text,unit", [
    ("760 mm Hg", "mmHg"),
    ("760 mmHg", "mmHg"),
    ("760 MM HG", "mmHg"),
    ("1.5 torr", "mmHg"),
    ("0.79 g/cm3", "g/cm³"),
    ("0.79 G/CM3", "g/cm³"),
    ("0.79 g/cc", "g/cm³"),
])
def test_unit_spellings_are_normalised(text, unit):
    """
    A caller grouping by unit must not see three spellings of one unit.

    Matching is case-insensitive for the known spellings, which is why
    ``MM HG`` resolves.
    """
    assert parse_value(text, "Vapor Pressure").unit == unit


@pytest.mark.unit
def test_a_single_letter_is_a_unit_only_as_a_whole_token():
    """
    The ``c`` of ``cP`` is not Celsius.

    This is the bug that a permissive unit pattern plus a ``'c': '°C'`` alias
    produced: viscosity in centipoise came back as a temperature in kelvin.
    """
    parsed = parse_value("2.47cP at 20 °C", "Viscosity")

    assert parsed.unit == "cP"
    assert parsed.unit_si == "Pa·s"
    assert parsed.value_si == pytest.approx(2.47e-3)


@pytest.mark.unit
def test_molar_m_is_not_guessed():
    """
    A bare ``M`` is left as written rather than read as mol/L.

    It is indistinguishable from metres, from the M of a molecular weight, and
    from a stray capital, so guessing cost more than it fixed.
    """
    assert "m" not in UNIT_ALIASES
    assert parse_value("0.5 mmol/L", "Solubility").unit == "mmol/L"


@pytest.mark.unit
def test_bracketed_unit():
    """PubChem brackets a unit when the depositor's software did."""
    assert parse_value("0.05 [mmHg]", "Vapor Pressure").unit == "mmHg"
    assert parse_value("1.2 [ug/mL] (est)", "Solubility").unit == "µg/mL"


@pytest.mark.unit
def test_unrecognised_unit_is_reported_as_written():
    """
    An unknown unit is kept, not dropped.

    Losing it would make the number meaningless; keeping it lets the caller
    decide.
    """
    parsed = parse_value("5 furlongs", "Some Heading")

    assert parsed.value == 5.0
    assert parsed.unit == "furlongs"
    assert parsed.unit_si is None
    assert parsed.value_si is None


@pytest.mark.unit
def test_english_words_are_not_units():
    """A permissive pattern would capture the next word of the sentence."""
    parsed = parse_value("5 parts water", "Solubility")

    assert parsed.value == 5.0
    assert parsed.unit is None


# --------------------------------------------------------------------------
# SI conversion
# --------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("text,heading,si,unit_si", [
    ("135 °C", "Melting Point", 408.15, "K"),
    ("275 °F", "Melting Point", 408.15, "K"),
    ("408.15 K", "Melting Point", 408.15, "K"),
    ("760 mmHg", "Vapor Pressure", 101325.0, "Pa"),
    ("1 atm", "Vapor Pressure", 101325.0, "Pa"),
    ("0.79 g/cm³", "Density", 790.0, "kg/m³"),
    ("1 g/L", "Solubility", 1.0, "kg/m³"),
    ("1000 mg/L", "Solubility", 1.0, "kg/m³"),
    ("1 mg/mL", "Solubility", 1.0, "kg/m³"),
    ("2.47 cP", "Viscosity", 2.47e-3, "Pa·s"),
    ("1.95 mm²/s", "Viscosity", 1.95e-6, "m²/s"),
    ("1 mol/L", "Solubility", 1000.0, "mol/m³"),
])
def test_si_conversion(text, heading, si, unit_si):
    """
    The point of the SI columns: two records in different units compare.

    Aspirin's 135 °C and 275 °F entries both land on 408.15 K, which is the
    cross-check that the temperature conversion is right in both directions.
    """
    parsed = parse_value(text, heading)

    assert parsed.value_si == pytest.approx(si, rel=1e-6), text
    assert parsed.unit_si == unit_si, text


@pytest.mark.unit
def test_range_ends_are_converted_separately():
    """
    A temperature conversion is affine, so a range cannot just be scaled.

    138-140 °C is 411.15-413.15 K; scaling by the °C factor alone would give
    nonsense.
    """
    parsed = parse_value("138-140 °C", "Melting Point")

    assert parsed.value_si is None
    assert parsed.value_min_si == pytest.approx(411.15)
    assert parsed.value_max_si == pytest.approx(413.15)
    assert parsed.unit_si == "K"


@pytest.mark.unit
@pytest.mark.parametrize("unit", ["%", "ppm", "ppb"])
def test_composition_units_are_not_converted(unit):
    """
    A percentage is not a concentration without a density.

    Reporting None here is deliberate: a guess would be indistinguishable from
    a real measurement, and wrong.
    """
    parsed = parse_value(f"0.9 {unit}", "Solubility")

    assert parsed.value == 0.9
    assert parsed.unit == unit
    assert parsed.unit_si is None
    assert parsed.value_si is None


@pytest.mark.unit
def test_every_si_target_is_a_real_si_unit():
    """Guards against a typo in the table silently creating a new "unit"."""
    assert {si for si, _, _ in UNIT_TO_SI.values()} == {
        "K", "Pa", "kg/m³", "mol/m³", "Pa·s", "m²/s", "N/m"}


@pytest.mark.unit
def test_every_alias_target_is_convertible_or_deliberate():
    """
    An alias should point at a spelling the conversion table knows.

    The exceptions are the composition units, which are normalised for
    grouping but deliberately not converted.
    """
    unconvertible = {target for target in UNIT_ALIASES.values()
                     if target not in UNIT_TO_SI}
    assert unconvertible <= {"%", "% w/w", "% v/v", "ppm", "ppb",
                             "g/100g", "hPa", "MPa", "kg/L", "cSt"} | set(UNIT_TO_SI)


# --------------------------------------------------------------------------
# Conditions, operators, qualitative terms
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_measurement_temperature_is_a_condition_not_the_value():
    """
    ``"2.47 cP at 20 °C"`` reports a viscosity, not a temperature.

    Reading the condition's number as the value is the single easiest way to
    get this wrong, so the clause is removed before a value is looked for.
    """
    parsed = parse_value("2.47 cP at 20 °C", "Viscosity")

    assert parsed.value == 2.47
    assert parsed.temperature_c == 20.0


@pytest.mark.unit
def test_melting_point_value_is_not_read_as_a_condition():
    """The converse: 138 °C under "Melting Point" is the value."""
    parsed = parse_value("138 °C", "Melting Point")

    assert parsed.value == 138.0
    assert parsed.temperature_c is None


@pytest.mark.unit
@pytest.mark.parametrize("text,celsius", [
    ("1 g/L at 25 °C", 25.0),
    ("1 g/L @ 25 °C", 25.0),
    ("1 g/L at 77 °F", 25.0),
    ("1 g/L at 298.15 K", 25.0),
])
def test_condition_temperature_is_normalised_to_celsius(text, celsius):
    """So that two records stating the same condition compare equal."""
    assert parse_value(text, "Solubility").temperature_c == pytest.approx(celsius)


@pytest.mark.unit
@pytest.mark.parametrize("text,operator", [
    ("greater than or equal to 100 mg/mL", ">="),
    ("greater than 100 mg/mL", ">"),
    ("less than or equal to 1 mg/mL", "<="),
    ("less than 1 mg/mL", "<"),
    ("≥ 100 mg/mL", ">="),
    ("> 100 mg/mL", ">"),
])
def test_comparison_operators(text, operator):
    """
    A bound is not a measurement, and the difference matters downstream.

    The longest spelling has to win, or "greater than or equal to" reads as
    "greater than".
    """
    parsed = parse_value(text, "Solubility")

    assert parsed.operator == operator, text
    assert parsed.value == pytest.approx(100.0) or parsed.value == pytest.approx(1.0)


@pytest.mark.unit
@pytest.mark.parametrize("text,term", [
    ("Insoluble in water", "insoluble"),
    ("Practically insoluble", "practically insoluble"),
    ("Miscible with ethanol", "miscible"),
    ("Vapor pressure at 20°C: negligible", "negligible"),
])
def test_qualitative_terms(text, term):
    """
    "Insoluble" is a measurement result, so the row is not discarded.

    The longest spelling wins here too: "practically insoluble" must not be
    reported as "soluble".
    """
    assert parse_value(text, "Solubility").qualitative == term



@pytest.mark.unit
@pytest.mark.parametrize("text,term,value", [
    # PubChem's ICSC strings, 2026-09-24 (CIDs 23978, 14917, 2244, 8003,
    # 6344, 6569, 6342).
    ("Solubility in water: none", "none", None),
    ("Solubility in water: very good", "very good", None),
    ("Solubility in water, g/100ml at 15 °C: 0.25 (poor)", "poor", 0.25),
    ("Solubility in water, g/100ml at 20 °C: 0.004 (very poor)", "very poor", 0.004),
    ("Solubility in water, g/100ml at 20 °C: 1.3 (moderate)", "moderate", 1.3),
    ("Solubility in water, g/100ml at 20 °C: 29 (good)", "good", 29.0),
    ("Solubility in water, g/100ml at 20 °C: 1390 (very good)", "very good", 1390.0),
])
def test_icsc_solubility_words(text, term, value):
    """
    The ICSC cards grade water solubility in words. ``"Solubility in water:
    none"`` parsed to nothing at all, and ``"0.25 (poor)"`` lost the word.
    """
    parsed = parse_value(text, "Solubility")

    assert parsed.qualitative == term
    assert parsed.value == value


@pytest.mark.unit
@pytest.mark.parametrize("text", [
    "None reported",
    "Solubility in water: none found in the literature",
    "Good solubility in ethanol",
    "Solubility in ethanol: good",
    "Poor, 0.5 g/L",
])
def test_icsc_words_elsewhere_are_not_qualitative(text):
    """``none``, ``poor`` and ``good`` count only as the whole ICSC value."""
    assert parse_value(text, "Solubility").qualitative is None

@pytest.mark.unit
def test_slash_note_becomes_a_condition():
    """PubChem marks an estimate by wrapping the note in slashes."""
    parsed = parse_value("2.7X10+0 at 25 °C /Estimated/", "Vapor Pressure")

    assert parsed.value == pytest.approx(2.7)
    assert parsed.conditions == "Estimated"


@pytest.mark.unit
def test_pressure_condition_is_recorded():
    """A boiling point is only meaningful with the pressure it was taken at."""
    parsed = parse_value("100 °C at 760 mmHg", "Boiling Point")

    assert parsed.value == 100.0
    assert parsed.conditions is not None
    assert "760" in parsed.conditions


# --------------------------------------------------------------------------
# The heading
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_heading_supplies_celsius_for_a_bare_number():
    """
    A melting point of "135" is 135 °C; depositors often omit the unit.

    Without the heading the number would be unitless and unusable.
    """
    assert parse_value("135", "Melting Point").unit == "°C"
    assert parse_value("135", "Boiling Point").value_si == pytest.approx(408.15)


@pytest.mark.unit
def test_heading_marks_a_dimensionless_quantity():
    """
    A logP of 1.19 is complete as it stands.

    Without the heading, "1.19" under a temperature heading and "1.19" under
    LogP would be indistinguishable, and one of them would get a wrong unit.
    """
    parsed = parse_value("1.19", "LogP")

    assert parsed.value == 1.19
    assert parsed.unit is None
    assert parsed.unit_si is None


@pytest.mark.unit
def test_labelled_value():
    """``"pKa = 3.49 at 25 °C"`` states the value after a label."""
    parsed = parse_value("pKa = 3.49 at 25 °C", "Dissociation Constants")

    assert parsed.value == 3.49
    assert parsed.temperature_c == 25.0
    assert parsed.unit is None


@pytest.mark.unit
def test_unit_stated_in_the_label():
    """
    ``"Vapor pressure, kPa at 20°C: 24"`` puts the unit before the number.

    The value is 24 kPa, and the 20 belongs to the condition — a shape that
    catches a parser reading left to right.
    """
    parsed = parse_value("Vapor pressure, kPa at 20°C: 24", "Vapor Pressure")

    assert parsed.value == 24.0
    assert parsed.unit == "kPa"
    assert parsed.value_si == pytest.approx(24000.0)
    assert parsed.temperature_c == 20.0


@pytest.mark.unit
@pytest.mark.parametrize("text,value,unit,si,celsius", [
    # Benzene's ICSC card: the unit in the label has digits in it.
    ("Solubility in water, g/100ml at 25 °C: 0.18", 0.18, "g/100mL", 1.8, 25.0),
    ("Solubility in water, g/100ml at 20 °C: 0.07 (very poor)", 0.07, "g/100mL", 0.7, 20.0),
    ("Solubility in water, g/100ml: 0.18", 0.18, "g/100mL", 1.8, None),
    ("Solubility in water, mg/100ml at 20 °C: 35", 35.0, "mg/100mL", 0.35, 20.0),
    ("Solubility in water, g/l at 20 °C: 1.8", 1.8, "g/L", 1.8, 20.0),
])
def test_unit_with_digits_stated_in_the_label(text, value, unit, si, celsius):
    """
    The ICSC cards write solubility per 100 mL, with the unit in the label.

    The label's tokens were read without digits, so ``g/100ml`` split into
    ``g/`` and ``ml`` and the value came back with no unit.
    """
    parsed = parse_value(text, "Solubility")

    assert parsed.value == value
    assert parsed.unit == unit
    assert parsed.value_si == pytest.approx(si)
    assert parsed.unit_si == "kg/m³"
    assert parsed.temperature_c == celsius


@pytest.mark.unit
def test_a_formula_in_the_label_is_not_a_unit():
    """Letting digits into label tokens must not make ``C6H6`` a unit."""
    parsed = parse_value("Solubility of C6H6 in water: 1.8", "Solubility")

    assert parsed.value == 1.8
    assert parsed.unit is None


@pytest.mark.unit
def test_heading_is_optional():
    """
    Passing no heading declines the hints rather than failing.

    The parser is usable on a string of unknown provenance.
    """
    parsed = parse_value("0.79 g/cm³")

    assert parsed.value == 0.79
    assert parsed.unit == "g/cm³"
    assert parse_value("135").unit is None


@pytest.mark.unit
def test_parsed_value_is_immutable():
    """It is a record of what a string said, so it does not get edited."""
    parsed = parse_value("135 °C", "Melting Point")

    with pytest.raises(Exception):
        parsed.value = 999
