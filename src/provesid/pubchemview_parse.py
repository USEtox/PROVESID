"""
One parser for the free-text values PubChem's PUG-View returns.

PUG-View reports every experimental property as prose written by whoever
deposited it: ``"138-140 °C"``, ``"8.5X10-5 mm Hg at 25 °C"``,
``"greater than or equal to 100 mg/mL"``, ``"Vapor pressure, kPa at 20°C: 24"``.
A caller that wants to compare two compounds has to turn that back into numbers,
and doing it per call site is how PROVESID ended up with two parsers that
disagreed with each other.

This module holds the only one.
[`parse_value`][provesid.pubchemview_parse.parse_value] takes the string and
the heading it was found under, and returns a
[`ParsedValue`][provesid.pubchemview_parse.ParsedValue] carrying the number (or
the range), the unit as written, the same quantity in SI units, the temperature
the measurement was made at, and any comparison operator or qualitative term
that replaced the number.

Nothing here makes a network request, so it can be tested against real strings
without touching PubChem.

Examples:
    >>> parse_value("138-140 °C", "Melting Point")
    ParsedValue(text='138-140 °C', value=None, value_min=138.0, value_max=140.0, unit='°C', ...)
    >>> round(parse_value("8.5X10-5 mm Hg at 25 °C", "Vapor Pressure").value_si, 7)
    0.0113324
    >>> parse_value("Insoluble in water", "Solubility").qualitative
    'insoluble'
"""

import re
from dataclasses import dataclass
from typing import Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------------

UNIT_ALIASES: Dict[str, str] = {
    # temperature. Single letters are handled by _clean_unit against a
    # whole-token pattern, never by this table: an entry for 'c' would turn the
    # leading letter of 'cP' into a temperature.
    '°c': '°C', 'degc': '°C', 'deg c': '°C', '˚c': '°C',
    '°f': '°F', 'degf': '°F', 'deg f': '°F',
    'kelvin': 'K',
    # pressure
    'mmhg': 'mmHg', 'mm hg': 'mmHg', 'torr': 'mmHg',
    'pa': 'Pa', 'kpa': 'kPa', 'hpa': 'hPa', 'mpa': 'MPa',
    'atm': 'atm', 'bar': 'bar', 'mbar': 'mbar', 'psi': 'psi',
    # density
    'g/cm3': 'g/cm³', 'g/cm³': 'g/cm³', 'g/cu cm': 'g/cm³', 'g/cc': 'g/cm³',
    'g/ml': 'g/mL', 'kg/m3': 'kg/m³', 'kg/m³': 'kg/m³', 'kg/l': 'kg/L',
    # dynamic viscosity
    'cp': 'cP', 'mpa.s': 'mPa·s', 'mpa·s': 'mPa·s', 'mpa s': 'mPa·s',
    'pa.s': 'Pa·s', 'pa·s': 'Pa·s', 'pa s': 'Pa·s',
    # kinematic viscosity
    'cst': 'cSt', 'mm2/s': 'mm²/s', 'mm²/s': 'mm²/s',
    'cm2/s': 'cm²/s', 'cm²/s': 'cm²/s', 'm2/s': 'm²/s',
    # mass concentration
    'g/l': 'g/L', 'mg/l': 'mg/L', 'ug/l': 'µg/L', 'μg/l': 'µg/L', 'µg/l': 'µg/L',
    'ng/l': 'ng/L', 'mg/ml': 'mg/mL', 'g/ml ': 'g/mL',
    'ug/ml': 'µg/mL', 'μg/ml': 'µg/mL', 'µg/ml': 'µg/mL',
    'ng/ml': 'ng/mL', 'mg/100ml': 'mg/100mL', 'g/100ml': 'g/100mL',
    'g/100 ml': 'g/100mL', 'g/100g': 'g/100g',
    # amount concentration
    # 'M' for molar is deliberately absent: it is indistinguishable from
    # metres, from the M of a molecular weight, and from a stray capital in
    # prose, and mapping it cost more wrong answers than it fixed.
    'mol/l': 'mol/L', 'mmol/l': 'mmol/L',
    'µmol/l': 'µmol/L', 'umol/l': 'µmol/L',
    # composition
    '%': '%', '% w/w': '% w/w', '% v/v': '% v/v', 'ppm': 'ppm', 'ppb': 'ppb',
    # surface tension
    'dyn/cm': 'dyn/cm', 'mn/m': 'mN/m', 'n/m': 'N/m',
}
"""Spellings that mean the same unit, mapped to the one this module reports.
PubChem's depositors write the same unit half a dozen ways, and a caller
grouping by unit should not see ``mmHg``, ``mm Hg`` and ``torr`` as three
different things.
"""

UNIT_TO_SI: Dict[str, Tuple[str, float, float]] = {
    # temperature -> kelvin
    '°C': ('K', 1.0, 273.15),
    '°F': ('K', 5.0 / 9.0, 255.3722222222222),
    'K': ('K', 1.0, 0.0),
    # pressure -> pascal
    'Pa': ('Pa', 1.0, 0.0),
    'hPa': ('Pa', 100.0, 0.0),
    'kPa': ('Pa', 1000.0, 0.0),
    'MPa': ('Pa', 1e6, 0.0),
    'mmHg': ('Pa', 133.322387415, 0.0),
    'atm': ('Pa', 101325.0, 0.0),
    'bar': ('Pa', 1e5, 0.0),
    'mbar': ('Pa', 100.0, 0.0),
    'psi': ('Pa', 6894.757293168, 0.0),
    # density and mass concentration -> kg/m³ (1 g/L is exactly 1 kg/m³)
    'g/cm³': ('kg/m³', 1000.0, 0.0),
    'g/mL': ('kg/m³', 1000.0, 0.0),
    'kg/L': ('kg/m³', 1000.0, 0.0),
    'kg/m³': ('kg/m³', 1.0, 0.0),
    'g/L': ('kg/m³', 1.0, 0.0),
    'mg/L': ('kg/m³', 1e-3, 0.0),
    'µg/L': ('kg/m³', 1e-6, 0.0),
    'ng/L': ('kg/m³', 1e-9, 0.0),
    'mg/mL': ('kg/m³', 1.0, 0.0),
    'µg/mL': ('kg/m³', 1e-3, 0.0),
    'ng/mL': ('kg/m³', 1e-6, 0.0),
    'g/100mL': ('kg/m³', 10.0, 0.0),
    'mg/100mL': ('kg/m³', 1e-2, 0.0),
    # amount concentration -> mol/m³
    'mol/L': ('mol/m³', 1000.0, 0.0),
    'mmol/L': ('mol/m³', 1.0, 0.0),
    'µmol/L': ('mol/m³', 1e-3, 0.0),
    # dynamic viscosity -> pascal second
    'cP': ('Pa·s', 1e-3, 0.0),
    'mPa·s': ('Pa·s', 1e-3, 0.0),
    'Pa·s': ('Pa·s', 1.0, 0.0),
    # kinematic viscosity -> m²/s
    'cSt': ('m²/s', 1e-6, 0.0),
    'mm²/s': ('m²/s', 1e-6, 0.0),
    'cm²/s': ('m²/s', 1e-4, 0.0),
    'm²/s': ('m²/s', 1.0, 0.0),
    # surface tension -> newton per metre
    'dyn/cm': ('N/m', 1e-3, 0.0),
    'mN/m': ('N/m', 1e-3, 0.0),
    'N/m': ('N/m', 1.0, 0.0),
}
"""Conversion of a reported unit to the SI unit this module reports alongside
it, as ``(si_unit, scale, offset)`` with ``si = value * scale + offset``.
The offset is what makes a temperature range need its ends converted
separately rather than scaled.

Deliberately absent: ``%``, ``ppm``, ``ppb`` and the ``w/w`` and ``v/v``
variants. They are compositions, not concentrations, and turning one into
kg/m³ needs a density this module does not have. A guess there would be
worse than the honest ``None`` a caller can test for.
"""

CELSIUS_HEADINGS = (
    'melting point', 'boiling point', 'flash point', 'autoignition',
    'decomposition', 'freezing point', 'sublimation',
)
"""Headings whose bare numbers are degrees Celsius. PubChem's depositors often
leave the unit off a melting point, and a melting point without a unit is
not ambiguous in practice.
"""

DIMENSIONLESS_HEADINGS = (
    'logp', 'log p', 'kow', 'dissociation constant', 'ph', 'refractive index',
    'relative density', 'specific gravity', 'logs', 'pka',
)
"""Headings whose values are dimensionless, so that a bare number is a complete
answer rather than a value missing its unit.
"""

QUALITATIVE_TERMS = (
    'practically insoluble', 'very slightly soluble', 'slightly soluble',
    'sparingly soluble', 'freely soluble', 'very soluble', 'readily soluble',
    'miscible', 'immiscible', 'insoluble', 'soluble', 'negligible',
    'decomposes', 'stable', 'not applicable',
)
"""Words that stand in place of a number. PubChem records "insoluble" as often
as it records a solubility, and dropping the row loses the information that
somebody measured it and found it negligible.
"""

_OPERATORS = (
    ('greater than or equal to', '>='),
    ('less than or equal to', '<='),
    ('greater than', '>'),
    ('less than', '<'),
    ('at least', '>='),
    ('>=', '>='), ('<=', '<='), ('≥', '>='), ('≤', '<='),
    ('>', '>'), ('<', '<'), ('ca.', '~'), ('approx.', '~'), ('~', '~'),
)
"""Comparison operators, longest spelling first so that "greater than or equal
to" is not read as "greater than".
"""


# ---------------------------------------------------------------------------
# The parsed value
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ParsedValue:
    """
    A PUG-View value string turned into numbers, with the original kept.

    Attributes:
        text: The string exactly as PubChem returned it. Always populated, so
            nothing this parser fails to understand is ever lost.
        value: The single number the string reports, or None when it reports a
            range, a qualitative term, or nothing numeric at all.
        value_min: Low end of the reported quantity. For a single value this
            equals ``value``, so a numeric filter can use ``value_min`` and
            ``value_max`` without special-casing ranges.
        value_max: High end of the reported quantity.
        unit: The unit as PubChem wrote it, with the spelling normalised
            through [`UNIT_ALIASES`][provesid.pubchemview_parse.UNIT_ALIASES] —
            ``torr`` and ``mm Hg`` both report as ``mmHg``. None when the
            quantity is dimensionless or no unit was found.
        value_si: ``value`` expressed in
            [`unit_si`][provesid.pubchemview_parse.ParsedValue]. None when
            there is no single value, or when the unit has no unambiguous SI
            equivalent.
        value_min_si: ``value_min`` expressed in
            [`unit_si`][provesid.pubchemview_parse.ParsedValue].
        value_max_si: ``value_max`` expressed in
            [`unit_si`][provesid.pubchemview_parse.ParsedValue].
        unit_si: The SI unit the ``*_si`` fields are expressed in: ``K``,
            ``Pa``, ``kg/m³``, ``mol/m³``, ``Pa·s``, ``m²/s`` or ``N/m``. None
            for a dimensionless quantity and for units that cannot be converted
            without information the string does not carry — a percentage or a
            ppm needs a density to become a concentration.
        temperature_c: The temperature the measurement was made at, in degrees
            Celsius, read from an ``at 25 °C`` clause. This is a *condition*,
            not the value: for a melting point of ``138 °C`` the value is 138
            and this field is None.
        operator: ``'>'``, ``'<'``, ``'>='``, ``'<='`` or ``'~'`` when the
            string bounds the quantity rather than stating it, as in
            ``"greater than 100 mg/mL"``. The bound itself is in ``value``.
        qualitative: The term that stood in place of a number, lowercased —
            ``'insoluble'``, ``'miscible'``, ``'negligible'``. Set whether or
            not a number was also found.
        conditions: Any remaining qualifying text, such as a pressure the
            measurement was made at or a ``/Estimated/`` note.

    Examples:
        >>> v = parse_value("2.47 cP at 20 °C", "Viscosity")
        >>> v.value, v.unit, round(v.value_si, 5), v.unit_si, v.temperature_c
        (2.47, 'cP', 0.00247, 'Pa·s', 20.0)
    """

    text: str
    value: Optional[float] = None
    value_min: Optional[float] = None
    value_max: Optional[float] = None
    unit: Optional[str] = None
    value_si: Optional[float] = None
    value_min_si: Optional[float] = None
    value_max_si: Optional[float] = None
    unit_si: Optional[str] = None
    temperature_c: Optional[float] = None
    operator: Optional[str] = None
    qualitative: Optional[str] = None
    conditions: Optional[str] = None

    @property
    def is_numeric(self) -> bool:
        """
        Whether a number was recovered from the string.

        Returns:
            True when either a single value or a range was parsed.

        Examples:
            >>> parse_value("138-140 °C", "Melting Point").is_numeric
            True
            >>> parse_value("Insoluble in water", "Solubility").is_numeric
            False
        """
        return self.value is not None or self.value_min is not None

    @property
    def is_range(self) -> bool:
        """
        Whether the string reported a range rather than a single number.

        Returns:
            True when the two ends differ.

        Examples:
            >>> parse_value("138-140 °C", "Melting Point").is_range
            True
            >>> parse_value("135 °C", "Melting Point").is_range
            False
        """
        return (self.value is None and self.value_min is not None
                and self.value_max is not None)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_NUMBER = r'[-+]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?'
"""A number, including scientific notation. PubChem's ``X10`` notation is
rewritten to ``e`` before this is applied, so it need not be handled here.
"""

_TEMPERATURE_UNIT = r'°\s*[CFK]\b|deg(?:rees?)?\s*[CF]\b|[CFK](?![A-Za-zµμ])'
"""Temperature units, which are the only ones written as a single letter. The
lookahead is what stops the ``c`` of ``cP`` being read as Celsius.
"""

_KNOWN_UNIT_SPELLINGS = sorted(
    set(UNIT_ALIASES) | set(UNIT_ALIASES.values()) | set(UNIT_TO_SI),
    key=len, reverse=True)
"""Every spelling this module recognises, longest first so that ``g/cm³``
wins over ``g`` and ``mm Hg`` over ``mm``. Regex alternation is first-match,
not longest-match, so the ordering is what makes the pattern correct.
"""

_UNIT = (r'(?:' + _TEMPERATURE_UNIT + r'|(?:(?i:'
         + r'|'.join(re.escape(spelling) for spelling in _KNOWN_UNIT_SPELLINGS)
         # A known spelling has to end where the token ends, or the 'pa' of
         # 'parts' is read as pascals.
         + r'))(?![A-Za-zµμ])|[A-Za-zµμ%][A-Za-zµμ%0-9/·²³]{0,11})')
"""A unit token: a known spelling, or — so that an unrecognised unit is
still reported as written rather than dropped — a generic run of
unit-ish characters.
The known spellings match case-insensitively, so ``mm Hg``, ``MM HG`` and
``mmhg`` all resolve; the scoped flag keeps that leniency away from
_TEMPERATURE_UNIT, where a case-insensitive single ``c`` would start reading
stray letters as Celsius.
"""

_RANGE_SEPARATOR = r'\s*(?:to|-|–|—)\s*'

_TEMPERATURE_CLAUSE = re.compile(
    r'(?:at|@)\s*(' + _NUMBER + r')\s*°?\s*([CFK])\b', re.IGNORECASE)
"""``at 25 °C``, ``@ 25 C``."""

_TEMPERATURE_LABEL = re.compile(
    r'\(\s*(' + _NUMBER + r')\s*°?\s*([CFK])\s*\)\s*:', re.IGNORECASE)
"""``(77 °F): 0.3%`` — a parenthesised temperature used as a label for the
value that follows. The trailing colon is what distinguishes it from a
parenthesised value such as ``(135 °C)``, so it is required.
"""

_TEMPERATURE_CONDITIONS = (_TEMPERATURE_CLAUSE, _TEMPERATURE_LABEL)
"""Both spellings of a temperature stated as a condition, tried in order."""

_PRESSURE_CLAUSE = re.compile(
    r'(?:at|@)\s*(' + _NUMBER + r')\s*(mm\s?Hg|torr|kPa|Pa|atm|bar|psi)\b',
    re.IGNORECASE)
"""A pressure stated as a condition rather than as the value."""

_SLASH_NOTE = re.compile(r'/([^/]{2,60})/')
"""``/Estimated/``, ``/from tables/`` — PubChem's own convention for a note."""


def _normalise_text(text: str) -> str:
    """
    Rewrite the notations that would otherwise need their own pattern.

    ``8.5X10-5`` becomes ``8.5e-5``, ``4,600`` becomes ``4600``, non-breaking
    spaces become spaces, and unicode minus signs become ASCII hyphens. Doing
    this once here is what keeps the pattern list below short.

    Args:
        text: Raw value string.

    Returns:
        The string with equivalent notations unified.
    """
    text = text.replace(' ', ' ').replace('−', '-')
    # 8.5X10-5 / 8.5x10+3 -> 8.5e-5 / 8.5e+3
    text = re.sub(r'(\d)\s*[Xx]\s*10\s*([+-]?\d+)', r'\1e\2', text)
    # Thousands separators: "In water, 4,600 mg/L" is one number, not two.
    # Only a comma between a digit and exactly three digits qualifies, so the
    # comma after "water" is left alone.
    while re.search(r'\d,\d{3}(?!\d)', text):
        text = re.sub(r'(\d),(\d{3})(?!\d)', r'\1\2', text)
    # 1.35 X 10-5 with the exponent as a superscript is rare enough to leave.
    return text.strip()


def _clean_unit(raw: Optional[str]) -> Optional[str]:
    """
    Normalise a captured unit token, or reject it.

    Args:
        raw: The text captured where a unit was expected.

    Returns:
        The normalised unit, or None when the token was empty or is one of the
        English words that a permissive pattern captures instead of a unit.
    """
    if not raw:
        return None
    token = re.sub(r'\s+', ' ', raw).strip(' .,;:()[]')
    if not token:
        return None

    # A permissive unit pattern will happily capture the next word of a
    # sentence; these are the ones that actually turn up.
    if token.lower() in {
        'at', 'in', 'on', 'to', 'from', 'with', 'and', 'or', 'of', 'per',
        'approx', 'ca', 'about', 'above', 'below', 'water', 'ethanol', 'air',
        'deg', 'degree', 'degrees', 'wt', 'vol', 'parts', 'times',
    }:
        return None

    return UNIT_ALIASES.get(token.lower(), token)


def _to_si(value: Optional[float], unit: Optional[str]) -> Optional[float]:
    """
    Convert one number to the SI unit for its quantity.

    Args:
        value: The number, or None.
        unit: The normalised unit, or None.

    Returns:
        The converted number, or None when either input is None or the unit has
        no unambiguous SI equivalent.
    """
    if value is None or unit is None:
        return None
    conversion = UNIT_TO_SI.get(unit)
    if conversion is None:
        return None
    _, scale, offset = conversion
    return value * scale + offset


def _heading_unit(heading: Optional[str]) -> Optional[str]:
    """
    The unit a bare number carries by virtue of the heading it appears under.

    Args:
        heading: The PUG-View section heading, e.g. ``"Melting Point"``.

    Returns:
        ``'°C'`` for the temperature headings, None otherwise — including for
        the dimensionless headings, where None is the correct unit rather than
        a missing one.
    """
    if not heading:
        return None
    lowered = heading.lower()
    if any(name in lowered for name in DIMENSIONLESS_HEADINGS):
        return None
    if any(name in lowered for name in CELSIUS_HEADINGS):
        return '°C'
    return None


def _find_qualitative(text: str) -> Optional[str]:
    """
    Find the qualitative term a value string uses in place of a number.

    Args:
        text: Value string.

    Returns:
        The matched term, lowercased, or None. The longest spelling wins, so
        ``"practically insoluble"`` is not reported as ``"soluble"``.
    """
    lowered = text.lower()
    for term in QUALITATIVE_TERMS:
        if re.search(r'\b' + re.escape(term) + r'\b', lowered):
            return term
    return None


def _find_operator(text: str) -> Optional[str]:
    """
    Find a comparison operator bounding the value.

    Args:
        text: Value string.

    Returns:
        One of ``'>'``, ``'<'``, ``'>='``, ``'<='``, ``'~'``, or None.
    """
    lowered = text.lower()
    for spelling, operator in _OPERATORS:
        if spelling in lowered:
            return operator
    return None


def _find_temperature_c(text: str) -> Optional[float]:
    """
    Read the temperature a measurement was made at, in degrees Celsius.

    Args:
        text: Value string.

    Returns:
        The temperature in °C, or None when the string states no condition.
        Fahrenheit and kelvin conditions are converted, so the field is
        comparable across records.
    """
    match = None
    for pattern in _TEMPERATURE_CONDITIONS:
        match = pattern.search(text)
        if match:
            break
    if not match:
        return None
    try:
        number = float(match.group(1))
    except ValueError:
        return None

    scale = match.group(2).upper()
    if scale == 'C':
        return number
    if scale == 'F':
        return (number - 32.0) * 5.0 / 9.0
    return number - 273.15


def _find_conditions(text: str) -> Optional[str]:
    """
    Collect the qualifying text that is neither the value nor the temperature.

    The temperature is reported separately, in
    [`ParsedValue.temperature_c`][provesid.pubchemview_parse.ParsedValue], so
    it is not repeated here; a pressure a measurement was taken at has no field
    of its own and lands here.

    Args:
        text: Value string.

    Returns:
        The conditions joined with ``"; "``, or None when there are none.
    """
    parts = []

    pressure = _PRESSURE_CLAUSE.search(text)
    if pressure:
        parts.append(f"{pressure.group(1)} {pressure.group(2)}")

    parts.extend(note.strip() for note in _SLASH_NOTE.findall(text))

    return '; '.join(parts) if parts else None


def _parse_numbers(text: str, heading: Optional[str]) -> Tuple[
        Optional[float], Optional[float], Optional[float], Optional[str]]:
    """
    Recover the number or range, and the unit, from a value string.

    The patterns are tried in order of how specific they are. A range is looked
    for before a single number, because ``"138-140"`` also matches a single
    number followed by a negative one, and a labelled form is looked for before
    the bare fallback so that ``"Vapor pressure, kPa at 20 °C: 24"`` finds 24
    rather than the 20 of its own condition.

    Args:
        text: Value string, already normalised by `_normalise_text`.
        heading: The section heading, used to decide what a bare number means.

    Returns:
        ``(value, value_min, value_max, unit)``. ``value`` is None for a range
        and ``value_min``/``value_max`` are None for a single value; the caller
        fills in the equal-ends case.
    """
    # The temperature clause is the single biggest source of wrong answers: its
    # number looks exactly like a value. Blank it out before looking for one.
    searchable = text
    for pattern in _TEMPERATURE_CONDITIONS:
        searchable = pattern.sub(' ', searchable)

    # 1. A range, with or without a unit: "138-140 °C", "20 to 25 g/L".
    match = re.match(r'^\s*\(?\s*(' + _NUMBER + r')' + _RANGE_SEPARATOR
                     + r'(' + _NUMBER + r')\s*(' + _UNIT + r')?',
                     searchable)
    if match:
        low, high = float(match.group(1)), float(match.group(2))
        unit = _clean_unit(match.group(3)) or _heading_unit(heading)
        if low > high:
            low, high = high, low
        return None, low, high, unit

    # 2. A number then a unit, optionally bracketed: "0.79 g/cm³", "0.05 [mmHg]".
    match = re.match(r'^\s*\(?\s*(' + _NUMBER + r')\s*\[?\s*(' + _UNIT + r')?\]?',
                     searchable)
    if match:
        value = float(match.group(1))
        unit = _clean_unit(match.group(2)) or _heading_unit(heading)
        return value, None, None, unit

    # 3. A labelled value: "pKa = 3.49", "Vapor pressure, kPa at 20 °C: 24".
    #    The unit, if stated, sits in the label rather than after the number.
    match = re.search(r'([A-Za-z][A-Za-z0-9\s,()/·°%µμ]{0,40}?)\s*[=:]\s*('
                      + _NUMBER + r')\s*(' + _UNIT + r')?', searchable)
    if match:
        value = float(match.group(2))
        unit = _clean_unit(match.group(3))
        if unit is None:
            # "Vapor pressure, kPa at 20°C: 24" -- the unit is in the label.
            # Digits may follow the first character, for the ICSC cards'
            # "Solubility in water, g/100ml at 25 °C: 0.18"; a formula such
            # as C6H6 is then a token too, but not a known unit.
            label_units = [_clean_unit(token) for token
                           in re.findall(r'[A-Za-zµμ°%][A-Za-z0-9/·°%²³]*',
                                         match.group(1))]
            known = [candidate for candidate in label_units
                     if candidate in UNIT_TO_SI]
            unit = known[-1] if known else None
        return value, None, None, unit or _heading_unit(heading)

    # 4. Any number at all, with whatever follows it.
    match = re.search(r'(' + _NUMBER + r')\s*(' + _UNIT + r')?', searchable)
    if match:
        value = float(match.group(1))
        unit = _clean_unit(match.group(2)) or _heading_unit(heading)
        return value, None, None, unit

    return None, None, None, None


def parse_value(text: Optional[str], heading: Optional[str] = None) -> ParsedValue:
    """
    Parse one PUG-View value string into numbers, units and conditions.

    This is the only parser in the package for these strings. Every PubChemView
    method that reports a value routes through it, so a fix here reaches all of
    them.

    Args:
        text: The value string, as PubChem returned it. None or empty gives a
            [`ParsedValue`][provesid.pubchemview_parse.ParsedValue] carrying
            nothing but the (empty) text.
        heading: The PUG-View heading the value was found under, e.g.
            ``"Melting Point"``. Used for two decisions a string alone cannot
            settle: a bare number under a temperature heading is degrees
            Celsius, and a bare number under a dimensionless heading such as
            ``"LogP"`` or ``"Dissociation Constants"`` is complete as it
            stands. Passing None is safe and simply declines both hints.

    Returns:
        A [`ParsedValue`][provesid.pubchemview_parse.ParsedValue]. It always
        carries ``text``; every other field is None when the string did not
        state it. Nothing is invented: an unrecognised unit is reported as
        written with no SI conversion, and a string with no number gives
        ``value is None`` rather than a zero.

    Examples:
        >>> parse_value("138-140 °C", "Melting Point").value_max_si
        413.15
        >>> parse_value("greater than or equal to 100 mg/mL", "Solubility").operator
        '>='
        >>> parse_value("2.7e+0 at 25 °C /Estimated/", "Vapor Pressure").conditions
        'Estimated'
    """
    if not text:
        return ParsedValue(text=text or '')

    normalised = _normalise_text(text)

    qualitative = _find_qualitative(normalised)
    operator = _find_operator(normalised)
    temperature_c = _find_temperature_c(normalised)
    conditions = _find_conditions(normalised)

    value, value_min, value_max, unit = _parse_numbers(normalised, heading)

    # A single value is also a degenerate range, so a numeric filter can use the
    # bounds without asking which shape it got.
    if value is not None:
        value_min = value_max = value

    unit_si = UNIT_TO_SI[unit][0] if unit in UNIT_TO_SI else None

    return ParsedValue(
        text=text,
        value=value,
        value_min=value_min,
        value_max=value_max,
        unit=unit,
        value_si=_to_si(value, unit),
        value_min_si=_to_si(value_min, unit),
        value_max_si=_to_si(value_max, unit),
        unit_si=unit_si,
        temperature_c=temperature_c,
        operator=operator,
        qualitative=qualitative,
        conditions=conditions,
    )
