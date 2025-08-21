
def _has_casrn_format(s: str):
    return len(s.split("-")) == 3 and all([i.isdigit() for i in s.split("-")])

def check_CASRN(cas_rn: str):
    """
    Check if a string is in the CASRN format and then check if it is a valid CASRN
    """
    # Check if the CASRN has the correct format
    if not _has_casrn_format(cas_rn):
        return False

    # Split the CASRN into its parts
    parts = cas_rn.split("-")
    if len(parts) != 3:
        return False

    # Extract the digits and the check digit
    digits = "".join(parts[:-1])
    check_digit = int(parts[-1])

    # Calculate the check digit
    calculated_check_digit = 0
    for i, digit in enumerate(reversed(digits)):
        calculated_check_digit += (i + 1) * int(digit)

    # Validate the check digit
    return calculated_check_digit % 10 == check_digit