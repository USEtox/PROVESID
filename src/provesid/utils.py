
import os

from platformdirs import user_data_dir


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

def data_path():
    """
    Get the path to the data directory
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def user_dataset_path(*parts: str, ensure_exists: bool = True) -> str:
    """Return the OS-specific persistent dataset directory for PROVESID.

    The default root comes from :mod:`platformdirs` and resolves to a
    per-user data directory that is shared across virtual environments
    on the same machine.

    Power users can override the root directory by setting
    ``PROVESID_DATA_DIR``.

    Args:
        *parts: Optional subdirectories appended to the root directory.
        ensure_exists: When True (default), create the directory.

    Returns:
        Absolute path to the requested dataset directory.
    """
    override = os.environ.get("PROVESID_DATA_DIR")
    if override:
        root = os.path.abspath(os.path.expanduser(os.path.expandvars(override)))
    else:
        root = user_data_dir(appname="provesid", appauthor="USEtox")

    target = os.path.join(root, *parts) if parts else root
    if ensure_exists:
        os.makedirs(target, exist_ok=True)
    return target