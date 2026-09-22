"""
Small helpers shared across PROVESID: CAS number checking and the
directories where datasets and cached responses live.

Examples:
    >>> from provesid.utils import check_CASRN
    >>> check_CASRN("50-78-2"), check_CASRN("50-78-3")
    (True, False)
"""

import os

from platformdirs import user_cache_dir, user_data_dir


def _has_casrn_format(s: str):
    return len(s.split("-")) == 3 and all([i.isdigit() for i in s.split("-")])

def check_CASRN(cas_rn: str):
    """
    Check if a string is in the CASRN format and then check if it is a valid CASRN.

    The format is three hyphen-separated runs of digits; the check digit is the
    last, and must equal the sum of the other digits, each weighted by its
    position from the right, modulo 10. The lengths of the runs are not
    checked.

    Args:
        cas_rn: The candidate CAS number.

    Returns:
        (bool): True when the format is right and the check digit agrees.

    Examples:
        >>> check_CASRN("50-78-2")
        True
        >>> check_CASRN("001-16-2")     # a malformed number PubChem lists for aspirin
        False
        >>> check_CASRN("aspirin")
        False
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
    Get the path to the data directory shipped inside the package.

    This holds the small files that ship with PROVESID (the REACH workbook,
    the CAS Common Chemistry Swagger file). The large offline databases live
    under [`user_dataset_path`][provesid.utils.user_dataset_path] instead.

    Returns:
        (str): Absolute path to ``provesid/data``.

    Examples:
        >>> os.path.basename(data_path())
        'data'
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def user_dataset_path(*parts: str, ensure_exists: bool = True) -> str:
    """Return the OS-specific persistent dataset directory for PROVESID.

    The default root comes from `platformdirs` and resolves to a
    per-user data directory that is shared across virtual environments
    on the same machine.

    Power users can override the root directory by setting
    ``PROVESID_DATA_DIR``.

    Args:
        *parts: Optional subdirectories appended to the root directory.
        ensure_exists: When True (default), create the directory.

    Returns:
        Absolute path to the requested dataset directory.

    Examples:
        >>> user_dataset_path("chebifier", ensure_exists=False).endswith("chebifier")
        True
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


def user_cache_path(*parts: str, ensure_exists: bool = True) -> str:
    """Return the OS-specific persistent cache directory for PROVESID.

    This is where [`provesid.cache`][provesid.cache] keeps API responses. It is
    deliberately *not* the system temp directory: most Linux distributions
    clear ``/tmp`` on boot, which silently threw away every cached response
    between sessions even though the caching layer advertises itself as
    persistent. The root comes from `platformdirs` and resolves to a per-user
    cache directory shared across virtual environments on the same machine.

    Cached responses are disposable --- unlike the datasets under
    [`user_dataset_path`][provesid.utils.user_dataset_path], everything here
    can be re-fetched --- which is why the two live under different roots and
    can be cleaned independently.

    Power users can override the root directory by setting
    ``PROVESID_CACHE_DIR``.

    Args:
        *parts: Optional subdirectories appended to the root directory, e.g.
            the service name.
        ensure_exists: When True (default), create the directory.

    Returns:
        Absolute path to the requested cache directory.

    Examples:
        >>> user_cache_path("pubchem", ensure_exists=False).endswith("pubchem")
        True
    """
    override = os.environ.get("PROVESID_CACHE_DIR")
    if override:
        root = os.path.abspath(os.path.expanduser(os.path.expandvars(override)))
    else:
        root = user_cache_dir(appname="provesid", appauthor="USEtox")

    target = os.path.join(root, *parts) if parts else root
    if ensure_exists:
        os.makedirs(target, exist_ok=True)
    return target
