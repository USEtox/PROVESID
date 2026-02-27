---
name: dev-principles
description: PROVESID development principles and conventions. Load automatically when writing code, adding features, implementing APIs, handling data retrieval, writing tests, or making any architectural decision for this project.
user-invocable: false
---

# PROVESID Development Principles

Apply these rules consistently across all contributions to the PROVESID package.

---

## 1. Active Development — No Backward Compatibility

The codebase is under **heavy active development**. Breaking changes between versions are expected and acceptable. Do **not** add deprecation shims or compatibility layers unless explicitly requested. Focus on correctness and clean design over preserving old interfaces.

---

## 2. Offline Databases Are the Priority

Serve data from **local/offline databases** (SQLite, cached flat files, embedded datasets) wherever possible. Online API calls are a secondary fallback. When adding new data retrieval logic, always check whether an offline source can satisfy the request first.

---

## 3. Online API Calls: Preserve Core, Extend With Wrappers

- **Do not alter** core API call logic, endpoint structure, or raw response format. Treat the upstream API contract as immutable.
- **Add new methods** on top of existing ones to provide cleaner interfaces, formatted outputs, and user-friendly abstractions.
- Clearly name and document wrapper methods as convenience layers over the raw API methods.

---

## 4. Readability and Elegance Over Micro-Optimizations

Prefer clear variable names, well-structured functions, and idiomatic Python. Only deviate when there is a **significant, measurable computational benefit** (e.g., orders-of-magnitude speedup on large datasets). Document any such trade-offs with inline comments.

---

## 5. Tests Run Locally Only — No GitHub Actions CI

Tests are executed **locally** with `pytest`. Do **not** configure, modify, or add GitHub Actions workflows for running tests. Do not add or modify `.github/workflows/` CI/CD test jobs.

---

## 6. Mandatory Detailed Docstrings

Every public method, function, and class **must** have a comprehensive docstring. Use Google-style or NumPy-style consistently within each module. Cover:

- One-line summary
- Extended description for non-trivial behavior
- `Args:` / `Parameters:` — all parameters with types and descriptions
- `Returns:` — return type and description
- `Raises:` — any exceptions that may be raised
- Example usage where helpful

Private methods (prefixed with `_`) should also be documented when their logic is non-obvious.

---

## 7. Example Scripts and Jupyter Notebooks for All Features

- **Standalone scripts** go in `examples/` at the project root.
- **Jupyter notebooks** go in `examples/<feature-area>/` subfolders (e.g., `examples/pubchem/`, `examples/chebi/`), one subfolder per integration or feature.
- Examples must be self-contained, clearly commented, and demonstrate real usage scenarios.
- When adding a new feature, add or update the corresponding example file(s).

---

## 8. Graceful API Rate-Limit and Error Handling

- Detect **HTTP 429** and other rate-limit responses; apply exponential back-off with a configurable retry count.
- Catch network errors, timeouts, and unexpected status codes; raise informative, package-specific exceptions — never let raw `requests` or `httpx` exceptions propagate.
- Use the standard `logging` module for warnings and retries. Never use `print()` for operational messages.
- Centralize rate-limit handling (shared utility or base class) rather than duplicating it across modules.

---

## 9. Offline → Online Fallback

Data retrieval must follow a **two-stage lookup**:

1. **Primary:** Query the local/offline database or cache.
2. **Fallback:** If not found offline, transparently fall back to the corresponding online API, optionally cache the result, and return it to the caller.

Fallback behavior must be:
- Opt-in or opt-out via a parameter (e.g., `use_online_fallback: bool = True`).
- Clearly documented in the method's docstring.
- Logged at `DEBUG` level to allow observability without production noise.
