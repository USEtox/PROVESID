# Cache examples

PROVESID caches every online response on disk, forever, in a per-user cache
directory — `~/.cache/provesid/<service>/` on Linux,
`~/Library/Caches/provesid/<service>/` on macOS,
`%LOCALAPPDATA%\USEtox\provesid\Cache\<service>\` on Windows. Set
`PROVESID_CACHE_DIR` to move it.

| Script | What it shows |
|---|---|
| `cache_layout_and_versioning_demo.py` | Where the cache lives, addressing one service with `service=`, and how a version bump retires entries of a stale shape. Fully offline; uses a throwaway directory. |
| `cache_demo.py` | The cache end to end: unlimited storage, size warnings, export/import, statistics. |

Both run without a network connection:

```bash
python examples/cache/cache_layout_and_versioning_demo.py
python examples/cache/cache_demo.py
```

See `docs/guide/caching.md` for the full reference.
