from importlib import import_module

_OPTIONAL_SUBMODULES = (
    'genomics',
    'germplasm',
    'phenotype',
    'gwas',
    'wheat',
)

for _sub in _OPTIONAL_SUBMODULES:
    try:
        import_module(f'{__name__}.{_sub}')
    except Exception:
        # Optional/extra dependencies (e.g. hail) are not required for every runtime.
        pass
