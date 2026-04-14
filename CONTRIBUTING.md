# Contributing to AGFusion

## Releasing a New Version

AGFusion uses [bump2version](https://github.com/c4urself/bump2version) to manage versioning. It automatically updates the version in `setup.py`, commits the change, and creates a git tag. Pushing a tag triggers the PyPI deployment workflow automatically.

**Install bump2version:**

```
pip install bump2version
```

### Stable releases

```
# Patch release (e.g. 1.4.3 -> 1.4.4)
bump2version patch

# Minor release (e.g. 1.4.3 -> 1.5.0)
bump2version minor

# Major release (e.g. 1.4.3 -> 2.0.0)
bump2version major
```

### Pre-releases (alpha / beta / release candidate)

Use `--new-version` to set the exact pre-release version following [PEP 440](https://peps.python.org/pep-0440/) conventions (`a` = alpha, `b` = beta, `rc` = release candidate):

```
# Start an alpha (e.g. next minor alpha)
bump2version --new-version 1.6.0a1 minor

# Increment alpha number
bump2version --new-version 1.6.0a2 minor

# Promote to beta
bump2version --new-version 1.6.0b1 minor

# Promote to release candidate
bump2version --new-version 1.6.0rc1 minor

# Final stable release
bump2version --new-version 1.6.0 minor
```

Pre-releases are published to PyPI and are installable with:

```
pip install --pre agfusion
```

### Push to trigger deployment

After bumping, push the commit and tag to GitHub — the CI workflow will publish to PyPI automatically:

```
git push && git push --tags
```
