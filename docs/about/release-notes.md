# Release Management



## Versioning

Minor version numbers (0.0.x) are used for changes that are API compatible. You should be able to upgrade between minor point releases without any other code changes.

Medium version numbers (0.x.0) may include API changes, in line with the deprecation policy. You should read the release notes carefully before upgrading between medium point releases.

Major version numbers (x.0.0) are reserved for substantial project milestones.


##  Release process

The release manager is selected on every quarterly maintenance cycle.

1. The manager should be selected by [@rrmerugu](https://github.com/rrmerugu).
2. The manager will then have the maintainer role added to PyPI package.
3. The previous manager will then have the maintainer role removed from the PyPI package.

Our PyPI releases will be handled by either the current release manager, or by [@rrmerugu](https://github.com/rrmerugu). 
Every release should have an open issue tagged with the Release label and marked against the appropriate milestone.

The following template should be used for the description of the issue, and serves as a release checklist.



```
Release manager is @***.
Pull request is #***.

During development cycle:

- [ ] Upload the new content to be translated to [transifex](https://github.com/rsquarelabs/rsquarelabs-core/about/project-management/#translations).


Checklist:

- [ ] Create pull request for [release notes](https://github.com/rsquarelabs/rsquarelabs-core/blob/master/docs/topics/release-notes.md) based on the [*.*.* milestone](https://github.com/rsquarelabs/rsquarelabs-core/milestones/***).
- [ ] Ensure the pull request increments the version to `*.*.*` in [`rsquarelabs_core/__init__.py`](https://github.com/rsquarelabs/rsquarelabs-core/blob/master/rsquarelabs_core/__init__.py).
- [ ] Confirm with @rrmerugu that release is finalized and ready to go.
- [ ] Ensure that release date is included in pull request.
- [ ] Merge the release pull request.
- [ ] Push the package to PyPI with `./setup.py publish`.
- [ ] Tag the release, with `git tag -a *.*.* -m 'version *.*.*'; git push --tags`.
- [ ] Deploy the documentation with `mkdocs gh-deploy`.
- [ ] Make a release announcement on the [discussion group](https://groups.google.com/forum/?fromgroups#!forum/django-rest-framework).
- [ ] Make a release announcement on twitter.
- [ ] Close the milestone on GitHub.

To modify this process for future releases make a pull request to the [project management](http://rsquarelabs.github.io/rsquarelabs-core/about/project-management/) documentation.

```