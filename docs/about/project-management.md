# Project management

This document outlines our project management processes for R2-Core Framework.


The aim is to ensure that the project has a high ["bus factor"](https://en.wikipedia.org/wiki/Bus_factor), and can continue to remain well
 supported for the foreseeable future. Suggestions for improvements to our process are welcome.


## Maintenance team

We have a quarterly maintenance cycle where new members may join the maintenance team. 
We currently cap the size of the team at 5 members, and may encourage folks to step out of the 
team for a cycle to allow new members to participate.


### Current Team

The current maintenance team includes:

1. [@rrmerugu](https://github.com/rrmerugu)
2. [@nitish515](https://github.com/nitish515)


### Maintenance cycles

Each maintenance cycle is initiated by an issue being opened with the `Process` label.

To be considered for a maintainer role simply comment against the issue.
Existing members must explicitly opt-in to the next cycle by check-marking their name.
The final decision on the incoming team will be made by `@rrmerugu`

Members of the maintenance team will be added as collaborators to the repository.


###  Responsibilities of team members.

Team members have the following responsibilities.

1. Close invalid or resolved tickets.
2. Add triage labels and milestones to tickets.
3. Merge finalized pull requests.
4. Build and deploy the documentation, using `mkdocs gh-deploy`.
5. Build and update the included translation packs.


Further notes for maintainers:

1. Code changes should come in the form of a pull request - do not push directly to master.
2. Maintainers should typically not merge their own pull requests.
3. Each issue/pull request should have exactly one label once triaged.
4. Search for un-triaged issues with [is:open no:label](https://github.com/rsquarelabs/rsquarelabs-core/issues?q=is%3Aopen+no%3Alabel).


It should be noted that participating actively in this project clearly 
** does not require being part of the maintenance team**. Almost every import part of issue triage and project 
improvement can be actively worked on regardless of your collaborator status on the repository.



## Project requirements

All our test requirements are pinned to exact versions, in order to ensure that our test runs are reproducible. 
We maintain the requirements in the `requirements` directory. The requirements files are referenced from the `tox.ini`
configuration file, ensuring we have a single source of truth for package versions used in testing.

Package upgrades should generally be treated as isolated pull requests. You can check if there are any packages available at 
a newer version, by using the `pip list --outdated`.


## Project ownership

The PyPI package is owned by [@rrmerugu](https://github.com/rrmerugu). As a backup [@shivamsk](https://github.com/shivamsk) also has ownership of the package.

If [@rrmerugu](https://github.com/rrmerugu) ceases to participate in the project then [@shivamsk](https://github.com/shivamsk) has responsibility for handing over ownership duties.


