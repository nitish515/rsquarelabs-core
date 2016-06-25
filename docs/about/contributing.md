# Contributing


## Development
If you want to **use the demo** app to work on this package:

In the project [repository](https://github.com/rsquarelabs/framework) you can find the demo app(at /demo). It is a project with Django & Django Rest Framework that will allow you to work with this project.

From the root of the repository:

```bash
# Create the virtual environment
virtualenv venv

# Install requirements
venv/bin/pip install -r requirements/dev-requirements.txt

# Activate the environment
source venv/bin/activate

# Access the commands
python sbin/r2_gromacs.py init # for gromacs module
python sbin/r2_server_start # start the webclient in localhost
```


**Note:** You do not need a database or to run migrate.

 
## Contributing to the project

1. Fork the `rsquarelabs/framework` repository!
2. Clone the repository to your local machine.
3. Select an issue from [rsquarelabs/framework](https://github.com/rsquarelabs/framework/issues?q=is%3Aopen) to work on or submit a proposal of your own.

2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request
6. Make sure tests are passing


## Pull Request 

It's a good idea to make pull requests early on. A pull request represents the start of a discussion, 
and doesn't necessarily need to be the final, finished submission.

It's also always best to make a new branch before starting work on a pull request. 
This means that you'll be able to later switch back to working on another separate issue without interfering with an ongoing pull requests.

It's also useful to remember that if you have an outstanding pull request then pushing new commits to your 
GitHub repo will also automatically update the pull requests.

GitHub's documentation for working on pull requests is [available here](https://help.github.com/articles/using-pull-requests/).

Always run the tests before submitting pull requests, and ideally run tox in order to check that your 
`modifications are compatible` with both Python 2 and Python 3.

Once you've made a pull request take a look at the `Travis build status` in the GitHub interface and make sure 
the tests are running as you'd expect.

## Code of Conduct

As contributors and maintainers of this project, and in the interest of
fostering an open and welcoming community, we pledge to respect all people who
contribute through reporting issues, posting feature requests, updating
documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free
experience for everyone, regardless of level of experience, gender, gender
identity and expression, sexual orientation, disability, personal appearance,
body size, race, ethnicity, age, religion, or nationality.

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery
* Personal attacks
* Trolling or insulting/derogatory comments
* Public or private harassment
* Publishing other's private information, such as physical or electronic
  addresses, without explicit permission
* Other unethical or unprofessional conduct

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

By adopting this Code of Conduct, project maintainers commit themselves to
fairly and consistently applying these principles to every aspect of managing
this project. Project maintainers who do not follow or enforce the Code of
Conduct may be permanently removed from the project team.

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting a project maintainer at ravi@rsquarelabs.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. Maintainers are
obligated to maintain confidentiality with regard to the reporter of an
incident.


This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.3.0, available at
[http://contributor-covenant.org/version/1/3/0/][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/3/0/