# framework

[![Build Status](https://travis-ci.org/rsquarelabs/framework.svg?branch=dev)](https://travis-ci.org/rsquarelabs/framework)
[![Requirements Status](https://requires.io/github/rsquarelabs/framework/requirements.svg?branch=dev)](https://requires.io/github/rsquarelabs/framework/requirements/?branch=dev)

![framework hero ](docs/images/rsquarelabs-hero.jpg)

This is the library of automation pipeline modules

[Website](http://rsquarelabs.org) |
[Documentation](https://github.com/rsquarelabs/framework/wiki/) |
[Installation](https://github.com/rsquarelabs/framework/wiki#installation) |
[Mailing List](https://groups.google.com/d/forum/framework) |
[Gitter Chat](https://gitter.im/rsquarelabs/framework)


## Summary
- [**Install**](#install)
- [**Features**](#features)
- [**Usage**](#usage)
- [**Why framework?**](#why-framework)
- [**Community**](#community)
- [**Roadmap**](#roadmap)
- [**Support**](#support)
- [**License**](#license)


## Install
```
pip install framework
```
**we currently support python 2.7 only**

## Usage
```
# first: import the class that does protein-ligand minimisation
from rsquarelabs_core.engines.gromacs import ProteinLigMin

# Create and object with input files 
obj = ProteinLigMin(
    ligand_file='ligand.gro',
    ligand_topology_file='ligand.itp',
    protein_file='protein.pdb',
    working_dir='./',
    verbose=True,
    quiet=False
)

# call the method you want to start with
obj.create_topology()
obj.prepare_system()
obj.write_em_mdp()
obj.add_ions()
obj.write_emreal_mdp()
obj.minimize()

 

```

## Features
1. Start a light weight webserver
2. commands for /usr/local/bin/

## Why framework
1. scaffolding the project
2. project management
3. Tracking the project via webclient


## Community
Want to join an open source project? Now it's your chance!

Don't know what you want to help out with? Well here are some areas that we could use help with:

- framework scientific thoughts,
- framework [technical stack](_docs/notes/technical-stack.md)
- framework code - see [milestones](https://github.com/rsquarelabs/framework/milestones) and [issues](https://github.com/rsquarelabs/framework/issues) for more
- framework [documentation](https://github.com/rsquarelabs/framework/wiki), [tutorials](https://github.com/rsquarelabs/framework/wiki/Tutorials) and [examples](https://github.com/rsquarelabs/framework/wiki/Examples)
- framework [website](http://rsquarelabs.org)
- If you want to join our community as a contributor, please leave a message as [@rrmerugu](https://twitter.com/rrmerugu)



## RoadMap
You can find a detailed Roadmap of framework [here](https://github.com/rsquarelabs/framework/milestones).

## Support
We support universties and research labs to install, configure and setup rl-core

## License
