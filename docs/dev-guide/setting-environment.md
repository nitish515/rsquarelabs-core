# Development Env Setup


This is the first time setup doc for the developers.

## First time installation 

```bash
git clone https://github.com/rsquarelabs/framework
virtualenv venv
source venv/bin/activate
pip2.7 install -r requirements.txt
```

Install Gromacs (www.gromacs.org)

In Mac `brew install gromacs`

In ubuntu/debian `sudo apt-get install gromacs`

In Fedora/CentOS `sudo yum install gromacs`


## File system

- **rsquarelabs-core** - main module 
- **rsquarelabs-core/engines** - packages of the modules (gromacs, db_engine(sqlite)) 
- **rsquarelabs-core/utils** - utilities that can be used accross the core package
- **rsquarelabs-core/websuite** webserver that starts the webclient and can be accessed at http://localhost:9090
- **tests** - tests for all the modules/packages
- **sbin** - the python files in this module as supposed to be linked to /usr/local/bin/ so that they can be accessed from anywhere 



## Accessing the command

Once this package is installed via `pip install rsquarelabs-core` , all the python scripts in sbin are supposed to be accessbile from anywhere, 
but to simulate that behaviour in the development, we can create alias for these commands in `.bashrc` or add a softlink, here is how to
add these commands in `.bashrc` ie., `gedit ~/.bashrc`

add the following lines 

```bash
alias r2_gromacs='python /home/<USERNAME>/<R2-CORE-framework-PATH>/sbin/r2_gromacs.py'
alias r2_server_start='python /home/<USERNAME>/<R2-CORE-framework-PATH>/sbin/r2_server_start.py'
```

restart the terminal or do `source ~/.bashrc` to enable these commands for first time





